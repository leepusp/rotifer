"""
Showing how a long SQL operation is getting on.

Loading a release takes minutes to hours, and for most of that time a
cursor has nothing to say. These helpers give it something: a bar that
counts rows as they are inserted, and, for the loads that hand the
file to another program, a watcher that asks the database how far it
has got.

Everything here is inert when ``progress`` is off, so a caller need
not decide twice: the code reads the same either way.

See Also
--------
rotifer.db.uniprot.clickhouse.MappingCursor.load : counts rows on the server
rotifer.db.uniprot.sqlite3.MappingCursor.load : counts rows, or bytes written
"""

import os
import threading
import types

import rotifer
logger = rotifer.logging.getLogger(__name__)


def estimate_rows(path, sample=1 << 23, compressed=False):
    """
    Guess how many rows a delimited file holds, without reading it.

    A bar needs a total to show a proportion, and counting the lines
    of a 90 GB file to find one would cost more than the bar is worth.
    The first few megabytes give an average line length, and the size
    of the file does the rest.

    Parameters
    ----------
    path : str
        The file to look at.
    sample : int, optional
        How many bytes to read. The default is 8 MiB, which is plenty
        for an average and quick even on a cold cache.
    compressed : bool, default False
        Whether the file is gzip compressed, in which case its size on
        disk says nothing about how many rows it holds and None is
        returned.

    Returns
    -------
    int or None
        The estimate, or None when there is no honest one to give.
    """
    if compressed:
        return None
    try:
        size = os.path.getsize(path)
        if not size:
            return 0
        with open(path, 'rb') as fh:
            head = fh.read(min(sample, size))
    except OSError:
        logger.debug(f'Could not sample {path}', exc_info=1)
        return None
    lines = head.count(b"\n")
    if not lines:
        return None
    return int(size / (len(head) / lines))


class Progress:
    """
    A progress bar that costs nothing when it is switched off.

    Parameters
    ----------
    total : int, optional
        What counts as finished. Without one the bar shows a running
        count and a rate, which is all there is to say when the end is
        unknown.
    unit : str, default 'rows'
        What is being counted.
    desc : str, optional
        Label shown beside the bar.
    enabled : bool, default True
        Whether to show anything at all.
    position : int, optional
        Which line to draw on. A delegator takes 0 and its backends 1,
        so the overall count stays above the backend working on it
        rather than the two overwriting each other.
    leave : bool, default True
        Whether the bar stays on screen once finished. A backend's
        does not: several of them run in turn for one query, and only
        the total is worth keeping.
    """

    def __init__(self, total=None, unit='rows', desc=None, enabled=True,
                 position=None, leave=True):
        self.total = total
        self._bar = None
        if not enabled:
            return
        try:
            from tqdm import tqdm
            self._bar = tqdm(total=total, unit=unit, unit_scale=True,
                             desc=desc, initial=0, position=position, leave=leave)
        except Exception:
            logger.debug('Could not open a progress bar', exc_info=1)

    def update(self, amount):
        """
        Move the bar on.

        Parameters
        ----------
        amount : int
            How much further along.
        """
        if self._bar is not None and amount:
            self._bar.update(amount)

    def position(self, value):
        """
        Move the bar to an absolute point.

        For a watcher that asks how far something else has got, rather
        than being told how much it just did.

        Parameters
        ----------
        value : int
            Where the bar should now be.
        """
        if self._bar is None:
            return
        self.update(max(0, value - self._bar.n))

    def close(self):
        """Take the bar down."""
        if self._bar is not None:
            self._bar.close()
            self._bar = None

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()
        return False


class Watcher:
    """
    Report on work another process is doing.

    A load that hands the file to another program has nothing to count
    as it goes, so instead something is asked periodically how far it
    has got: the number of rows in the table, or the size of the file
    being written.

    The question is asked on a thread and its failures are ignored: a
    bar that cannot be drawn must never be the reason a load stops.

    Parameters
    ----------
    poll : callable
        Called with no arguments, returns how far along things are.
    total : int, optional
        What counts as finished.
    interval : float, default 2.0
        Seconds between questions.
    unit : str, default 'rows'
        What is being counted.
    desc : str, optional
        Label shown beside the bar.
    enabled : bool, default True
        Whether to watch at all.
    """

    def __init__(self, poll, total=None, interval=2.0, unit='rows',
                 desc=None, enabled=True):
        self._poll = poll
        self._interval = interval
        self._stop = threading.Event()
        self._thread = None
        self._progress = Progress(total=total, unit=unit, desc=desc, enabled=enabled)
        self._enabled = enabled and self._progress._bar is not None

    def _run(self):
        while not self._stop.wait(self._interval):
            try:
                self._progress.position(int(self._poll() or 0))
            except Exception:
                logger.debug('Progress poll failed', exc_info=1)

    def __enter__(self):
        if self._enabled:
            self._thread = threading.Thread(target=self._run, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, *exception):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self._interval + 1)
        # One last look, so the bar ends where the work did rather
        # than wherever the last poll happened to catch it.
        if self._enabled and not any(exception):
            try:
                self._progress.position(int(self._poll() or 0))
            except Exception:
                logger.debug('Final progress poll failed', exc_info=1)
        self._progress.close()
        return False
