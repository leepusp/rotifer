#!/usr/bin/env python3
"""Tests for the progress helpers the SQL backends use.

A bar is decoration, so the only thing that must be true of it is that it
never becomes the reason a load fails. That is what most of these check: the
poll that raises, the file that vanishes, the terminal that cannot draw. The
rest check the estimate, which is the one part that could be quietly wrong --
a bar showing a confident and false proportion is worse than one showing none.

Run under pytest, or standalone: python test/db/test_sql_progress.py
"""

import os
import time

import pytest

from rotifer.db.sql import progress as sqlprog


def write_lines(path, n, width=40):
    with open(path, 'wt') as fh:
        for i in range(n):
            fh.write(('x' * width + str(i))[:width] + "\n")
    return path


# ------------------------------------------------------- the estimate

def test_the_estimate_is_close_enough_to_be_worth_showing(tmp_path):
    """It is read off a sample rather than counted, since counting the lines
    of a 90 GB file costs more than the bar is worth."""
    path = write_lines(str(tmp_path / 'rows.tsv'), 50000)
    estimate = sqlprog.estimate_rows(path)
    assert abs(estimate - 50000) / 50000 < 0.01


def test_a_short_file_is_estimated_exactly(tmp_path):
    """Shorter than the sample means the whole file was read, so there is
    nothing to extrapolate and no room to be wrong."""
    path = write_lines(str(tmp_path / 'rows.tsv'), 10)
    assert sqlprog.estimate_rows(path) == 10


def test_an_empty_file_is_no_rows(tmp_path):
    path = str(tmp_path / 'empty.tsv')
    open(path, 'w').close()
    assert sqlprog.estimate_rows(path) == 0


def test_a_compressed_file_gets_no_estimate(tmp_path):
    """Its size on disk says nothing about how many rows it holds, so the
    honest answer is none: a bar with no total still shows a rate."""
    path = write_lines(str(tmp_path / 'rows.tsv'), 100)
    assert sqlprog.estimate_rows(path, compressed=True) is None


def test_a_missing_file_gets_no_estimate(tmp_path):
    assert sqlprog.estimate_rows(str(tmp_path / 'nope.tsv')) is None


def test_a_file_of_one_long_line_gets_no_estimate(tmp_path):
    """No newline in the sample means no average line length, and dividing by
    zero is not an estimate."""
    path = str(tmp_path / 'oneline.tsv')
    with open(path, 'wt') as fh:
        fh.write('x' * 4096)
    assert sqlprog.estimate_rows(path) is None


# ---------------------------------------------------------- the bar

def test_a_disabled_bar_does_nothing_and_says_nothing(capsys):
    with sqlprog.Progress(total=10, enabled=False) as bar:
        bar.update(5)
        bar.position(9)
    assert capsys.readouterr().err == ''


def test_updates_and_positions_agree(tmp_path):
    """update() says how much was just done, position() says where things now
    stand: a watcher asking another process knows the second, not the first."""
    bar = sqlprog.Progress(total=100, enabled=True)
    bar.update(10)
    bar.position(30)
    assert bar._bar.n == 30
    bar.position(20)          # never goes backwards
    assert bar._bar.n == 30
    bar.close()


def test_closing_twice_is_harmless():
    bar = sqlprog.Progress(total=1, enabled=True)
    bar.close()
    bar.close()


# ------------------------------------------------------- the watcher

def test_the_watcher_follows_what_it_is_told(tmp_path):
    state = {'n': 0}
    with sqlprog.Watcher(lambda: state['n'], total=100, interval=0.05,
                         enabled=True) as watcher:
        for value in (10, 40, 90):
            state['n'] = value
            time.sleep(0.12)
    assert watcher._progress._bar is None      # closed on the way out


def test_a_poll_that_raises_does_not_stop_the_work():
    """The bar is decoration. Whatever it is watching must finish whether or
    not the question can be answered."""
    def broken():
        raise RuntimeError('no answer')

    done = []
    with sqlprog.Watcher(broken, total=10, interval=0.05, enabled=True):
        time.sleep(0.15)
        done.append(True)
    assert done == [True]


def test_a_disabled_watcher_never_polls():
    """Nothing should be asked when nobody is watching: the question may cost
    a query."""
    asked = []
    with sqlprog.Watcher(lambda: asked.append(1) or 0, interval=0.01,
                         enabled=False):
        time.sleep(0.1)
    assert asked == []


def test_the_watcher_finishes_where_the_work_did():
    """The last poll happens on the way out, so the bar ends where things
    actually got to rather than wherever the interval last caught them."""
    state = {'n': 0}
    with sqlprog.Watcher(lambda: state['n'], total=100, interval=5,
                         enabled=True) as watcher:
        state['n'] = 100
    assert watcher._progress.total == 100


def test_an_exception_inside_still_closes_the_bar():
    with pytest.raises(ValueError):
        with sqlprog.Watcher(lambda: 0, total=10, interval=0.05, enabled=True):
            raise ValueError('the work failed')


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))


# --------------------------------------------- who draws the bar, and where

def test_the_bars_stack_rather_than_overwrite():
    """A delegator and the backend it is driving both draw. Without a line
    each they would write over one another, so the total sits above the
    backend working towards it."""
    outer = sqlprog.Progress(total=10, enabled=True, position=0)
    inner = sqlprog.Progress(total=10, enabled=True, position=1, leave=False)
    assert outer._bar.pos == 0
    assert abs(inner._bar.pos) == 1        # tqdm keeps it negated internally
    inner.close(); outer.close()


def test_a_backend_bar_does_not_stay_behind():
    """Several backends run in turn for one query; only the total is worth
    keeping once they are done."""
    inner = sqlprog.Progress(total=1, enabled=True, position=1, leave=False)
    assert inner._bar.leave is False
    inner.close()


def test_dictionary_access_is_answered_quietly(capsys):
    """A lookup is not a job to watch. fetchall and fetchone report; mc[x]
    does not."""
    from rotifer.db import uniprot

    class Quiet(uniprot.BaseUniProtDelegatorCursor):
        def __init__(self):
            self._backends = {}
            super().__init__(readers=[])

        def reset_cursors(self):
            self.cursors = {}

    cursor = Quiet()
    cursor.progress = True
    cursor['P00750']
    assert capsys.readouterr().err == ''
    # and the setting is put back, since it is shared with the backends
    assert cursor.progress is True


def test_a_sql_backend_draws_for_dictionary_access_too(tmp_path):
    """The exception to the rule above. A SQL backend does its batching in
    __getitem__, and a query there can be as long as any other, so it reports
    like the rest. Only the delegator's dictionary access stays quiet, and it
    silences its backends by clearing the setting they share."""
    from rotifer.db.uniprot import sqlite3 as rus
    cursor = rus.MappingCursor(str(tmp_path / 'q.sqlite3'), progress=True)
    cursor.create()
    assert cursor.progress is True
    cursor['P00750']                  # draws, and does not raise
    cursor.progress = False
    cursor['P00750']                  # and honours being switched off
