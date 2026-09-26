"""
Abstract base class for all database cursors.

This module defines :class:`BaseCursor`, the root of the cursor class
hierarchy used across :mod:`rotifer.db`. It standardizes identifier
handling and the bookkeeping of entries that could not be retrieved.
Concrete cursors implement the actual data access methods
(``__getitem__``, :meth:`~BaseCursor.fetchone` and
:meth:`~BaseCursor.fetchall`).
"""

import types
import typing
import pandas as pd

# Rotifer
import rotifer
import rotifer.core.functions as rcf
logger = rotifer.logging.getLogger(__name__)

class BaseCursor:
    """
    Generic database cursor abstract interface.

    Subclasses must implement ``__getitem__``, :meth:`fetchone`,
    :meth:`fetchall` and :meth:`getids`. The base class provides
    identifier normalization (:meth:`parse_ids`) and a registry of
    missing entries (:attr:`missing`, :meth:`missing_ids`,
    :meth:`update_missing`, :meth:`remove_missing`).

    Parameters
    ----------
    progress : bool, default False
        Whether to print a progress bar.

    Attributes
    ----------
    giveup : set of str
        Substrings of error messages that mark an entry as
        unrecoverable, so no further retrieval attempts are made.
    maxgetitem : int
        Maximum number of accessions accepted per ``__getitem__``
        call.
    """
    def __init__(self, progress=False, *args, **kwargs):
        self.progress = progress
        self.__name__ = str(type(self)).split("'")[1]
        self._missing = dict() # Keys are accessions, values are lists of the fields in _missing_fields
        self.giveup = set() # List of errors that will prevent further attempts to use failed accessions
        # Errors saying the entry will never be found, by this or any
        # other cursor. Matching one implies giving up here too, so
        # these need not be repeated in giveup. Reserve it for
        # statements about the data rather than about one source: a
        # cursor refusing a kind of accession, or rejecting a
        # malformed request, is not speaking for the others.
        self.final_errors = set()
        self.maxgetitem = 1 # Maximum number of arguments accepted by __getitem__()

    @property
    def progress_label(self):
        """
        Name this cursor answers to on a progress bar.

        A delegator draws its own bar above the backend working under
        it, and an unlabelled bar says nothing about which backend
        that is -- which is the one thing worth knowing while several
        are tried in turn. What distinguishes them is the module they
        come from, so that is what the bar is named after, without the
        ``rotifer.db.`` prefix every one of them shares.

        Returns
        -------
        str
            The module name, e.g. ``uniprot.clickhouse``.
        """
        name = type(self).__module__
        prefix = 'rotifer.db.'
        return name[len(prefix):] if name.startswith(prefix) else name

    def parse_ids(self, accessions, as_string=True):
        """
        Convert a collection of accessions into a set.

        Cursors accept identifiers as strings, lists, tuples, sets,
        pandas Series or other iterables. This method converges any
        of those inputs into a standard representation, a Python set.
        Comma separated strings are split into individual accessions.

        Parameters
        ----------
        accessions : str or iterable
            One or more database identifiers. A string containing
            commas is split on the commas.
        as_string : bool, default True
            Convert individual accessions to strings.

        Returns
        -------
        set
            The normalized, deduplicated identifiers.

        Examples
        --------
        >>> from rotifer.db.core import BaseCursor
        >>> cursor = BaseCursor()
        >>> sorted(cursor.parse_ids("acc1,acc2"))
        ['acc1', 'acc2']
        >>> sorted(cursor.parse_ids(["acc1", "acc1", "acc2"]))
        ['acc1', 'acc2']
        """
        from copy import deepcopy
        targets = deepcopy(accessions)
        if isinstance(targets,str):
            targets = targets.split(",") # Useful for Entrez... remove?
        elif not isinstance(targets,typing.Iterable):
            targets = [ targets ]
        if as_string:
            targets = [ str(x) for x in targets ]
        targets = set(targets)
        return targets

    #: Fields kept for every entry of the registry of missing entries.
    _missing_fields = ("error", "class", "retry", "final")

    @staticmethod
    def _missing_record(entry):
        """
        Pad a registry entry to the current number of fields.

        Entries used to hold three fields, and code outside this class
        still builds them that way. A record short of ``final`` is
        read as not final, which is what it meant before the field
        existed.

        Parameters
        ----------
        entry : list
          One entry of the registry.

        Returns
        -------
        list
          The same entry, extended to four fields when needed.
        """
        entry = list(entry)
        while len(entry) < len(BaseCursor._missing_fields):
            entry.append(False)
        return entry

    @property
    def missing(self):
        """
        Registry of entries that could not be retrieved.

        Returns
        -------
        pandas.DataFrame
            One row per missing accession, indexed by accession, with
            columns ``error`` (last error message), ``class`` (the
            cursor that failed) and ``retry`` (whether another attempt
            may succeed).
        """
        return pd.DataFrame(self._missing, index="error class retry".split(" ")).T

    def missing_ids(self, retry=None, final=None):
        """
        Retrieve accessions not found in the target database.

        The two filters answer different questions. `retry` asks
        whether this cursor might still recover an accession by trying
        again, which is what its own retry loop consults. `final` asks
        whether the answer is binding on every other cursor as well,
        which is what a delegator consults before handing the
        accession to the next backend.

        An accession absent from one database is normally neither: not
        worth retrying here, but well worth asking the next backend
        for.

        Parameters
        ----------
        retry : bool, optional
            When set, filter accessions based on whether they might
            still be recovered (``retry=True``) or not
            (``retry=False``). By default all missing accessions are
            returned.

        Returns
        -------
        set of str
            The missing accessions.
        """
        selected = []
        for accession, entry in self._missing.items():
            entry = self._missing_record(entry)
            if not isinstance(retry, types.NoneType) and entry[2] != retry:
                continue
            if not isinstance(final, types.NoneType) and bool(entry[3]) != final:
                continue
            selected.append(accession)
        return set(sorted(selected))

    def update_missing(self, accessions=[], error=None, retry=None, final=False, data=None, *args, **kwargs):
        """
        Update or add entries to the registry of missing entries.

        Parameters
        ----------
        accessions : list, tuple or set
            Database identifiers to register as missing.
        error : str, optional
            A string describing the latest error. When not given,
            the previously recorded error is kept, or the text
            ``"Unknown error"`` is used for new entries.
        retry : bool, optional
            Whether the error is recoverable. When not given, the
            value is inferred by matching `error` against the
            ``giveup`` patterns.
        data : dict, optional
            A dictionary matching the internal registry layout, with
            accessions as keys and three element lists (error,
            calling class, retry flag) as values. When given, the
            `accessions`, `error` and `retry` parameters are ignored.

        Returns
        -------
        bool
            Whether at least one of the registered entries may still
            be recovered by another retrieval attempt.
        """
        if isinstance(data, types.NoneType):
            definitive = False
            gaveup = False
            if not isinstance(error,types.NoneType):
                definitive = any([ x in error for x in self.final_errors ])
                gaveup = definitive or any([ x in error for x in self.giveup ])
            if isinstance(retry,types.NoneType):
                retry = not gaveup
            err = [error, rcf.who_is_calling(self), retry, bool(final) or definitive]
            targets = self.parse_ids(accessions)
            for x in targets:
                if error == None:
                    if x in self._missing:
                        err[0] = self._missing[x][0]
                    else:
                        err[0] = "Unknown error"
                self._missing[x] = err
        else:
            for k,v in data.items():
                entry = self._missing_record(v)
                # A final verdict is permanent by definition, so a
                # later report about the same entry may replace the
                # message but never downgrade it
                previous = self._missing.get(k)
                if previous and self._missing_record(previous)[3]:
                    entry[3] = True
                self._missing[k] = entry
            retry = any([ self._missing_record(v)[2] for k,v in data.items() ])
        return retry

    def remove_missing(self, accessions=None):
        """
        Unregister missing accessions.

        Parameters
        ----------
        accessions : list of str, optional
            Accessions to remove from the registry. When not given,
            all entries are removed.

        Returns
        -------
        dict or None
            When all entries are removed, the previous registry
            content is returned. Otherwise nothing is returned.
        """
        if isinstance(accessions,types.NoneType):
            old = self._missing.copy()
            self._missing = dict()
            return old
        else:
            for k in self.parse_ids(accessions):
                self._missing.pop(k, None)

    def content_id(self):
        """
        Identify the body of data this cursor reads.

        Two cursors that return this same value are serving the same
        data, so consulting both answers nothing the first did not.
        A delegator uses that to skip a backend whose contents another
        one has already offered, which is worth doing when the skipped
        backend is the expensive one.

        This is a statement about data, not about reachability or
        freshness: it says two sources hold the same thing, never that
        what they hold is complete. A backend that cannot cheaply
        prove which data it holds returns None, and None is never
        equal to anything, so it is always consulted.

        Cursors that can answer should build the value from whatever
        identifies the data itself rather than the copy of it: two
        backends serving one release must agree, however differently
        they store it.

        Returns
        -------
        hashable or None
            None by default, i.e. never skipped.
        """
        return None

    def getids(self, obj):
        """
        Extract accessions from the objects generated by the cursor.

        Subclasses must override this method to recognize their own
        result objects.

        Parameters
        ----------
        obj : object
            A result object produced by the cursor.

        Returns
        -------
        set of str
            The accessions found in `obj`.
        """
        return NotImplementedError(f'Method getids() must be implemented by descendants')

    def __getitem__(self, accession, *args, **kwargs):
        """
        Fetch data for one entry, dictionary style.

        Subclasses must override this method.

        Parameters
        ----------
        accession : str
            Database entry identifier.

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.
        """
        raise NotImplementedError(f'Method __getitem__() must be implemented by descendants')

    def fetchone(self, accessions, *args, **kwargs):
        """
        Iterate over entries as they are retrieved.

        Subclasses must override this method. Input order is not
        preserved.

        Parameters
        ----------
        accessions : list of str
            Database entry identifiers.

        Yields
        ------
        object
            One result object per retrieved entry.

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.
        """
        raise NotImplementedError(f'Method fetchone() must be implemented by descendants')

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch data for all entries at once.

        Subclasses must override this method. Input order is not
        preserved.

        Parameters
        ----------
        accessions : list of str
            Database entry identifiers.

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.
        """
        raise NotImplementedError(f'Method fetchall() must be implemented by descendants')

if __name__ == '__main__':
    pass
