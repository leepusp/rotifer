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
        self._missing = dict() # Keys are accessions, values are lists of three elements
        self.giveup = set() # List of errors that will prevent further attempts to use failed accessions
        self.maxgetitem = 1 # Maximum number of arguments accepted by __getitem__()

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

    def missing_ids(self, retry=None):
        """
        Retrieve accessions not found in the target database.

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
        if isinstance(retry, types.NoneType):
            return set(sorted(list(self._missing.keys())))
        else:
            return set(sorted([ x for x in self._missing.keys() if self._missing[x][2] == retry ]))

    def update_missing(self, accessions=[], error=None, retry=None, data=None, *args, **kwargs):
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
            if isinstance(retry,types.NoneType):
                retry = True
                if not isinstance(error,types.NoneType):
                    for x in self.giveup:
                        if x in error:
                            retry = False
                            break
            err = [error, rcf.who_is_calling(self), retry]
            targets = self.parse_ids(accessions)
            for x in targets:
                if error == None:
                    if x in self._missing:
                        err[0] = self._missing[x][0]
                    else:
                        err[0] = "Unknown error"
                self._missing[x] = err
        else:
            self._missing.update(data)
            retry = any([ v[2] for k,v in data.items() ])
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
