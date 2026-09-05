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

    Parameters
    ----------
    batch_size: int, default 1
      Number of accessions per batch
    threads: integer, default 3
      Number of simultaneous threads to run
    progress: boolean, deafult False
      Whether to print a progress bar
    tries: int, default 3
      Number of attempts to download each batch

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

    def parse_ids(self, accessions, as_string=True):
        """
        Convert a list of accessions into a set object

        Usage
        -----
        setobj = cursor.parse_ids(["acc1","acc2"])

        Parameters
        ----------
        as_string : boolean, default True
          Convert individual accessions to strings

        Rational
        --------
        Cursors usually receive a list of valid database
        identifiers as input but such a list can be given
        as strings, lists, Pandas series or other objects.

        This method ensures that the input converges into
        a standard representation, i.e. a Python set.
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
        data = { k: self._missing_record(v) for k,v in self._missing.items() }
        return pd.DataFrame(data, index=list(self._missing_fields)).T

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
        retry: boolean, default None
          Filter accessions based on whether this cursor might
          still recover them (retry=True) or not (retry=False)
        final: boolean, default None
          Filter accessions based on whether some cursor reported
          that they will never be recovered, by any backend
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
        accessions: list, tuple or set
          Database identifiers
        error: string, default None
          A string describing the latest error
        retry: boolean, default None
          Whether this cursor might still recover the entry by
          trying again. When not given, it is inferred by matching
          the error against the giveup patterns.
        final: boolean, default False
          Whether the entry will never be recovered, by this or any
          other cursor. Only set this when the cursor has the
          authority to say so: an accession simply absent from one
          database is not final, because another backend may still
          hold it, and a delegator skips every remaining backend for
          entries that are.
        data: dictionary, default None
          A dictionary that matches the internal
          registry, with accessions as keys and
          lists of fields as values. Entries of three
          fields are accepted and read as not final.

          If using this parameter, error, retry and final
          are ignored.
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
        accessions: list of strings or None
          If set to None, all entries in the cache will be removed
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
        Extract accessions from the objects generated by parser().

        Returns
        -------
        A set of strings.
        """
        return NotImplementedError(f'Method getids() must be implemented by descendants')

    def __getitem__(self, accession, *args, **kwargs):
        """
        Fetch data from the database.

        Parameters
        ----------
        accession: string
          Database entry identifier.

        """
        raise NotImplementedError(f'Method __getitem__() must be implemented by descendants')

    def fetchone(self, accessions, *args, **kwargs):
        """
        Asynchronously fetch sequences data from a database.
        Note: input order is not preserved.

        Parameters
        ----------
        accessions: list of strings
          Database entry identifiers.

        Returns
        -------
        A generator for Bio.SeqRecord objects
        """
        raise NotImplementedError(f'Method fetchone() must be implemented by descendants')

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch all data.
        Note: input order is not preserved.

        Parameters
        ----------
        accessions: list of strings
         Database entry identifiers 
        """
        raise NotImplementedError(f'Method fetchall() must be implemented by descendants')

if __name__ == '__main__':
    pass
