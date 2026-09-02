__doc__ = """
Read UniProt data from a local copy of the UniProt FTP site.

This module gives cursor-style access to the flat files distributed by
UniProt and mirrored locally, for example by ``rrsw`` using the
configuration in ``etc/rotifer/rrsw/uniprot.yml``. Nothing here ever
contacts the network: see :mod:`rotifer.db.uniprot.webapi` for the
REST API client.

The mirror is described by a single root directory, whose layout
follows the UniProt FTP site::

    <root>/knowledgebase/idmapping/idmapping.dat
    <root>/knowledgebase/idmapping/idmapping_selected.tab

Cursors accept either the root directory or the full path of a data
file, so that uncommon layouts and archived releases can be used
without reconfiguring the package.

Warning
-------
``idmapping.dat`` is a very large file: the 2026_01 release is about
90 GB uncompressed and holds a few billion rows. The cursors here
answer queries by scanning it, which costs minutes per call. They are
meant for occasional lookups and, above all, for feeding a database
that can index the data. To query the same content interactively, load
it into ClickHouse with :mod:`rotifer.db.uniprot.clickhouse` and use
the cursors defined there.

Notes
-----
Uncompressed files are scanned in parallel: the file is cut into one
byte range per worker process, each range is aligned to line
boundaries and searched independently. Gzip compressed copies cannot
be cut this way and are scanned by a single process, which is several
times slower.
"""

# Dependencies
import os
import types
import typing
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Defaults
_defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"uniprot"),
    'chunksize': 5000000,
    'threads': max(1, (os.cpu_count() or 2) // 2),
}
config = loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

# Size of the blocks read from disk while scanning, in bytes
_BLOCK = 1 << 26

#: Position of each column in the tab separated files scanned here.
_FIELDS = {'accession': 0, 'id_type': 1, 'id': 2}

def _select(block, targets, field):
    """
    Pick the lines of a block whose selected column is a target.

    Parameters
    ----------
    block : bytes
        Complete lines, separated by newlines and without a trailing
        one.
    targets : set of bytes
        Encoded identifiers to search.
    field : int
        Zero based position of the column to match.

    Returns
    -------
    list of bytes
        The matching lines.

    Note
    ----
    The first and last columns are extracted with ``partition`` and
    ``rpartition``, which stop at the first separator found, instead
    of splitting every line into all of its columns.
    """
    if field == 0:
        return [ x for x in block.split(b"\n") if x.partition(b"\t")[0] in targets ]
    elif field == 2:
        return [ x for x in block.split(b"\n") if x.rpartition(b"\t")[2] in targets ]
    else:
        found = []
        for line in block.split(b"\n"):
            columns = line.split(b"\t")
            if len(columns) > field and columns[field] in targets:
                found.append(line)
        return found

def _scan_range(args):
    """
    Search one byte range of an uncompressed file.

    The range is aligned to line boundaries so that ranges neither
    overlap nor leave a line out: a range starting inside a line skips
    that line, because the range before it owns it, and a range whose
    end falls inside a line reads on until that line is complete.

    Parameters
    ----------
    args : tuple
        ``(path, start, end, targets, field)``, packed into a single
        argument so that the function can be used with
        :meth:`concurrent.futures.Executor.map`.

    Returns
    -------
    list of bytes
        The matching lines, without their trailing newline.
    """
    path, start, end, targets, field = args
    found = []
    with open(path, "rb") as fh:
        if start:
            # Reading from one byte before the range makes this test
            # exact: when that byte is a newline the range already
            # starts on a line boundary and nothing is skipped.
            fh.seek(start - 1)
            fh.readline()
            start = fh.tell()
        if start >= end:
            return found

        remaining = end - start
        tail = b""
        while remaining > 0:
            block = fh.read(min(_BLOCK, remaining))
            if not block:
                break
            remaining -= len(block)
            block = tail + block
            cut = block.rfind(b"\n")
            if cut < 0:
                tail = block
                continue
            tail = block[cut+1:]
            found += _select(block[:cut], targets, field)

        # A line crossing the end of the range belongs to this range,
        # so read the rest of it. An empty tail means the range ended
        # exactly on a line boundary and the next line is not ours.
        if tail:
            tail += fh.readline()
            found += _select(tail.rstrip(b"\n"), targets, field)

    return found

def _scan_stream(stream, targets, field):
    """
    Search a whole stream sequentially.

    Used for gzip compressed files, which cannot be cut into
    independent byte ranges.

    Parameters
    ----------
    stream : file-like
        A binary stream positioned at the start of the data.
    targets : set of bytes
        Encoded identifiers to search.
    field : int
        Zero based position of the column to match.

    Returns
    -------
    list of bytes
        The matching lines, without their trailing newline.
    """
    found = []
    tail = b""
    while True:
        block = stream.read(_BLOCK)
        if not block:
            break
        block = tail + block
        cut = block.rfind(b"\n")
        if cut < 0:
            tail = block
            continue
        tail = block[cut+1:]
        found += _select(block[:cut], targets, field)
    if tail:
        found += _select(tail.rstrip(b"\n"), targets, field)
    return found

class BaseUniProtFileCursor(rotifer.db.core.BaseCursor):
    """
    Shared path handling for cursors reading local UniProt files.

    This class is not meant to be used directly: it locates the data
    file a subclass declares in ``_datafile``, opens it, transparently
    handling gzip compressed copies, and scans it for identifiers.
    Failed lookups are tracked through the inherited
    :attr:`~rotifer.db.core.BaseCursor.missing` registry.

    Parameters
    ----------
    path : str, optional
        Either the root directory of the local UniProt mirror or the
        full path of the data file to read. Defaults to the
        ``local_database_path`` configuration entry.
    threads : int, optional
        Number of worker processes used to scan uncompressed files.
        Defaults to the ``threads`` configuration entry, itself half
        of the number of available CPUs. Set to 1 to scan in the
        calling process.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    path : str
        The value received at construction time.
    datafile : str
        Full path of the file the cursor reads.

    See Also
    --------
    rotifer.db.uniprot.io.IdMappingCursor : identifier mapping cursor
    """

    #: Path of the data file, relative to the root of the mirror.
    _datafile = None

    def __init__(self, path=config['local_database_path'], threads=config['threads'], progress=False, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.path = path
        self.threads = max(1, int(threads or 1))
        self.datafile = self._find_datafile(path)

    def _find_datafile(self, path):
        """
        Locate the cursor's data file.

        Parameters
        ----------
        path : str
            Root directory of the local mirror or the full path of a
            data file.

        Returns
        -------
        str or None
            Path of the first candidate that exists, including the
            gzip compressed copy of the expected file. None when no
            candidate is found, in which case an error is logged.
        """
        if isinstance(path, types.NoneType):
            logger.error(f'No path given and no local_database_path configured for {self.__name__}')
            return None

        # A file was given: use it as is
        if os.path.isfile(path):
            return path

        candidates = []
        if not isinstance(self._datafile, types.NoneType):
            candidates.append(os.path.join(path, self._datafile))
            candidates.append(os.path.join(path, self._datafile) + ".gz")
        for candidate in candidates:
            if os.path.isfile(candidate):
                return candidate

        logger.error(f'No such file: {" or ".join(candidates) if candidates else path}')
        return None

    @property
    def compressed(self):
        """
        Whether the data file is gzip compressed.

        Returns
        -------
        bool
        """
        return isinstance(self.datafile, str) and self.datafile.endswith(".gz")

    def open(self, mode="rt"):
        """
        Open the data file.

        Parameters
        ----------
        mode : str, default 'rt'
            File mode, passed on to the underlying opener.

        Returns
        -------
        file-like
            A stream, decompressed on the fly when the file is gzip
            compressed.

        Raises
        ------
        FileNotFoundError
            If the cursor could not locate its data file.
        """
        if isinstance(self.datafile, types.NoneType):
            raise FileNotFoundError(f'{self.__name__}: no data file found under {self.path}')
        if self.compressed:
            import gzip
            return gzip.open(self.datafile, mode)
        return open(self.datafile, mode)

    def scan(self, targets, column):
        """
        Find every line whose column `column` is one of `targets`.

        Uncompressed files are cut into one byte range per worker and
        searched in parallel. Gzip compressed files are searched
        sequentially, since their byte ranges cannot be decoded
        independently.

        Parameters
        ----------
        targets : set of str
            Identifiers to search.
        column : str
            Name of the column to match, a key of :data:`_FIELDS`.

        Returns
        -------
        list of list of str
            The matching rows, split into their columns.

        Raises
        ------
        ValueError
            If `column` is not a column of the file.
        """
        if column not in _FIELDS:
            raise ValueError(f'Unknown column {column}: expected one of {", ".join(_FIELDS)}')
        field = _FIELDS[column]
        encoded = { x.encode() for x in targets }

        if self.compressed:
            with self.open("rb") as fh:
                found = _scan_stream(fh, encoded, field)
        else:
            size = os.path.getsize(self.datafile)
            workers = min(self.threads, max(1, size // _BLOCK)) or 1
            if workers == 1:
                found = _scan_range((self.datafile, 0, size, encoded, field))
            else:
                step = size // workers
                ranges = [
                    (self.datafile, i * step, (i+1) * step if i < workers - 1 else size, encoded, field)
                    for i in range(workers)
                ]
                found = []
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    for part in pool.map(_scan_range, ranges):
                        found += part

        return [ x.decode().split("\t") for x in found ]

class IdMappingCursor(rotifer.db.methods.IdMappingCursor, BaseUniProtFileCursor):
    """
    Fetch UniProt identifier mappings from a local ``idmapping.dat``.

    ``idmapping.dat`` lists, for every UniProtKB accession, the
    identifier of the same protein in each database UniProt
    cross-references. It is a tab separated file with three columns
    and no header: the UniProtKB accession, the name of the
    cross-referenced database and the identifier in that database.

    Queries are answered by scanning the whole file once per call, so
    a single call should carry as many identifiers as possible.

    Parameters
    ----------
    path : str, optional
        Root directory of the local UniProt mirror or the full path
        of an ``idmapping.dat`` file. Defaults to the
        ``local_database_path`` configuration entry.
    column : str, default 'accession'
        Which column the queried identifiers are matched against.
        Use ``accession`` to search UniProtKB accessions and ``id``
        to search the identifiers of cross-referenced databases.
    id_type : str or list of str, optional
        Restrict results to these cross-referenced databases, e.g.
        ``RefSeq`` or ``['EMBL-CDS', 'GeneID']``. By default every
        database is reported.
    threads : int, optional
        Number of worker processes used to scan the file.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    columns : list of str
        ``['accession', 'id_type', 'id']``.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.IdMappingCursor : same data, indexed and fast
    rotifer.db.uniprot.webapi.idmapping : UniProt's online mapping service

    Examples
    --------
    Fetch every cross-reference of two UniProtKB accessions:

    >>> from rotifer.db.uniprot import io as ruio
    >>> ic = ruio.IdMappingCursor("/scratch/global/databases/uniprot")  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP

    Find the UniProtKB accession of a RefSeq protein:

    >>> ic = ruio.IdMappingCursor(column='id', id_type='RefSeq')  # doctest: +SKIP
    >>> ic.fetchall(["YP_031579.1"])  # doctest: +SKIP
    """

    _datafile = os.path.join("knowledgebase","idmapping","idmapping.dat")

    def __init__(
            self,
            path = config['local_database_path'],
            column = 'accession',
            id_type = None,
            threads = config['threads'],
            progress = False,
            *args, **kwargs
        ):
        super().__init__(path=path, threads=threads, progress=progress, *args, **kwargs)
        self.column = column
        self.id_type = id_type
        self.maxgetitem = 1000000

    def _id_types(self):
        """
        Normalize the ``id_type`` filter to a list.

        Returns
        -------
        list of str
            Empty when no filter is set.
        """
        if isinstance(self.id_type, types.NoneType):
            return []
        if isinstance(self.id_type, str) or not isinstance(self.id_type, typing.Iterable):
            return [str(self.id_type)]
        return [ str(x) for x in self.id_type ]

    def reader(self, chunksize=config['chunksize'], id_type=None):
        """
        Iterate over the whole file in chunks.

        This method does not filter by accession: it streams the
        entire mapping table and is the entry point used to load the
        data into other databases.

        Parameters
        ----------
        chunksize : int, optional
            Number of rows per chunk. Defaults to the ``chunksize``
            configuration entry.
        id_type : str or list of str, optional
            Restrict the chunks to these cross-referenced databases.
            Defaults to the cursor's ``id_type`` attribute.

        Yields
        ------
        pandas.DataFrame
            Chunks with the columns listed in :attr:`columns`.

        Examples
        --------
        Count the rows of each cross-referenced database:

        >>> from rotifer.db.uniprot import io as ruio
        >>> ic = ruio.IdMappingCursor()  # doctest: +SKIP
        >>> counts = sum(c.id_type.value_counts() for c in ic.reader())  # doctest: +SKIP
        """
        if isinstance(id_type, types.NoneType):
            id_type = self._id_types()
        elif isinstance(id_type, str):
            id_type = [id_type]

        stream = pd.read_csv(
            self.datafile,
            sep = "\t",
            names = self.columns,
            header = None,
            dtype = str,
            keep_default_na = False,
            na_filter = False,
            chunksize = chunksize,
            compression = "gzip" if self.compressed else None,
        )
        for chunk in stream:
            if id_type:
                chunk = chunk[chunk.id_type.isin(id_type)]
                if chunk.empty:
                    continue
            yield chunk.reset_index(drop=True)

    def __getitem__(self, accessions):
        """
        Fetch identifier mappings, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to search in the column named by the
            cursor's ``column`` attribute.

        Returns
        -------
        pandas.DataFrame
            Mapping rows for the identifiers found, with the columns
            listed in :attr:`columns`. Identifiers that produced no
            row are registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.

        Note
        ----
        This method scans the entire data file, which takes minutes
        for a full ``idmapping.dat``. Batch your queries.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        if isinstance(self.datafile, types.NoneType):
            self.update_missing(targets, error=f'No idmapping file found under {self.path}', retry=False)
            return self.empty()

        if self.progress:
            logger.warn(f'Scanning {self.datafile} for {len(targets)} identifier(s)...')

        rows = [ x for x in self.scan(targets, self.column) if len(x) == len(self.columns) ]
        df = pd.DataFrame(rows, columns=self.columns)

        id_type = self._id_types()
        if id_type and not df.empty:
            df = df[df.id_type.isin(id_type)]
        df = df.reset_index(drop=True)

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'Identifier not found in {self.datafile}', retry=False)

        return df

    def fetchone(self, accessions):
        """
        Iterate over identifier mappings, one query at a time.

        The file is scanned once for the whole batch and the result
        is then split, so this method costs the same as
        :meth:`fetchall`. Input order is not preserved.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to search.

        Yields
        ------
        pandas.DataFrame
            The mapping rows of one identifier.

        Examples
        --------
        >>> from rotifer.db.uniprot import io as ruio
        >>> ic = ruio.IdMappingCursor()  # doctest: +SKIP
        >>> for df in ic.fetchone(["Q6GZX4","Q6GZX3"]):  # doctest: +SKIP
        ...     print(df.accession.iloc[0], len(df))
        """
        found = self.__getitem__(accessions)
        if found.empty:
            return
        for _, block in found.groupby(self.column, sort=False):
            yield block.reset_index(drop=True)

    def fetchall(self, accessions):
        """
        Fetch the identifier mappings of every query at once.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to search.

        Returns
        -------
        pandas.DataFrame
            All mapping rows found, with the columns listed in
            :attr:`columns`.

        Examples
        --------
        >>> from rotifer.db.uniprot import io as ruio
        >>> ic = ruio.IdMappingCursor()  # doctest: +SKIP
        >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP
        """
        return self.__getitem__(accessions)

if __name__ == '__main__':
    pass
