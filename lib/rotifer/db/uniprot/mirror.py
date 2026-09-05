__doc__ = """
Read UniProt data from a local mirror of the UniProt FTP site.

This module gives cursor-style access to the flat files distributed by
UniProt and mirrored locally, for example by ``rrsw`` using the
configuration in ``etc/rotifer/rrsw/uniprot.yml``. Nothing here ever
contacts the network: see :mod:`rotifer.db.uniprot.webapi` for the
REST API client.

Cursors here find their data by walking the directory tree of the FTP
site, which is what the module is named after: a mirror is identified
by its root directory alone, and each cursor knows the path of the
file it reads relative to that root::

    <root>/knowledgebase/idmapping/idmapping.dat
    <root>/knowledgebase/idmapping/idmapping_selected.tab

A cursor accepts either that root directory or the full path of a data
file, so that uncommon layouts and archived releases can still be read
without reconfiguring the package. Because the layout is UniProt's
own, pointing a cursor at a fresh ``rrsw`` mirror is all that is
needed to follow a new release.

See Also
--------
rotifer.db.ncbi.mirror : the equivalent module for NCBI genome mirrors

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
Uncompressed files are scanned in parallel: the file is cut into fixed
size byte ranges, each aligned to line boundaries, and the ranges are
searched independently by a pool of worker processes. pyarrow parses
those ranges, roughly 1.2 to 1.8 times faster than the standard
library. It is a dependency of rotifer, so that is the normal path;
the standard library scanner remains as a fallback for an environment
missing it, and can be asked for with ``engine='python'``.
Gzip compressed copies cannot be cut this way and are scanned by a
single process, which is several times slower.
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
    'engine': 'auto',
}
config = loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

# Bytes handed to each scanning task. It bounds how much of the file a
# worker holds in memory at once, so it must stay small enough that
# threads * _CHUNK fits comfortably in RAM.
_CHUNK = 1 << 26

#: Position of each column in the tab separated files scanned here.
_FIELDS = {'accession': 0, 'id_type': 1, 'id': 2}

#: Per worker state, set once by :func:`_init_worker`.
_worker = {}

def _has_pyarrow():
    """
    Find whether pyarrow is importable.

    Returns
    -------
    bool
    """
    try:
        import pyarrow  # noqa: F401
        return True
    except ImportError:
        return False

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

def _select_arrow(data, targets, field, names):
    """
    Pick the rows of a block whose selected column is a target, with pyarrow.

    Parameters
    ----------
    data : bytes
        Complete lines, separated by newlines.
    targets : pyarrow.Array
        Identifiers to search.
    field : int
        Zero based position of the column to match.
    names : list of str
        Column names of the file.

    Returns
    -------
    list of list of str
        The matching rows, split into their columns.

    Note
    ----
    Quoting and escaping are switched off: these files are plain tab
    separated text, and UniProt identifiers do contain quotes, which a
    CSV aware parser would otherwise swallow.
    """
    import pyarrow as pa
    import pyarrow.csv as pcsv
    import pyarrow.compute as pc

    table = pcsv.read_csv(
        pa.BufferReader(pa.py_buffer(data)),
        read_options = pcsv.ReadOptions(column_names=names, use_threads=False),
        parse_options = pcsv.ParseOptions(delimiter="\t", quote_char=False, escape_char=False, newlines_in_values=False),
        convert_options = pcsv.ConvertOptions(column_types={ x: pa.string() for x in names }, strings_can_be_null=False),
    )
    table = table.filter(pc.is_in(table.column(field), value_set=targets))
    if not table.num_rows:
        return []
    columns = [ x.to_pylist() for x in table.columns ]
    return [ list(row) for row in zip(*columns) ]

def _read_aligned(path, start, end):
    """
    Read one byte range of a file, aligned to line boundaries.

    Ranges neither overlap nor leave a line out: a range starting
    inside a line skips that line, because the range before it owns
    it, and a range whose end falls inside a line reads on until that
    line is complete.

    Parameters
    ----------
    path : str
        Path of the file to read.
    start, end : int
        Byte offsets delimiting the range.

    Returns
    -------
    bytes
        Complete lines, without a trailing newline.
    """
    with open(path, "rb") as fh:
        if start:
            # Reading from one byte before the range makes this test
            # exact: when that byte is a newline the range already
            # starts on a line boundary and nothing is skipped.
            fh.seek(start - 1)
            fh.readline()
            start = fh.tell()
        if start >= end:
            return b""
        fh.seek(start)
        data = fh.read(end - start)
        # A line crossing the end of the range belongs to this range, so
        # read the rest of it. When the range already ends on a newline
        # there is nothing to finish and the next line is not ours:
        # reading one here would return it twice.
        if data and not data.endswith(b"\n"):
            data += fh.readline()
    return data.rstrip(b"\n")

def _scan_range(args):
    """
    Search one byte range of an uncompressed file.

    Parameters
    ----------
    args : tuple
        ``(path, start, end, targets, field)``, with `targets` a set
        of encoded identifiers.

    Returns
    -------
    list of bytes
        The matching lines, without their trailing newline.
    """
    path, start, end, targets, field = args
    data = _read_aligned(path, start, end)
    return _select(data, targets, field) if data else []

def _init_worker(path, targets, field, names, engine):
    """
    Prepare a worker process to scan one file.

    The query is sent once per worker instead of once per task, which
    matters when a large file is cut into many tasks and the query
    carries thousands of identifiers.

    Parameters
    ----------
    path : str
        Path of the file to scan.
    targets : list of str
        Identifiers to search.
    field : int
        Zero based position of the column to match.
    names : list of str
        Column names of the file.
    engine : str
        Either ``arrow`` or ``python``.
    """
    _worker['path'] = path
    _worker['field'] = field
    _worker['names'] = names
    _worker['engine'] = engine
    if engine == 'arrow':
        import pyarrow as pa
        _worker['targets'] = pa.array(sorted(targets), type=pa.string())
    else:
        _worker['targets'] = { x.encode() for x in targets }

def _scan_task(bounds):
    """
    Search one byte range, using the state left by :func:`_init_worker`.

    Parameters
    ----------
    bounds : tuple of int
        The ``(start, end)`` offsets of the range.

    Returns
    -------
    list of list of str
        The matching rows, split into their columns.
    """
    data = _read_aligned(_worker['path'], *bounds)
    if not data:
        return []
    if _worker['engine'] == 'arrow':
        return _select_arrow(data, _worker['targets'], _worker['field'], _worker['names'])
    return [ x.decode().split("\t") for x in _select(data, _worker['targets'], _worker['field']) ]

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
        block = stream.read(_CHUNK)
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
    engine : str, optional
        How each byte range is matched:

        ``auto``
            Use ``arrow`` when pyarrow is importable, ``python``
            otherwise. This is the default, and picks ``arrow``
            wherever rotifer's dependencies are satisfied.
        ``arrow``
            Parse with pyarrow, roughly 1.2 to 1.8 times faster than
            ``python``. Raises an error when pyarrow is missing.
        ``python``
            Match with the standard library alone.

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
    rotifer.db.uniprot.mirror.IdMappingCursor : identifier mapping cursor
    """

    #: Path of the data file, relative to the root of the mirror.
    _datafile = None

    def __init__(self, path=config['local_database_path'], threads=config['threads'], engine=config['engine'], progress=False, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.path = path
        self.threads = max(1, int(threads or 1))
        self.engine = engine
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
        Find every row whose column `column` is one of `targets`.

        Uncompressed files are cut into fixed size byte ranges, each
        aligned to line boundaries, and the ranges are searched in
        parallel by :attr:`threads` worker processes. Gzip compressed
        files are searched sequentially, since their byte ranges
        cannot be decoded independently.

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
            If `column` is not a column of the file, or if
            ``engine='arrow'`` was asked for and pyarrow is missing.
        """
        if column not in _FIELDS:
            raise ValueError(f'Unknown column {column}: expected one of {", ".join(_FIELDS)}')
        field = _FIELDS[column]
        names = [ x for x, _ in sorted(_FIELDS.items(), key=lambda kv: kv[1]) ]
        engine = self._engine()

        # Compressed files decode as one stream, so they cannot be cut up
        if self.compressed:
            with self.open("rb") as fh:
                found = _scan_stream(fh, { x.encode() for x in targets }, field)
            return [ x.decode().split("\t") for x in found ]

        size = os.path.getsize(self.datafile)
        bounds = [ (x, min(x + _CHUNK, size)) for x in range(0, size, _CHUNK) ] or [(0, 0)]
        workers = min(self.threads, len(bounds))

        # One task only, or a single worker: stay in this process and
        # skip the cost of starting a pool
        if workers <= 1:
            _init_worker(self.datafile, list(targets), field, names, engine)
            found = []
            for pair in bounds:
                found += _scan_task(pair)
            return found

        found = []
        with ProcessPoolExecutor(
                max_workers = workers,
                initializer = _init_worker,
                initargs = (self.datafile, list(targets), field, names, engine),
            ) as pool:
            for part in pool.map(_scan_task, bounds):
                found += part
        return found

    def _engine(self):
        """
        Decide which scanning engine to use.

        Returns
        -------
        str
            Either ``arrow`` or ``python``.

        Raises
        ------
        ValueError
            If an unknown engine was requested, or if ``arrow`` was
            requested and pyarrow is not installed.
        """
        if self.engine == 'auto':
            return 'arrow' if _has_pyarrow() else 'python'
        if self.engine == 'arrow':
            if not _has_pyarrow():
                raise ValueError("engine='arrow' requires pyarrow, which is not installed")
            return 'arrow'
        if self.engine == 'python':
            return 'python'
        raise ValueError(f"Unknown engine {self.engine}: expected 'auto', 'arrow' or 'python'")

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
    engine : str, optional
        Matching engine, one of ``auto``, ``arrow`` or ``python``.
        See :class:`BaseUniProtFileCursor`.
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

    >>> from rotifer.db.uniprot import mirror as rum
    >>> ic = rum.IdMappingCursor("/scratch/global/databases/uniprot")  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP

    Find the UniProtKB accession of a RefSeq protein:

    >>> ic = rum.IdMappingCursor(column='id', id_type='RefSeq')  # doctest: +SKIP
    >>> ic.fetchall(["YP_031579.1"])  # doctest: +SKIP
    """

    _datafile = os.path.join("knowledgebase","idmapping","idmapping.dat")

    def __init__(
            self,
            path = config['local_database_path'],
            column = 'accession',
            id_type = None,
            threads = config['threads'],
            engine = config['engine'],
            progress = False,
            *args, **kwargs
        ):
        super().__init__(path=path, threads=threads, engine=engine, progress=progress, *args, **kwargs)
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

        >>> from rotifer.db.uniprot import mirror as rum
        >>> ic = rum.IdMappingCursor()  # doctest: +SKIP
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
        >>> from rotifer.db.uniprot import mirror as rum
        >>> ic = rum.IdMappingCursor()  # doctest: +SKIP
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
        >>> from rotifer.db.uniprot import mirror as rum
        >>> ic = rum.IdMappingCursor()  # doctest: +SKIP
        >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP
        """
        return self.__getitem__(accessions)

class CrossReferenceCursor(IdMappingCursor):
    """
    Find the UniProtKB accessions of identifiers from other databases.

    This is :class:`IdMappingCursor` searching the third column of
    ``idmapping.dat`` instead of the first, which costs exactly the
    same, since either way the whole file is scanned.

    Parameters
    ----------
    path : str, optional
        Root directory of the local UniProt mirror or the full path
        of an ``idmapping.dat`` file.
    id_type : str or list of str, optional
        Restrict the search to these cross-referenced databases.
    threads : int, optional
        Number of worker processes used to scan the file.
    engine : str, optional
        Matching engine, one of ``auto``, ``arrow`` or ``python``.
    progress : bool, default False
        Whether to print progress messages.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.CrossReferenceCursor : same query, indexed and fast

    Examples
    --------
    >>> from rotifer.db.uniprot import mirror as rum
    >>> xc = rum.CrossReferenceCursor(id_type='RefSeq')  # doctest: +SKIP
    >>> xc.fetchall(["YP_031579.1"])  # doctest: +SKIP
    """
    def __init__(self, path=config['local_database_path'], id_type=None,
                 threads=config['threads'], engine=config['engine'], progress=False, *args, **kwargs):
        kwargs.pop('column', None)
        super().__init__(path=path, column='id', id_type=id_type, threads=threads,
                         engine=engine, progress=progress, *args, **kwargs)

class MappingCursor(IdMappingCursor):
    """
    Translate identifiers from one database into another.

    This is the query UniProt's online ID mapping service answers.
    Without an index it takes two passes over the file: the first
    finds the UniProtKB accessions of the queried identifiers, the
    second collects the identifiers those accessions have in the
    target database. Expect it to cost twice a plain lookup.

    Parameters
    ----------
    source : str
        Name of the database the queried identifiers belong to, as
        written in ``idmapping.dat``, e.g. ``EMBL-CDS``. Use
        ``UniProtKB-AC`` to start from UniProtKB accessions
        themselves, which skips the first pass.
    target : str
        Name of the database to translate into, e.g. ``RefSeq``. Use
        ``UniProtKB-AC`` to translate into UniProtKB accessions,
        which skips the second pass.
    path : str, optional
        Root directory of the local UniProt mirror or the full path
        of an ``idmapping.dat`` file.
    threads : int, optional
        Number of worker processes used to scan the file.
    engine : str, optional
        Matching engine, one of ``auto``, ``arrow`` or ``python``.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    columns : list of str
        ``['from', 'accession', 'to']``.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.MappingCursor : same query, as one server side join

    Examples
    --------
    >>> from rotifer.db.uniprot import mirror as rum
    >>> mc = rum.MappingCursor(source='EMBL-CDS', target='RefSeq')  # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"])  # doctest: +SKIP
    """

    _columns = ['from','accession','to']

    def __init__(self, source, target, path=config['local_database_path'],
                 threads=config['threads'], engine=config['engine'], progress=False, *args, **kwargs):
        kwargs.pop('column', None)
        kwargs.pop('id_type', None)
        super().__init__(path=path, column='from', id_type=None, threads=threads,
                         engine=engine, progress=progress, *args, **kwargs)
        self.source = source
        self.target = target

    def _rows(self, targets, column):
        """
        Scan the file and return well formed rows as a dataframe.

        Parameters
        ----------
        targets : set of str
            Identifiers to search.
        column : str
            Name of the column to match.

        Returns
        -------
        pandas.DataFrame
            Columns ``accession``, ``id_type`` and ``id``.
        """
        names = [ x for x, _ in sorted(_FIELDS.items(), key=lambda kv: kv[1]) ]
        if not targets:
            return pd.DataFrame([], columns=names)
        rows = [ x for x in self.scan(targets, column) if len(x) == len(names) ]
        return pd.DataFrame(rows, columns=names)

    def __getitem__(self, accessions):
        """
        Translate identifiers, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers of the database named by ``source``.

        Returns
        -------
        pandas.DataFrame
            Columns ``from``, ``accession`` and ``to``. Identifiers
            with no translation are registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        if isinstance(self.datafile, types.NoneType):
            self.update_missing(targets, error=f'No idmapping file found under {self.path}', retry=False)
            return self.empty()

        # First pass: the UniProtKB accessions of the queried identifiers
        if self.source == "UniProtKB-AC":
            pairs = pd.DataFrame({'from': sorted(targets), 'accession': sorted(targets)})
        else:
            if self.progress:
                logger.warn(f'Scanning {self.datafile} for {len(targets)} {self.source} identifier(s)...')
            found = self._rows(targets, 'id')
            found = found[(found.id_type == self.source) & found.id.isin(targets)]
            pairs = found[['id','accession']].rename(columns={'id':'from'}).drop_duplicates()

        # Second pass: what those accessions are called in the target database
        if self.target == "UniProtKB-AC":
            result = pairs.assign(to=pairs.accession)
        else:
            accessions_found = set(pairs.accession)
            if self.progress:
                logger.warn(f'Scanning {self.datafile} for the {self.target} identifiers of {len(accessions_found)} accession(s)...')
            found = self._rows(accessions_found, 'accession')
            found = found[found.id_type == self.target]
            result = pairs.merge(
                found[['accession','id']].rename(columns={'id':'to'}),
                on = 'accession',
                how = 'inner',
            )

        result = result[self.columns].drop_duplicates().reset_index(drop=True)

        missing = targets.difference(self.getids(result))
        if missing:
            self.update_missing(missing, error=f'No {self.target} identifier found for this {self.source} identifier', retry=False)

        return result

if __name__ == '__main__':
    pass
