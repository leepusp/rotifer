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

    progress : bool, default True
        Whether to report progress.

    Attributes
    ----------
    path : str
        The value received at construction time.
    datafile : str
        Full path of the file the cursor reads.

    See Also
    --------
    rotifer.db.uniprot.mirror.MappingCursor : identifier mapping cursor
    """

    #: Path of the data file, relative to the root of the mirror.
    _datafile = None

    def __init__(self, path=config['local_database_path'], threads=config['threads'], engine=config['engine'], progress=True, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.path = path
        self.threads = max(1, int(threads or 1))
        self.engine = engine
        self.datafile = self._find_datafile(path)

    def content_id(self):
        """
        Identify the file this cursor reads.

        The value names the data rather than this copy of it, so that
        a backend loaded from the same file reports the same thing and
        can stand in for this one. Size and modification time are what
        make it cheap: identifying a 90 GB file by its contents would
        mean reading all of it, which is the cost the whole comparison
        exists to avoid. Whoever loads the file elsewhere records
        these same three values.

        Returns
        -------
        str or None
            ``<path>:<size>:<mtime>``, or None when no data file was
            found, in which case this cursor is never skipped. The
            path is given in full, with symbolic links resolved: a
            bare name would call two mirrors of different releases the
            same file whenever their sizes and timestamps happened to
            agree, while an unresolved one would call a single file
            two different files whenever it was reached by a link.

        See Also
        --------
        rotifer.db.core.BaseCursor.content_id : what the value means
        """
        if isinstance(self.datafile, types.NoneType):
            return None
        try:
            info = os.stat(self.datafile)
        except OSError:
            logger.debug(f'Cannot stat {self.datafile}', exc_info=1)
            return None
        return f'{os.path.realpath(self.datafile)}:{info.st_size}:{int(info.st_mtime)}'

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

class MappingCursor(rotifer.db.methods.MappingCursor, BaseUniProtFileCursor):
    """
    Translate identifiers by scanning a local ``idmapping.dat``.

    One query answers what used to need three cursors: the
    identifiers given are matched in the databases named by
    ``source``, their UniProtKB accessions are found, and the
    identifiers those accessions carry in the databases named by
    ``target`` are returned. Pinning either end to
    :attr:`~rotifer.db.methods.MappingCursor.UNIPROTKB` gives the
    accession itself.

    Every query scans the whole file, and translating between two
    external databases scans it twice, once per end. A single call
    should therefore carry as many identifiers as possible, and this
    backend is best kept behind one that can answer from an index.

    Parameters
    ----------
    path : str, optional
        Root directory of the local UniProt mirror or the full path
        of an ``idmapping.dat`` file. Defaults to the
        ``local_database_path`` configuration entry.
    threads : int, optional
        Number of worker processes used to scan the file.
    engine : str, optional
        Matching engine, one of ``auto``, ``arrow`` or ``python``.
    progress : bool, default True
        Whether to report progress.

    Note
    ----
    Asking for :attr:`~rotifer.db.methods.MappingCursor.UNIPROTKB` on
    the source side takes the caller's word that the identifiers are
    accessions rather than scanning to confirm it, which would double
    the cost of the commonest query. An accession that does not exist
    therefore yields no rows rather than an error, unless the target
    end is the accession itself, in which case it is echoed back.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.MappingCursor : the same query, indexed

    Examples
    --------
    >>> from rotifer.db.uniprot import mirror as rum
    >>> mc = rum.MappingCursor()                                    # doctest: +SKIP
    >>> mc.fetchall(["Q6GZX4"], source=mc.UNIPROTKB)                # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"], source=['EMBL-CDS'],            # doctest: +SKIP
    ...             target=['RefSeq'])
    """

    _datafile = os.path.join("knowledgebase","idmapping","idmapping.dat")

    #: Columns of the file itself, which are not the columns this
    #: cursor returns: the file stores one row per cross-reference,
    #: while a mapping names both of its ends.
    _table_columns = ['accession','id_type','id']

    def __init__(
            self,
            path = config['local_database_path'],
            threads = config['threads'],
            engine = config['engine'],
            progress = False,
            *args, **kwargs
        ):
        kwargs.pop('column', None)
        kwargs.pop('id_type', None)
        super().__init__(path=path, threads=threads, engine=engine, progress=progress, *args, **kwargs)
        self.maxgetitem = 1000000

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

        Yields
        ------
        pandas.DataFrame
            Chunks with the columns listed in :attr:`_table_columns`.

        Examples
        --------
        >>> from rotifer.db.uniprot import mirror as rum
        >>> mc = rum.MappingCursor()  # doctest: +SKIP
        >>> counts = sum(c.id_type.value_counts() for c in mc.reader())  # doctest: +SKIP
        """
        if isinstance(id_type, str):
            id_type = [id_type]

        stream = pd.read_csv(
            self.datafile,
            sep = "\t",
            names = self._table_columns,
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
            The columns listed in :attr:`_table_columns`.
        """
        if not targets:
            return pd.DataFrame([], columns=self._table_columns)
        rows = [ x for x in self.scan(targets, column) if len(x) == len(self._table_columns) ]
        return pd.DataFrame(rows, columns=self._table_columns)

    def __getitem__(self, accessions, source=None, target=None):
        """
        Translate identifiers, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to translate.
        source, target : str or list of str, optional
            The two ends of the mapping. See
            :meth:`~rotifer.db.methods.MappingCursor.fetchall`.

        Returns
        -------
        pandas.DataFrame
            The columns listed in
            :attr:`~rotifer.db.methods.MappingCursor.columns`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        if isinstance(self.datafile, types.NoneType):
            self.update_missing(targets, error=f'No idmapping file found under {self.path}', retry=False)
            return self.empty()
        source = self.parse_databases(source)
        target = self.parse_databases(target)

        # First pass: what UniProtKB accession each queried identifier
        # belongs to. Accessions are the file's key rather than rows of
        # it, so that branch needs no scan.
        pairs = []
        if source is None or self.UNIPROTKB in source:
            pairs.append(pd.DataFrame({
                'source': sorted(targets),
                'source_type': self.UNIPROTKB,
                'accession': sorted(targets),
            }))
        others = None if source is None else [ x for x in source if x != self.UNIPROTKB ]
        if source is None or others:
            if self.progress:
                logger.warning(f'Scanning {self.datafile} for {len(targets)} identifier(s)...')
            found = self._rows(targets, 'id')
            found = found[found.id.isin(targets)]
            if others:
                found = found[found.id_type.isin(others)]
            pairs.append(found.rename(columns={'id':'source','id_type':'source_type'})
                              [['source','source_type','accession']])
        pairs = pd.concat(pairs, ignore_index=True).drop_duplicates() if pairs else pd.DataFrame(
            [], columns=['source','source_type','accession'])

        # Second pass: what those accessions are called in the target
        # databases. The accession itself is already in hand.
        results = []
        wants_accession = target is not None and self.UNIPROTKB in target
        others = None if target is None else [ x for x in target if x != self.UNIPROTKB ]
        if target is None or others:
            accessions_found = set(pairs.accession)
            if self.progress:
                logger.warning(f'Scanning {self.datafile} for the identifiers of {len(accessions_found)} accession(s)...')
            found = self._rows(accessions_found, 'accession')
            if others:
                found = found[found.id_type.isin(others)]
            results.append(pairs.merge(
                found.rename(columns={'id':'target','id_type':'target_type'})
                     [['accession','target','target_type']],
                on = 'accession', how = 'inner',
            ))
        if wants_accession:
            results.append(pairs.assign(target=pairs.accession, target_type=self.UNIPROTKB))

        if results:
            result = pd.concat(results, ignore_index=True)[self.columns]
            result = result.drop_duplicates().reset_index(drop=True)
        else:
            result = self.empty()

        missing = targets.difference(self.getids(result))
        if missing:
            self.update_missing(missing, error='No mapping found in the local mirror', retry=False)

        return result

    def fetchone(self, accessions, source=None, target=None):
        """
        Iterate over mappings.

        The whole query is answered by one pair of scans, so this
        yields a single dataframe rather than streaming batches:
        splitting the input would mean scanning the file again for
        each batch.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to translate.
        source, target : str or list of str, optional
            The two ends of the mapping.

        Yields
        ------
        pandas.DataFrame
        """
        found = self.__getitem__(accessions, source=source, target=target)
        if not found.empty:
            yield found

if __name__ == '__main__':
    pass
