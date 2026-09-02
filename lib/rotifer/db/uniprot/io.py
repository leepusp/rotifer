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
"""

# Dependencies
import os
import types
import typing
import tempfile
import subprocess
import pandas as pd

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
}
config = loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

class BaseUniProtFileCursor(rotifer.db.core.BaseCursor):
    """
    Shared path handling for cursors reading local UniProt files.

    This class is not meant to be used directly: it locates the data
    file a subclass declares in ``_datafile`` and opens it,
    transparently handling gzip compressed copies. Failed lookups are
    tracked through the inherited
    :attr:`~rotifer.db.core.BaseCursor.missing` registry.

    Parameters
    ----------
    path : str, optional
        Either the root directory of the local UniProt mirror or the
        full path of the data file to read. Defaults to the
        ``local_database_path`` configuration entry.
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

    def __init__(self, path=config['local_database_path'], progress=False, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.path = path
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

    def open(self):
        """
        Open the data file as a text stream.

        Returns
        -------
        file-like
            A text mode stream, decompressed on the fly when the file
            is gzip compressed.

        Raises
        ------
        FileNotFoundError
            If the cursor could not locate its data file.
        """
        if isinstance(self.datafile, types.NoneType):
            raise FileNotFoundError(f'{self.__name__}: no data file found under {self.path}')
        if self.compressed:
            import gzip
            return gzip.open(self.datafile, "rt")
        return open(self.datafile, "rt")

    def _scan_command(self, patternfile):
        """
        Build the shell pipeline used to scan the data file.

        Parameters
        ----------
        patternfile : str
            Path of a file with one fixed string per line, as
            accepted by ``grep -F -f``.

        Returns
        -------
        list of str
            A command suitable for :func:`subprocess.Popen`, running
            through ``/bin/sh`` so that compressed files can be piped
            through ``zcat``.

        Note
        ----
        ``LC_ALL=C`` makes ``grep`` treat the input as bytes, which is
        both correct for these ASCII files and much faster.
        """
        reader = f'zcat -f -- {self.datafile}' if self.compressed else f'cat -- {self.datafile}'
        return ["/bin/sh","-c", f'LC_ALL=C {reader} | LC_ALL=C grep -F -w -f {patternfile}']

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
            progress = False,
            *args, **kwargs
        ):
        super().__init__(path=path, progress=progress, *args, **kwargs)
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
        This method scans the entire data file, which takes several
        minutes for a full ``idmapping.dat``. Batch your queries.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        if isinstance(self.datafile, types.NoneType):
            self.update_missing(targets, error=f'No idmapping file found under {self.path}', retry=False)
            return self.empty()

        if self.progress:
            logger.warn(f'Scanning {self.datafile} for {len(targets)} identifier(s)...')

        rows = []
        with tempfile.NamedTemporaryFile("wt", suffix=".patterns", delete=False) as patternfile:
            patternfile.write("\n".join(sorted(targets)) + "\n")
            patternfile.flush()
            name = patternfile.name
        try:
            process = subprocess.Popen(
                self._scan_command(name),
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                text = True,
            )
            for line in process.stdout:
                row = line.rstrip("\n").split("\t")
                if len(row) == 3:
                    rows.append(row)
            process.stdout.close()
            error = process.stderr.read()
            process.stderr.close()
            # grep exits with status 1 when it matches nothing, which is not an error here
            if process.wait() > 1:
                logger.error(f'Failed to scan {self.datafile}: {error}')
        finally:
            os.unlink(name)

        df = pd.DataFrame(rows, columns=self.columns)

        # grep matches anywhere in the line: keep only the rows
        # where the query is in the column we were asked to search
        if not df.empty:
            df = df[df[self.column].isin(targets)]
            id_type = self._id_types()
            if id_type:
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
