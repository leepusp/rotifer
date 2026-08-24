__doc__ = """
Fetch gene neighborhoods and identical protein group reports from
local SQLite3 databases.

The databases queried here are built by other rotifer tools (for
example, by caching batches downloaded through
:mod:`rotifer.db.ncbi`) and follow a fixed schema: a ``features``
table of genome annotation rows and, optionally, a ``nr`` table of
non-redundant sequence clusters. This module never contacts the
network; all cursors here only read and write the SQLite3 file
given as ``path``.
"""

# Dependencies
import re
import os
import sys
import uuid
import types
import typing
import sqlite3
import numpy as np
import pandas as pd
from tqdm import tqdm

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.ncbi.utils as rdnu
import rotifer.devel.beta.sequence as rdbs
from rotifer.core import functions as rcf
from rotifer.genome.data import NeighborhoodDF
logger = rotifer.logging.getLogger(__name__)
config = rcf.loadConfig(__name__, defaults = {})

class BaseSQLite3Cursor(rotifer.db.core.BaseCursor):
    """
    Shared connection and query helpers for SQLite3 cursors.

    This class is not meant to be used directly: it is a parent
    class that opens the database file and provides the temporary
    query table used by subclasses to submit batches of accessions.
    Failed lookups are still tracked through the inherited
    :attr:`~rotifer.db.core.BaseCursor.missing` registry.

    Parameters
    ----------
    path : str
        Path to a SQLite3 database file. The file is created on
        first write if it does not exist.
    replace : bool, default False
        If True, delete any existing file at ``path`` before
        opening the connection.

    See Also
    --------
    rotifer.db.sql.sqlite3.GeneNeighborhoodCursor : gene neighborhood cursor
    rotifer.db.sql.sqlite3.IPGCursor : identical protein group cursor
    """
    def __init__(
            self,
            path,
            replace = False,
            *args, **kwargs
        ):

        super().__init__(*args, **kwargs)
        self.path = path
        self.replace = replace
        if os.path.exists(self.path) and self.replace:
            os.remove(self.path)
        self._dbconn = sqlite3.connect(self.path)
        self.uuid = str(uuid.uuid4())

    def stored(self, data, column='block_id', table='features'):
        """
        Find which rows of the input data are already stored.

        Parameters
        ----------
        data : str, list, pandas.Series or pandas.DataFrame
            Input data to scan for entries in the database.
        column : str or list of str
            Column(s) to use while searching. Ignored when ``data``
            is not a dataframe.
        table : str, default 'features'
            Name of the table to search for matches.

        Returns
        -------
        pandas.Series of bool
            True for rows whose value is already present in
            ``table``.

        Note
        ----
        The default ``column`` is documented as a single column
        name but is iterated as ``for col in column``; see
        ``docs/OPEN_QUESTIONS.md``.
        """
        ret = pd.Series([ False for x in range(1,len(data)) ])
        if not self.has_table(table):
            return ret
        if isinstance(data, str):
            data = [data]
        if isinstance(data, list):
            data = pd.Series(data, name='input')
        for col in column:
            inStore = pd.read_sql(f"""SELECT DISTINCT {col} from {table}""", self._dbconn)[col]
            if isinstance(data, pd.DataFrame):
                ret = ret | data[col].isin(inStore)
            else:
                ret = ret | data.isin(inStore)
        return ret

    def has_table(self, name):
        """
        Find whether a table exists in the database.

        Parameters
        ----------
        name : str
            Table name.

        Returns
        -------
        bool
        """
        sql = self._dbconn.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{name}'").fetchall()
        return len(sql) > 0

    @property
    def schema(self):
        """
        The SQL statement used to create the database's first table.

        Returns
        -------
        str
        """
        return self._dbconn.execute("""SELECT sql FROM sqlite_schema;""").fetchall()[0][0]

    def submit(self, accessions):
        """
        Register query accessions in a temporary table.

        Subclasses join this table against ``features`` (or another
        stored table) to restrict SQL queries to the requested
        accessions.

        Parameters
        ----------
        accessions : iterable
            Any values supported by SQLite3, such as strings,
            integers or floats.

        Note
        ----
        Accessions submitted by a previous call are cleared first;
        this cursor holds only one pending batch at a time.
        """
        if isinstance(accessions,str) or not isinstance(accessions,typing.Iterable):
            ids = set([accessions])
        else:
            ids = set(accessions)
        cursor = self._dbconn.cursor()
        cursor.execute('CREATE TEMPORARY TABLE IF NOT EXISTS queries (id TEXT, uuid TEXT)')
        self.cleanup()
        cursor.executemany("INSERT INTO queries VALUES (?,?)",[ (x,self.uuid) for x in ids ])
        cursor.execute("CREATE INDEX IF NOT EXISTS uidx ON queries (uuid)")
        self._dbconn.commit()

    def cleanup(self):
        """
        Remove this cursor's rows from the temporary query table.
        """
        self._dbconn.execute(f"DELETE FROM queries WHERE uuid = '{self.uuid}'")
        self._dbconn.commit()

class GeneNeighborhoodCursor(rotifer.db.methods.GeneNeighborhoodCursor, BaseSQLite3Cursor):
    """
    Fetch gene neighborhoods cached in a local SQLite3 database.

    All parameters accepted at initialization are also exposed as
    mutable attributes that can be changed between calls to tune the
    cursor's behaviour. Failed lookups are tracked through the
    inherited :attr:`~rotifer.db.core.BaseCursor.missing` registry.

    Parameters
    ----------
    path : str
        Path to a SQLite3 database file.
    replace : bool, default False
        If True, overwrite the database file.
    column : str, default 'pid'
        Name of the ``features`` table column to match against the
        queried accessions. See :class:`rotifer.genome.data.NeighborhoodDF`.
    before : int, default 7
        Keep at most this number of features, of the same type as
        the target, upstream of each target.
    after : int, default 7
        Keep at most this number of features, of the same type as
        the target, downstream of each target.
    min_block_distance : int, default 0
        Minimum distance between two consecutive blocks.
    strand : str, optional
        How to evaluate rows relative to the strand of the target.
        One of:

        - None : ignore strand
        - same : same strand as the target
        - + : positive strand features and targets only
        - - : negative strand features and targets only

    fttype : str, default 'same'
        How to process feature types when counting neighbors:

        - same : consider only features of the same type as the target
        - any : ignore feature type when setting neighborhood
          boundaries

    eukaryotes : bool, default False
        Whether the queried genomes are eukaryotic.
    exclude_type : list of str, default ['source', 'gene', 'mRNA']
        Feature types to ignore.
    autopid : bool, default False
        Automatically set protein identifiers.
    codontable : str or int, default 'Bacterial'
        Default codon table, used when not set in the data.
    progress : bool, default False
        Whether to print a progress bar.

    See Also
    --------
    rotifer.db.methods.GeneNeighborhoodCursor : shared gene neighborhood interface
    rotifer.db.sql.sqlite3.IPGCursor : identical protein group cursor

    Examples
    --------
    Using the dictionary-like interface, fetch the gene
    neighborhood around the gene encoding a target protein:

    >>> from rotifer.db.sql import sqlite3 as rdss
    >>> gnc = rdss.GeneNeighborhoodCursor("genomes.sqlite3")  # doctest: +SKIP
    >>> df = gnc["EEE9598493.1"]  # doctest: +SKIP

    Fetch all gene neighborhoods for a sample of proteins:

    >>> q = ['WP_012291365.1','WP_013208129.1','WP_122330970.1']
    >>> df = gnc.fetchall(q)  # doctest: +SKIP
    """
    def __init__(
            self,
            path,
            replace = False,
            identical = None,
            identical_column = 'c100i100',
            column = 'pid',
            before = 7,
            after = 7,
            min_block_distance = 0,
            strand = None,
            fttype = 'same',
            eukaryotes=False,
            exclude_type=['source','gene','mRNA'],
            autopid=False,
            codontable='Bacterial',
            progress=False,
            *args, **kwargs
        ):

        super().__init__(path=path, replace=replace, progress=progress, *args, **kwargs)
        self.column = column
        self.before = before
        self.after = after
        self.min_block_distance = min_block_distance
        self.strand = strand
        self.fttype = fttype
        self.eukaryotes = eukaryotes
        self.exclude_type = exclude_type
        self.autopid = autopid
        self.codontable = codontable

    def __getitem__(self, accession, ipgs=None):
        """
        Fetch gene neighborhoods, dictionary style.

        Parameters
        ----------
        accession : str or iterable of str
            Database identifiers.
        ipgs : pandas.DataFrame, optional
            Identical protein group report used to restrict results
            to nucleotides confirmed by NCBI's IPG database.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF
            Neighborhood rows for the found accessions. Accessions
            not found in the ``features`` table are registered in
            ``self.missing``.
        """
        if not self.has_table('features'):
            return NeighborhoodDF()
        if not isinstance(accession, typing.Iterable) or isinstance(accession,str):
            accession = [accession]

        # Register queries in the database and search neighborhoods
        self.submit(accession)
        sqlquery = f"""
            SELECT f3.nucleotide, f3.start, f3.end, f3.strand, t.block_id,
                   CASE WHEN t.ids LIKE "%" || f3.pid || "%" THEN 1 ELSE 0 END as query,
                   f3.pid, f3.type, f3.plen, f3.locus, f3.seq_type, f3.assembly, gene, f3.origin,
                   f3.topology, f3.product, f3.organism, f3.lineage, f3.classification,
                   f3.feature_order, f3.internal_id, f3.pid as replaced
            FROM (
                SELECT assembly, nucleotide, type, block_id, min(idup) as idup, max(iddown) as iddown, group_concat(id,char(1)) as ids
                FROM (
                    SELECT *, SUM(nooverlap) OVER (ORDER BY assembly, nucleotide, idup, iddown) as block_id
                    FROM (
                        SELECT *,
                               CASE WHEN 
                                 nucleotide = LAG(nucleotide) OVER (ORDER BY assembly, nucleotide, idup, iddown)
                                 and idup - LAG(iddown) OVER (ORDER BY assembly, nucleotide, idup, iddown) <= {self.min_block_distance}
                                THEN 0
                                ELSE 1
                               END AS nooverlap
                        FROM (
                            SELECT q.id, f1.assembly, f1.nucleotide, f1.type,
                                   f1.feature_order - {self.before} as foup, f1.feature_order + {self.after} as fodown,
                                   min(f2.internal_id) as idup, max(f2.internal_id) as iddown
                            FROM queries as q
                             inner join features as f1 on (q.id = f1.{self.column})
                             inner join features as f2 on (
                                f1.assembly = f2.assembly
                                and f1.nucleotide = f2.nucleotide
                                and f1.type == f2.type
                                and f2.feature_order >= foup
                                and f2.feature_order <= fodown
                             )
                            WHERE q.uuid = '{self.uuid}'
                            GROUP BY q.id, f1.assembly, f1.nucleotide, f1.type, foup, fodown
                            ORDER BY f1.assembly, f1.nucleotide, idup, iddown
                        ) as v
                    ) as w
                ) as z
                GROUP BY assembly, nucleotide, type, block_id
            ) as t
            inner join features as f3 on (
               t.assembly = f3.assembly 
               and t.nucleotide = f3.nucleotide
               and f3.internal_id >= idup 
               and f3.internal_id <= iddown
            )
            WHERE f3.type NOT IN ('{"','".join(self.exclude_type)}')
            ORDER BY f3.assembly, f3.nucleotide, f3.block_id, f3.start, f3.end
        """
        df = pd.read_sql(sqlquery, self._dbconn)
        self.cleanup()

        # Restrict results to nucleotides found in the IPG reports
        if not (isinstance(ipgs,type(None)) or ipgs.empty):
            df = df[df.nucleotide.isin(ipgs.nucleotide)]

        # Find missing entries, if any
        missing = set(accession).difference(self.getids(df, ipgs=ipgs))
        if len(missing):
            self.update_missing(missing, error=f'Entry not found in SQLite3 database {self.path}', retry=True)

        return NeighborhoodDF(df)

    def fetchone(self, proteins, ipgs=None):
        """
        Iterate over gene neighborhoods, one block at a time.

        Parameters
        ----------
        proteins : list of str
            Database identifiers.
        ipgs : pandas.DataFrame, optional
            Identical protein group report, passed through to
            :meth:`__getitem__` to avoid recomputing it.

        Yields
        ------
        rotifer.genome.data.NeighborhoodDF
            The rows of one neighborhood block.

        Examples
        --------
        >>> from rotifer.db.sql import sqlite3 as rdss
        >>> gnc = rdss.GeneNeighborhoodCursor("genomes.sqlite3", progress=True)  # doctest: +SKIP
        >>> for n in gnc.fetchone(['WP_063732599.1']):  # doctest: +SKIP
        ...     print(n.groupby('nucleotide').block_id.nunique())
        """
        if not isinstance(proteins,typing.Iterable) or isinstance(proteins,str):
            proteins = [proteins]
        if self.progress:
            logger.warn(f'Searching {len(proteins)} protein(s) in SQLite3 database at {self.path}')
            p = tqdm(total=len(proteins), initial=0)
        found = self.__getitem__(proteins, ipgs=ipgs)
        for bid, block in found.groupby('block_id'):
            done = self.getids(block, ipgs)
            done = proteins.intersection(done)
            if self.progress and len(done) > 0:
                p.update(len(done))
            yield block.copy()

    def fetchall(self, ids, ipgs=None):
        """
        Fetch all gene neighborhoods at once.

        Parameters
        ----------
        ids : list of str
            Database identifiers.
        ipgs : pandas.DataFrame, optional
            Identical protein group report, passed through to
            :meth:`__getitem__` to avoid recomputing it.

        Returns
        -------
        rotifer.genome.data.NeighborhoodDF

        Examples
        --------
        >>> from rotifer.db.sql import sqlite3 as rdss
        >>> gnc = rdss.GeneNeighborhoodCursor("genomes.sqlite3")  # doctest: +SKIP
        >>> n = gnc.fetchall(['WP_063732599.1'])  # doctest: +SKIP
        """
        return self.__getitem__(ids, ipgs=ipgs)

    def insert(self, data):
        """
        Store genome annotation data in the SQLite3 database.

        Rows already present in the ``features`` table (matched by
        ``block_id``) are skipped.

        Parameters
        ----------
        data : rotifer.genome.data.NeighborhoodDF
            Gene neighborhood dataframe.
        """
        data = data[~self.stored(data)]
        if len(data) > 0:
            data.to_sql("features", self._dbconn, if_exists = 'append', index=False)

class IPGCursor(rotifer.db.methods.SequenceCursor, BaseSQLite3Cursor):
    """
    Fetch identical protein group (IPG) reports from a local SQLite3 database.

    Depending on ``identical``, results are built either from a
    precomputed non-redundant clustering table (``nr``) or directly
    from the ``features`` table.

    Parameters
    ----------
    path : str
        Path to a local SQLite3 database.
    replace : bool, default False
        If True, overwrite the database file.
    identical : str or None, default 'c100'
        Name of the clustering column used to group identical
        proteins. If None, IPGs are derived from ``features`` alone.
    identical_column : str, default 'c100i100'
        Column added to the report identifying the cluster
        representative.

    See Also
    --------
    rotifer.db.sql.sqlite3.GeneNeighborhoodCursor : gene neighborhood cursor
    rotifer.db.ncbi.entrez.IPGCursor : equivalent cursor backed by NCBI Entrez

    Examples
    --------
    >>> from rotifer.db.sql import sqlite3 as rdss
    >>> ic = rdss.IPGCursor("genomes.sqlite3")  # doctest: +SKIP
    >>> df = ic.fetchall("YP_009724395.1")  # doctest: +SKIP
    """
    def __init__(
            self,
            path,
            replace=False,
            identical='c100',
            identical_column='c100i100',
            *args, **kwargs):
        self._columns = ['id','ipg_source','nucleotide','start','stop','strand','pid','description','ipg_organism','strain','assembly']
        self._added_columns = ['order','is_query','representative']
        super().__init__(path=path, replace=replace, *args, **kwargs)
        self.identical = identical
        self.identical_column = identical_column

    # Fetch identical sequences and merge
    def __getitem__(self, accessions):
        """
        Fetch identical protein group reports, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Database identifiers.

        Returns
        -------
        pandas.DataFrame
            Columns as listed in ``self._columns`` plus
            ``self._added_columns``. Empty if the SQL file for the
            selected ``identical`` mode could not be found.
        """
        self.submit(accessions)
        if self.identical == None:
            sql = "ipgs_from_features.sql"
        else:
            sql = "ipgs_from_nr.sql"
        sqlfile = rcf.findDataFiles(__name__ + "." + sql)
        if not len(sqlfile):
            logger.error(f"Could not load SQL file {sql}")
            return pd.DataFrame([], columns=self._columns + self._added_columns)
        sql = " ".join(open(sqlfile,"rt").readlines())
        sql = sql.format(uuid=self.uuid, path=self.path, identical=self.identical, identical_column=self.identical_column)
        result = pd.read_sql(sql, self._dbconn)
        self.cleanup()
        return result

    def fetchall(self, accessions):
        """
        Fetch identical protein group reports for all accessions at once.

        Parameters
        ----------
        accessions : str or iterable of str
            Database identifiers.

        Returns
        -------
        pandas.DataFrame
        """
        return self.__getitem__(accessions)
