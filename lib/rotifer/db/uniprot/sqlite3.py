__doc__ = """
Query UniProt identifier mappings stored in SQLite3.

This is the same table and the same queries as
:mod:`rotifer.db.uniprot.clickhouse`, in a file rather than on a
server. It exists so that a mapping can be carried with an analysis,
or answered where no server is reachable, without falling back on
scanning ``idmapping.dat``: a scan costs about ninety seconds whatever
is asked of it, and an indexed file answers in milliseconds.

The SQL is not written here. Both backends build it from
:mod:`rotifer.db.sql.mapping`, which knows the shape of a mapping
query and nothing about either engine; this module supplies a table
name, a way of binding values, and the little SQLite3 does
differently.

What SQLite3 does differently
-----------------------------
Values are bound by position, so a statement carries one placeholder
per value and there is a limit on how many, which is why identifiers
are asked for in batches.

There are no projections. The reverse lookup that ClickHouse answers
from a second copy of the data sorted by ``id`` is answered here from
an ordinary index on that column, which costs the same kind of thing:
a second structure, written during the load.

A whole release fits. SQLite3 handles a few billion rows of this shape
perfectly well when it is written for reading, which is what a mapping
table is: loaded once, then never updated. What it needs is a load
written for the purpose, so :meth:`MappingCursor.load` turns off the
journal, defers every index to the end and writes in one transaction.
An ``id_type`` filter is there for the cases where a subset is what is
wanted -- one project, the databases an analysis actually touches --
not because the whole is out of reach.

One release at a time, though. ClickHouse keeps several side by side
because a partition makes that free, and drops one in an instant;
SQLite3 has no such thing, so loading a release here replaces what was
there rather than adding to it. Pass ``replace=False`` to accumulate,
and expect the deletion of a release afterwards to cost a scan.

Examples
--------
Build a file holding only the RefSeq and KEGG mappings, then query it:

>>> from rotifer.db.uniprot import sqlite3 as rus     # doctest: +SKIP
>>> mc = rus.MappingCursor('uniprot.sqlite3')          # doctest: +SKIP
>>> mc.create()                                        # doctest: +SKIP
>>> mc.load('/scratch/global/databases/uniprot',       # doctest: +SKIP
...         release='2026_01', id_type=['RefSeq','KEGG'])
>>> mc.fetchall(['Q6GZX4'], source=mc.UNIPROTKB)       # doctest: +SKIP

See Also
--------
rotifer.db.sql.mapping : the SQL both backends share
rotifer.db.uniprot.clickhouse : the same queries, on a server
"""

# Dependencies
import os
import types
import pandas as pd

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.sql.mapping as sqlmap
from rotifer.db.sql.sqlite3 import BaseSQLite3Cursor
from rotifer.core import functions as rcf
logger = rotifer.logging.getLogger(__name__)

_defaults = {
    'batch_size': 500,
    'chunksize': 1000000,
}
config = rcf.loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)


class MappingCursor(rotifer.db.methods.MappingCursor, BaseSQLite3Cursor):
    """
    Translate identifiers through a mapping table held in SQLite3.

    The query is the one described in
    :mod:`rotifer.db.sql.mapping`: identifiers are matched in the
    databases named by ``source``, their UniProtKB accessions found,
    and the identifiers those accessions carry in the databases named
    by ``target`` returned. Either end may be
    :attr:`~rotifer.db.methods.MappingCursor.UNIPROTKB`, meaning the
    accession itself, and either may be omitted, meaning every
    database.

    Parameters
    ----------
    path : str
        Path of the SQLite3 file.
    release : str, optional
        Restrict results to one UniProt release. A file usually holds
        one, in which case there is nothing to restrict.
    batch_size : int, optional
        Identifiers per statement. SQLite3 limits how many values one
        statement may bind, so a long query is asked in batches.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    columns : list of str
        ``['source', 'source_type', 'accession', 'target',
        'target_type']``.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.MappingCursor : the same, on a server
    """

    #: The one table this cursor reads, and the role it plays.
    tables = {'mapping': 'idmapping'}

    #: Columns of the mapping table itself, which are not the columns
    #: this cursor returns.
    _table_columns = ['accession', 'id_type', 'id']

    def __init__(self, path, release=None, batch_size=config['batch_size'],
                 progress=False, *args, **kwargs):
        super().__init__(path=path, progress=progress, *args, **kwargs)
        self.release = release
        self.batch_size = batch_size
        self.maxgetitem = 1000000

    @property
    def source_version(self):
        """
        Release distinguishing one load of this table from another.

        Returns
        -------
        str
        """
        return self.release or ''

    def qualified(self, role=None):
        """
        Name one of this cursor's tables.

        SQLite3 has one database per file, so a name needs no
        qualifying, which is the only reason this differs from the
        server backend.

        Parameters
        ----------
        role : str, optional
            Which of the cursor's tables.

        Returns
        -------
        str
        """
        return self.table_name(role)

    #: Indexes the queries need, by name and by the columns they cover.
    _indexes = {
        'by_accession': '(accession, id_type)',
        'by_id': '(id, id_type)',
    }

    def create(self, replace=False, indexes=True, role=None):
        """
        Create the mapping table, and by default the indexes too.

        Two indexes, because the two directions read different
        columns: one on the accession, which answers "what is this
        accession called elsewhere", and one on the identifier, which
        answers the reverse. The second is what a ClickHouse
        projection is for, in the form SQLite3 offers.

        Parameters
        ----------
        replace : bool, default False
            Drop any existing table first. Every row it holds is lost.
        indexes : bool, default True
            Build the indexes as well. A bulk load passes False and
            builds them afterwards, which is much faster than
            maintaining them row by row.
        role : str, optional
            Which of the cursor's tables to build.

        Returns
        -------
        bool
            Whether the table exists afterwards.
        """
        table = self.table_name(role)
        if replace:
            self._dbconn.execute(f'DROP TABLE IF EXISTS {table}')
        self._dbconn.execute(f"""
            CREATE TABLE IF NOT EXISTS {table} (
                accession TEXT NOT NULL,
                id_type   TEXT NOT NULL,
                id        TEXT NOT NULL,
                release   TEXT NOT NULL DEFAULT ''
            )
        """)
        if indexes:
            self.create_indexes(role=role)
        self._dbconn.commit()
        return self.has_table(table)

    def create_indexes(self, role=None):
        """
        Build the indexes the queries read.

        Building them after the rows are in is far cheaper than
        keeping them up to date while the rows arrive, which is why a
        load leaves them until last.

        Parameters
        ----------
        role : str, optional
            Which of the cursor's tables.
        """
        table = self.table_name(role)
        for name, columns in self._indexes.items():
            self._dbconn.execute(
                f'CREATE INDEX IF NOT EXISTS {table}_{name} ON {table} {columns}')
        self._dbconn.commit()

    def drop_indexes(self, role=None):
        """
        Remove the indexes, so that rows can be written without them.

        Parameters
        ----------
        role : str, optional
            Which of the cursor's tables.
        """
        table = self.table_name(role)
        for name in self._indexes:
            self._dbconn.execute(f'DROP INDEX IF EXISTS {table}_{name}')
        self._dbconn.commit()

    def has_indexes(self, role=None):
        """
        Find whether the indexes are in place.

        Returns
        -------
        bool
        """
        table = self.table_name(role)
        found = { r[0] for r in self._dbconn.execute(
            "SELECT name FROM sqlite_master WHERE type='index'").fetchall() }
        return all(f'{table}_{name}' in found for name in self._indexes)

    def count(self, role=None):
        """
        Count the rows of the mapping table.

        Parameters
        ----------
        role : str, optional
            Which of the cursor's tables.

        Returns
        -------
        int
        """
        table = self.table_name(role)
        if not self.has_table(table):
            return 0
        return self._dbconn.execute(f'SELECT count(*) FROM {table}').fetchone()[0]

    def insert(self, data, release=None, role=None):
        """
        Insert mapping rows into the table.

        Parameters
        ----------
        data : pandas.DataFrame
            Rows with the columns listed in :attr:`_table_columns`.
        release : str, optional
            Value for the ``release`` column. Defaults to the
            cursor's.
        role : str, optional
            Which of the cursor's tables.
        """
        if data.empty:
            return
        table = self.table_name(role)
        rows = data[self._table_columns].copy()
        rows['release'] = (release if not isinstance(release, types.NoneType)
                           else self.release) or ''
        self._dbconn.executemany(
            f'INSERT INTO {table} (accession, id_type, id, release) VALUES (?, ?, ?, ?)',
            rows.itertuples(index=False, name=None),
        )
        self._dbconn.commit()

    def vacuum(self):
        """
        Compact the file after a load that replaced its contents.

        Deleting a release leaves the pages behind for reuse, which is
        what a database should do and not what a file meant to be
        copied around wants.
        """
        self._dbconn.execute('VACUUM')

    #: Pragmas set while bulk loading, and what each is for. The
    #: journal and the fsync are what a load spends its time on and
    #: neither earns its cost here: the table is written once from a
    #: file that still exists, so an interrupted load is thrown away
    #: and started again rather than rolled back.
    _load_pragmas = {
        'journal_mode': 'OFF',
        'synchronous': 'OFF',
        'temp_store': 'MEMORY',
        'cache_size': '-1048576',      # a gibibyte, negative means KiB
    }

    def _pragmas(self, settings):
        """
        Apply pragmas and return what they were.

        Parameters
        ----------
        settings : dict
            Pragma names and values.

        Returns
        -------
        dict
            The previous values, for putting back.
        """
        previous = {}
        for name, value in settings.items():
            try:
                was = self._dbconn.execute(f'PRAGMA {name}').fetchone()
                previous[name] = was[0] if was else None
                self._dbconn.execute(f'PRAGMA {name} = {value}')
            except Exception:
                logger.debug(f'Could not set PRAGMA {name}', exc_info=1)
        return previous

    def load(self, mirror, release=None, id_type=None, replace=True,
             chunksize=config['chunksize'], tune=True, role=None):
        """
        Fill the table from a local copy of ``idmapping.dat``.

        A release replaces what the table held rather than adding to
        it. ClickHouse keeps several side by side because partitioning
        makes that free and dropping one instant; here it would mean a
        larger file and a scan to undo, so one release at a time is
        the useful default. Pass ``replace=False`` to accumulate.

        The load is written for the purpose: no journal, no fsync, no
        indexes until the rows are in, and one transaction. Those are
        safe here for a reason worth stating -- the table is built once
        from a file that still exists, so a load that dies is discarded
        and repeated, never recovered.

        Parameters
        ----------
        mirror : str or rotifer.db.uniprot.mirror.MappingCursor
            Root of a local UniProt mirror, the path of an
            ``idmapping.dat`` file, or a cursor already reading one.
        release : str, optional
            Value stamped on every row loaded. Defaults to the
            cursor's.
        id_type : str or list of str, optional
            Load only these cross-referenced databases. By default
            every one of them is loaded.
        replace : bool, default True
            Empty the table first, so it holds this release alone.
        chunksize : int, optional
            Rows read from the file at a time.
        tune : bool, default True
            Set the bulk loading pragmas and defer the indexes. Turn
            it off to load into a database something else is using,
            where those settings would not be safe.
        role : str, optional
            Which of the cursor's tables to fill.

        Returns
        -------
        int
            Number of rows in the table afterwards.
        """
        from rotifer.db.uniprot import mirror as rum

        if isinstance(release, types.NoneType):
            release = self.release
        reader = mirror if isinstance(mirror, rum.MappingCursor) else \
            rum.MappingCursor(path=mirror, progress=self.progress)
        if isinstance(reader.datafile, types.NoneType):
            logger.error(f'No idmapping file found for {mirror}')
            return self.count(role)

        table = self.table_name(role)
        self.create(indexes=not tune, role=role)
        if replace:
            self._dbconn.execute(f'DELETE FROM {table}')
            self.forget_source(version='', table=table)

        previous = self._pragmas(self._load_pragmas) if tune else {}
        if tune:
            self.drop_indexes(role=role)
        if self.progress:
            logger.warning(f'Loading {reader.datafile} into {self.path}...')
        try:
            for chunk in reader.reader(chunksize=chunksize, id_type=id_type):
                self.insert(chunk, release=release, role=role)
        finally:
            if tune:
                if self.progress:
                    logger.warning('Building indexes...')
                self.create_indexes(role=role)
                self._pragmas({ k: v for k, v in previous.items()
                                if not isinstance(v, types.NoneType) })

        rows = self.count(role)
        if isinstance(id_type, types.NoneType):
            self.record_source(reader.datafile, rows, version=release, table=table)
        elif self.progress:
            logger.warning(
                'Loaded a subset, so no source is recorded: a file holding '
                'some of a release is not a copy of it'
            )
        return rows

    def __getitem__(self, accessions, source=None, target=None):
        """
        Translate identifiers, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to translate.
        source, target : str or list of str, optional
            The two ends of the mapping.

        Returns
        -------
        pandas.DataFrame
            The columns listed in
            :attr:`~rotifer.db.methods.MappingCursor.columns`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        table = self.table_name()
        if not self.has_table(table):
            self.update_missing(targets, error=f'No table {table} in {self.path}',
                                retry=False)
            return self.empty()
        source = self.parse_databases(source)
        target = self.parse_databases(target)

        stack = []
        try:
            for batch in self._batches(sorted(targets), self.batch_size):
                binder = sqlmap.QmarkBinder()

                def restrict(column, _batch=batch, _binder=binder):
                    return f'f.{column} IN {_binder.collection(_batch)}'

                sql = sqlmap.mapping_query(table, restrict, source, target, binder,
                                           release=self.release)
                found = pd.read_sql(sql, self._dbconn, params=binder.parameters)
                if not found.empty:
                    stack.append(found)
        except Exception as error:
            logger.error(f'Query to {self.path} failed: {error}')
            self.update_missing(targets, error=f'SQLite3 query failed: {error}', retry=True)
            return self.empty()

        df = pd.concat(stack, ignore_index=True) if stack else self.empty()
        if not df.empty:
            df = df[self.columns].drop_duplicates().reset_index(drop=True)

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'No mapping found in {self.path}',
                                retry=False)
        return df

    def _batches(self, values, size):
        """
        Cut a list of identifiers into statement sized pieces.

        Parameters
        ----------
        values : list
        size : int
            Identifiers per batch.

        Yields
        ------
        list
        """
        values = list(values)
        size = max(1, int(size or len(values) or 1))
        for start in range(0, len(values), size):
            yield values[start:start+size]

    def fetchone(self, accessions, source=None, target=None):
        """
        Iterate over mappings, one batch of identifiers at a time.

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
        targets = self.parse_ids(accessions)
        found = self.__getitem__(targets, source=source, target=target)
        if not found.empty:
            yield found


if __name__ == '__main__':
    pass
