__doc__ = """
Query UniProt identifier mappings stored in ClickHouse.

UniProt's ``idmapping.dat`` relates every UniProtKB accession to the
identifier of the same protein in each database UniProt
cross-references. The file is far too large to be searched
interactively (about 90 GB and a few billion rows in the 2026_01
release), so this module keeps it in a ClickHouse table and provides
cursors for the queries that table was designed to answer:

:class:`IdMappingCursor`
    Given UniProtKB accessions, return their cross-references.
:class:`CrossReferenceCursor`
    Given identifiers of other databases, return the UniProtKB
    accessions they belong to.
:class:`MappingCursor`
    Translate identifiers from one database to another, the
    equivalent of UniProt's online ID mapping service.

The table itself is created and populated through
:meth:`IdMappingCursor.create` and :meth:`IdMappingCursor.load`, from
a local copy of the flat file read by
:class:`rotifer.db.uniprot.io.IdMappingCursor`. Its schema lives in
``share/rotifer/db/uniprot/clickhouse/idmapping.sql``.

Configuration
-------------
Connection parameters are read from ``~/.rotifer/etc/db/uniprot/clickhouse.yml``
and default to a server on ``localhost``. Note that the default
``port`` below is ClickHouse's HTTP port, 8123, which is not
necessarily the port a given server listens on.
"""

# Dependencies
import os
import types
import typing
import subprocess
import pandas as pd

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
from rotifer.core import functions as rcf
logger = rotifer.logging.getLogger(__name__)

# Defaults
_defaults = {
    'host': 'localhost',
    'port': 8123,
    'user': 'default',
    'password': '',
    'database': 'uniprot',
    'table': 'idmapping',
    'release': '',
    'batch_size': 10000,
    'chunksize': 5000000,
    'executable': 'clickhouse',
}
config = rcf.loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

class BaseClickHouseCursor(rotifer.db.core.BaseCursor):
    """
    Shared connection and query helpers for ClickHouse cursors.

    This class is not meant to be used directly: it holds the
    connection parameters, opens the client on first use and wraps the
    few server calls the subclasses need. Failed lookups are tracked
    through the inherited
    :attr:`~rotifer.db.core.BaseCursor.missing` registry.

    Parameters
    ----------
    host : str, optional
        Host name of the ClickHouse server. Defaults to the ``host``
        configuration entry.
    port : int, optional
        Port of the server's HTTP interface. Defaults to the ``port``
        configuration entry.
    user : str, optional
        User name. Defaults to the ``user`` configuration entry.
    password : str, optional
        Password. Defaults to the ``password`` configuration entry.
    database : str, optional
        Name of the database. Defaults to the ``database``
        configuration entry.
    table : str, optional
        Name of the table. Defaults to the ``table`` configuration
        entry.
    secure : bool, default False
        Whether to connect over HTTPS.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    client : clickhouse_connect.driver.client.Client
        The connection, opened on first access.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.IdMappingCursor : identifier mapping cursor
    """
    def __init__(
            self,
            host = config['host'],
            port = config['port'],
            user = config['user'],
            password = config['password'],
            database = config['database'],
            table = config['table'],
            secure = False,
            progress = False,
            *args, **kwargs
        ):
        super().__init__(progress=progress, *args, **kwargs)
        self.host = host
        self.port = port
        self.user = user
        self.password = password
        self.database = database
        self.table = table
        self.secure = secure
        self._client = None

    @property
    def client(self):
        """
        The ClickHouse client, opened on first access.

        Returns
        -------
        clickhouse_connect.driver.client.Client

        Raises
        ------
        ImportError
            If the ``clickhouse_connect`` package is not installed.
        """
        if isinstance(self._client, types.NoneType):
            import clickhouse_connect
            self._client = clickhouse_connect.get_client(
                host = self.host,
                port = self.port,
                username = self.user,
                password = self.password,
                database = self.database,
                secure = self.secure,
            )
        return self._client

    @property
    def qualified_name(self):
        """
        The table name, qualified by its database.

        Returns
        -------
        str
            For example, ``uniprot.idmapping``.
        """
        return f'{self.database}.{self.table}'

    def query(self, sql, parameters=None):
        """
        Run a query and return its result as a dataframe.

        Parameters
        ----------
        sql : str
            A SELECT statement. Use ClickHouse's server side binding
            syntax, e.g. ``{ids:Array(String)}``, to refer to
            `parameters`.
        parameters : dict, optional
            Values bound to the placeholders in `sql`.

        Returns
        -------
        pandas.DataFrame
        """
        return self.client.query_df(sql, parameters=parameters)

    def command(self, sql, parameters=None):
        """
        Run a statement that does not return a result set.

        Parameters
        ----------
        sql : str
            Any statement, such as CREATE, ALTER or DROP.
        parameters : dict, optional
            Values bound to the placeholders in `sql`.

        Returns
        -------
        object
            Whatever the driver returns for the statement.
        """
        return self.client.command(sql, parameters=parameters)

    def has_table(self, name=None, database=None):
        """
        Find whether a table exists.

        Parameters
        ----------
        name : str, optional
            Table name. Defaults to the cursor's table.
        database : str, optional
            Database name. Defaults to the cursor's database.

        Returns
        -------
        bool
        """
        name = name or self.table
        database = database or self.database
        found = self.query(
            "SELECT count() AS n FROM system.tables WHERE database = {db:String} AND name = {tb:String}",
            parameters = {'db': database, 'tb': name},
        )
        return bool(found.n.iloc[0])

    @property
    def schema(self):
        """
        The SQL statement that created the cursor's table.

        Returns
        -------
        str
            Empty when the table does not exist.
        """
        if not self.has_table():
            return ""
        return self.client.command(f'SHOW CREATE TABLE {self.qualified_name}')

    def count(self):
        """
        Count the rows of the cursor's table.

        Returns
        -------
        int
            Zero when the table does not exist.
        """
        if not self.has_table():
            return 0
        return int(self.query(f'SELECT count() AS n FROM {self.qualified_name}').n.iloc[0])

    def _batches(self, targets, batch_size):
        """
        Split a set of identifiers into batches.

        Parameters
        ----------
        targets : iterable
            Identifiers to split.
        batch_size : int
            Maximum number of identifiers per batch.

        Yields
        ------
        list
            One batch of identifiers.
        """
        targets = list(targets)
        for start in range(0, len(targets), batch_size):
            yield targets[start:start + batch_size]

class BaseIdMappingCursor(rotifer.db.methods.IdMappingCursor, BaseClickHouseCursor):
    """
    Shared behaviour of the cursors reading the identifier mapping table.

    This class is not meant to be used directly: it builds, loads and
    describes the table that :class:`IdMappingCursor`,
    :class:`CrossReferenceCursor` and :class:`MappingCursor` query.
    Subclasses only have to set :attr:`column`, the name of the column
    their queries search.

    Parameters
    ----------
    id_type : str or list of str, optional
        Restrict results to these cross-referenced databases, e.g.
        ``RefSeq`` or ``['EMBL-CDS', 'GeneID']``. By default every
        database is reported.
    release : str, optional
        Restrict results to one UniProt release, e.g. ``2026_01``.
        Since the table is partitioned by release, setting this makes
        queries read a single partition. By default every release
        stored in the table is searched.
    batch_size : int, optional
        Number of identifiers sent to the server per query. Defaults
        to the ``batch_size`` configuration entry.
    **kwargs
        Connection parameters, passed to
        :class:`BaseClickHouseCursor`.

    See Also
    --------
    rotifer.db.uniprot.io.IdMappingCursor : the flat file this table is loaded from
    """

    #: Name of the column the cursor's queries search.
    column = 'accession'

    def __init__(
            self,
            id_type = None,
            release = config['release'],
            batch_size = config['batch_size'],
            *args, **kwargs
        ):
        super().__init__(*args, **kwargs)
        self.id_type = id_type
        self.release = release
        self.batch_size = batch_size
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

    def _filters(self, parameters):
        """
        Build the SQL conditions common to every cursor.

        Parameters
        ----------
        parameters : dict
            Query parameters, updated in place with the values bound
            by the conditions returned.

        Returns
        -------
        list of str
            Conditions to append to a WHERE clause.
        """
        conditions = []
        id_type = self._id_types()
        if id_type:
            conditions.append("id_type IN {id_type:Array(String)}")
            parameters['id_type'] = id_type
        if self.release:
            conditions.append("release = {release:String}")
            parameters['release'] = self.release
        return conditions

    def id_types(self):
        """
        List the cross-referenced databases present in the table.

        Returns
        -------
        pandas.DataFrame
            Columns ``id_type`` and ``rows``, sorted by decreasing
            number of rows.

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> ruch.IdMappingCursor().id_types()  # doctest: +SKIP
        """
        parameters = {}
        conditions = []
        if self.release:
            conditions.append("release = {release:String}")
            parameters['release'] = self.release
        where = f'WHERE {" AND ".join(conditions)}' if conditions else ""
        return self.query(
            f'SELECT id_type, count() AS rows FROM {self.qualified_name} {where} GROUP BY id_type ORDER BY rows DESC',
            parameters = parameters,
        )

    def create(self, replace=False, release=None):
        """
        Create the identifier mapping table.

        The schema is read from
        ``share/rotifer/db/uniprot/clickhouse/idmapping.sql`` and is
        located through :func:`rotifer.core.functions.findDataFiles`,
        so a copy under ``~/.rotifer/share`` takes precedence.

        Parameters
        ----------
        replace : bool, default False
            If True, drop any existing table before creating it.
        release : str, optional
            Value stored in the ``release`` column of rows inserted
            without one. Defaults to the cursor's ``release``.

        Returns
        -------
        bool
            Whether the table exists after the call.

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> ic = ruch.IdMappingCursor(release='2026_01')  # doctest: +SKIP
        >>> ic.create()  # doctest: +SKIP
        """
        sqlfile = rcf.findDataFiles(__name__ + ".idmapping.sql")
        if not sqlfile:
            logger.error("Could not find the SQL file describing the idmapping table")
            return False
        sql = open(sqlfile, "rt").read()
        sql = sql.format(
            database = self.database,
            table = self.table,
            release = release if not isinstance(release, types.NoneType) else self.release,
        )

        if replace:
            self.command(f'DROP TABLE IF EXISTS {self.qualified_name}')

        # Comments are stripped before the statements are split apart,
        # so that a semicolon inside a comment is not mistaken for the
        # end of a statement
        body = "\n".join([ x for x in sql.split("\n") if not x.strip().startswith("--") ])
        for statement in body.split(";"):
            if not statement.strip():
                continue
            self.command(statement)

        return self.has_table()

    def load(self, source, release=None, method='client', chunksize=config['chunksize'], executable=config['executable']):
        """
        Load a copy of ``idmapping.dat`` into the table.

        Parameters
        ----------
        source : str or rotifer.db.uniprot.io.IdMappingCursor
            The flat file to load: either its path, the root of a
            local UniProt mirror, or a cursor already pointing at one.
        release : str, optional
            Value stored in the ``release`` column of every row
            loaded. Defaults to the cursor's ``release``.
        method : str, default 'client'
            How to send the data:

            ``client``
                Pipe the file through the ``clickhouse client``
                command line program. This is by far the fastest
                option and the one to use for a full release, but it
                requires the program to be installed and able to
                reach the server.
            ``python``
                Read the file in chunks with
                :meth:`rotifer.db.uniprot.io.IdMappingCursor.reader`
                and insert each chunk through the driver. Slower, but
                it needs nothing besides this package and it honours
                the source cursor's ``id_type`` filter, which makes
                it convenient for loading a subset.

        chunksize : int, optional
            Rows per chunk when ``method='python'``. Defaults to the
            ``chunksize`` configuration entry.
        executable : str, optional
            Name or path of the ClickHouse program used when
            ``method='client'``. Defaults to the ``executable``
            configuration entry.

        Returns
        -------
        int
            Number of rows in the table after the load.

        Raises
        ------
        ValueError
            If `method` is not ``client`` or ``python``.

        Note
        ----
        Loading a full release moves a few billion rows and takes
        hours. The table is partitioned by release, so an interrupted
        load is cleaned up with
        ``ALTER TABLE ... DROP PARTITION '<release>'`` before trying
        again.

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> ic = ruch.IdMappingCursor(release='2026_01')  # doctest: +SKIP
        >>> ic.create()  # doctest: +SKIP
        >>> ic.load("/scratch/global/databases/uniprot")  # doctest: +SKIP
        """
        from rotifer.db.uniprot import io as ruio

        if isinstance(release, types.NoneType):
            release = self.release
        if isinstance(source, ruio.IdMappingCursor):
            reader = source
        else:
            reader = ruio.IdMappingCursor(path=source, progress=self.progress)
        if isinstance(reader.datafile, types.NoneType):
            logger.error(f'No idmapping file found for {source}')
            return self.count()

        if not self.has_table():
            self.create(release=release)

        if method == 'client':
            insert = (
                f'INSERT INTO {self.qualified_name} '
                f"SELECT c1, c2, c3, '{release}' "
                f"FROM input('c1 String, c2 String, c3 String') FORMAT TabSeparated"
            )
            command = [
                executable, "client",
                "--host", str(self.host),
                "--user", str(self.user),
                "--query", insert,
            ]
            if self.password:
                command += ["--password", str(self.password)]
            reading = f'zcat -f -- {reader.datafile}' if reader.compressed else f'cat -- {reader.datafile}'
            if self.progress:
                logger.warn(f'Loading {reader.datafile} into {self.qualified_name}, release {release}...')
            pipeline = subprocess.run(
                ["/bin/sh","-c", f'{reading} | ' + " ".join([ f"'{x}'" if " " in str(x) else str(x) for x in command ])],
                capture_output = True,
                text = True,
            )
            if pipeline.returncode != 0:
                logger.error(f'Failed to load {reader.datafile}: {pipeline.stderr}')

        elif method == 'python':
            if self.progress:
                logger.warn(f'Loading {reader.datafile} into {self.qualified_name} in chunks of {chunksize} rows...')
            for chunk in reader.reader(chunksize=chunksize):
                chunk = chunk.copy()
                chunk['release'] = release
                self.insert(chunk)

        else:
            raise ValueError(f'Unknown load method {method}: use "client" or "python"')

        return self.count()

    def insert(self, data, release=None):
        """
        Insert identifier mappings into the table.

        Parameters
        ----------
        data : pandas.DataFrame
            Rows to insert. Must have the columns listed in
            :attr:`~rotifer.db.methods.IdMappingCursor.columns`; a
            ``release`` column is added when missing.
        release : str, optional
            Value for the ``release`` column of rows that lack one.
            Defaults to the cursor's ``release``.
        """
        if data.empty:
            return
        data = data[self.columns].copy()
        data['release'] = release if not isinstance(release, types.NoneType) else self.release
        self.client.insert_df(table=self.table, df=data, database=self.database)

    def drop_release(self, release):
        """
        Delete every row of one UniProt release.

        Because the table is partitioned by release, this drops a
        whole partition and is nearly instantaneous.

        Parameters
        ----------
        release : str
            The release to remove, e.g. ``2024_06``.
        """
        self.command(f'ALTER TABLE {self.qualified_name} DROP PARTITION {{release:String}}', parameters={'release': release})

    def __getitem__(self, accessions):
        """
        Fetch identifier mappings, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to search in the column named by
            :attr:`column`.

        Returns
        -------
        pandas.DataFrame
            Mapping rows for the identifiers found, with the columns
            listed in
            :attr:`~rotifer.db.methods.IdMappingCursor.columns`.
            Identifiers that produced no row are registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()

        stack = []
        try:
            for batch in self._batches(targets, self.batch_size):
                parameters = {'targets': batch}
                conditions = [f'{self.column} IN {{targets:Array(String)}}'] + self._filters(parameters)
                stack.append(self.query(
                    f'SELECT accession, id_type, id FROM {self.qualified_name}'
                    f' WHERE {" AND ".join(conditions)}'
                    f' ORDER BY accession, id_type, id',
                    parameters = parameters,
                ))
        except Exception as error:
            # An unreachable or broken server must not abort the caller:
            # registering the query as missing lets a delegator hand it
            # to the next backend, and retry stays True because another
            # attempt may well succeed.
            logger.error(f'Query to {self.qualified_name} at {self.host} failed: {error}')
            self.update_missing(targets, error=f'ClickHouse query failed: {error}', retry=True)
            return self.empty()

        df = pd.concat(stack, ignore_index=True) if stack else self.empty()

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'Identifier not found in {self.qualified_name}', retry=False)

        return df

    def fetchone(self, accessions):
        """
        Iterate over identifier mappings, one batch at a time.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers to search.

        Yields
        ------
        pandas.DataFrame
            The mapping rows of one batch of at most ``batch_size``
            identifiers. Input order is not preserved.
        """
        targets = self.parse_ids(accessions)
        for batch in self._batches(targets, self.batch_size):
            found = self.__getitem__(batch)
            if not found.empty:
                yield found

class IdMappingCursor(BaseIdMappingCursor):
    """
    Fetch the cross-references of UniProtKB accessions.

    This is the forward lookup the table is sorted by, so it is the
    cheapest query available: ClickHouse reads only the granules that
    hold the requested accessions.

    Parameters
    ----------
    id_type : str or list of str, optional
        Restrict results to these cross-referenced databases.
    release : str, optional
        Restrict results to one UniProt release.
    batch_size : int, optional
        Number of accessions sent to the server per query.
    **kwargs
        Connection parameters, passed to
        :class:`BaseClickHouseCursor`.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.CrossReferenceCursor : the reverse lookup
    rotifer.db.uniprot.clickhouse.MappingCursor : translate between two databases
    rotifer.db.uniprot.io.IdMappingCursor : same data, read from the flat file

    Examples
    --------
    Every cross-reference of two accessions:

    >>> from rotifer.db.uniprot import clickhouse as ruch
    >>> ic = ruch.IdMappingCursor()  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP

    Only their RefSeq proteins:

    >>> ic = ruch.IdMappingCursor(id_type='RefSeq')  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP
    """
    column = 'accession'

class CrossReferenceCursor(BaseIdMappingCursor):
    """
    Fetch the UniProtKB accessions of identifiers from other databases.

    This is the reverse lookup, answered by the table's ``by_id``
    projection, a second copy of the data sorted by identifier that
    ClickHouse selects on its own. Naming the database in `id_type`
    makes the query cheaper still.

    Parameters
    ----------
    id_type : str or list of str, optional
        Restrict the search to these cross-referenced databases.
    release : str, optional
        Restrict results to one UniProt release.
    batch_size : int, optional
        Number of identifiers sent to the server per query.
    **kwargs
        Connection parameters, passed to
        :class:`BaseClickHouseCursor`.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.IdMappingCursor : the forward lookup

    Examples
    --------
    Which UniProtKB entry is a RefSeq protein part of?

    >>> from rotifer.db.uniprot import clickhouse as ruch
    >>> xc = ruch.CrossReferenceCursor(id_type='RefSeq')  # doctest: +SKIP
    >>> xc.fetchall(["YP_031579.1"])  # doctest: +SKIP
    """
    column = 'id'

class MappingCursor(BaseIdMappingCursor):
    """
    Translate identifiers from one database into another.

    This is the query UniProt's online ID mapping service answers:
    the identifiers of `from_type` are looked up in the table, and
    every identifier of `to_type` sharing their UniProtKB accession is
    returned. The two lookups are done in one server side join, so the
    intermediate accessions never travel over the network.

    Parameters
    ----------
    from_type : str
        Name of the database the queried identifiers belong to, as
        written in ``idmapping.dat``, e.g. ``EMBL-CDS``. Use
        ``UniProtKB-AC`` to start from UniProtKB accessions
        themselves.
    to_type : str
        Name of the database to translate into, e.g. ``RefSeq``. Use
        ``UniProtKB-AC`` to translate into UniProtKB accessions.
    release : str, optional
        Restrict results to one UniProt release.
    batch_size : int, optional
        Number of identifiers sent to the server per query.
    **kwargs
        Connection parameters, passed to
        :class:`BaseClickHouseCursor`.

    Attributes
    ----------
    columns : list of str
        ``['from', 'accession', 'to']``, overriding the three columns
        of the mapping table.

    See Also
    --------
    rotifer.db.uniprot.webapi.idmapping : the same query, run by UniProt's servers

    Examples
    --------
    Map GenBank CDS identifiers to RefSeq proteins:

    >>> from rotifer.db.uniprot import clickhouse as ruch
    >>> mc = ruch.MappingCursor(from_type='EMBL-CDS', to_type='RefSeq')  # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"])  # doctest: +SKIP
    """

    _columns = ['from','accession','to']
    column = 'from'

    def __init__(self, from_type, to_type, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.from_type = from_type
        self.to_type = to_type

    def __getitem__(self, accessions):
        """
        Translate identifiers, dictionary style.

        Parameters
        ----------
        accessions : str or iterable of str
            Identifiers of the database named by ``from_type``.

        Returns
        -------
        pandas.DataFrame
            Columns ``from``, ``accession`` and ``to``: the queried
            identifier, the UniProtKB accession that links it to the
            result, and the identifier in the database named by
            ``to_type``. Identifiers with no translation are
            registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()

        # UniProtKB accessions are the join key, not rows of the table
        source = "accession" if self.from_type == "UniProtKB-AC" else "id"
        target = "accession" if self.to_type == "UniProtKB-AC" else "id"

        stack = []
        try:
            stack = self._mapping_batches(targets, source, target)
        except Exception as error:
            logger.error(f'Query to {self.qualified_name} at {self.host} failed: {error}')
            self.update_missing(targets, error=f'ClickHouse query failed: {error}', retry=True)
            return self.empty()

        df = pd.concat(stack, ignore_index=True) if stack else self.empty()

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'No {self.to_type} identifier found for this {self.from_type} identifier', retry=False)

        return df

    def _mapping_batches(self, targets, source, target):
        """
        Run the join, one batch of identifiers at a time.

        Parameters
        ----------
        targets : set of str
            Identifiers to translate.
        source, target : str
            Names of the columns holding the queried and the returned
            identifiers, either ``accession`` or ``id``.

        Returns
        -------
        list of pandas.DataFrame
            One dataframe per batch.
        """
        stack = []
        for batch in self._batches(targets, self.batch_size):
            parameters = {'targets': batch}
            left = [f'f.{source} IN {{targets:Array(String)}}']
            if source == "id":
                left.append("f.id_type = {from_type:String}")
                parameters['from_type'] = self.from_type
            right = []
            if target == "id":
                right.append("t.id_type = {to_type:String}")
                parameters['to_type'] = self.to_type
            if self.release:
                left.append("f.release = {release:String}")
                right.append("t.release = {release:String}")
                parameters['release'] = self.release
            where = " AND ".join(left + right)
            stack.append(self.query(
                f'SELECT DISTINCT f.{source} AS `from`, f.accession AS accession, t.{target} AS `to`'
                f' FROM {self.qualified_name} AS f'
                f' INNER JOIN {self.qualified_name} AS t ON f.accession = t.accession'
                f' WHERE {where}'
                f' ORDER BY `from`, accession, `to`',
                parameters = parameters,
            ))
        return stack

if __name__ == '__main__':
    pass
