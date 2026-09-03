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
:class:`rotifer.db.uniprot.mirror.IdMappingCursor`. Its schema lives in
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
import uuid
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
    'batch_size': 5000,
    'submit_threshold': 1000,
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
        # Each cursor owns a differently named temporary table, so that
        # cursors sharing a session cannot overwrite each other's query
        self._query_table = '_rotifer_query_' + uuid.uuid4().hex[:12]

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
            # The session is deliberately not bound to self.database:
            # the driver refuses to connect at all when the database
            # does not exist yet, which would make it impossible to
            # create one. Every statement here names its table in
            # full, so the session database is never consulted.
            self._client = clickhouse_connect.get_client(
                host = self.host,
                port = self.port,
                username = self.user,
                password = self.password,
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
            A SELECT statement. Queries that carry a list of
            identifiers must use client side binding, ``%(name)s``,
            so that the values travel in the request body: server
            side binding puts them in the URL, which the server
            rejects as "Field value too long" beyond roughly 6000
            identifiers.
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

    def submit(self, accessions):
        """
        Send a list of identifiers to a temporary table.

        Queries then refer to that table instead of listing the
        identifiers in the SQL text, which is what makes an
        arbitrarily long query possible: values written into the
        statement are bounded by ClickHouse's ``max_query_size``, 256
        KiB by default, or about 18000 accessions.

        The table is temporary and lives in the connection's session,
        so it disappears when the connection does and is invisible to
        every other client.

        Parameters
        ----------
        accessions : iterable of str
            The identifiers to make available to the next query.

        Returns
        -------
        str
            Name of the table holding them.

        Note
        ----
        A ClickHouse session serves one query at a time, so a cursor
        that has submitted a list must not be shared between threads
        until the query using it has finished.

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> c = ruch.IdMappingCursor()  # doctest: +SKIP
        >>> table = c.submit(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP
        >>> c.query(f"SELECT count() FROM {table}")  # doctest: +SKIP
        """
        import pandas as pd
        self.cleanup()
        self.command(f'CREATE TEMPORARY TABLE {self._query_table} (id String) ENGINE = Memory')
        ids = pd.DataFrame({'id': [ str(x) for x in accessions ]})
        self.client.insert_df(table=self._query_table, df=ids)
        return self._query_table

    def cleanup(self):
        """
        Drop this cursor's temporary table of identifiers.

        Safe to call when nothing was submitted.
        """
        self.command(f'DROP TEMPORARY TABLE IF EXISTS {self._query_table}')

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
        Number of identifiers written into one query, for queries
        small enough to list them. Defaults to the ``batch_size``
        configuration entry.
    submit_threshold : int, optional
        Queries carrying at least this many identifiers send them to a
        temporary table with :meth:`~BaseClickHouseCursor.submit`
        instead of listing them in the SQL, which removes the limit on
        how many can be asked for at once. Below it, listing them is
        one round trip instead of three and therefore quicker.
        Defaults to the ``submit_threshold`` configuration entry.
    **kwargs
        Connection parameters, passed to
        :class:`BaseClickHouseCursor`.

    See Also
    --------
    rotifer.db.uniprot.mirror.IdMappingCursor : the flat file this table is loaded from
    """

    #: Name of the column the cursor's queries search.
    column = 'accession'

    def __init__(
            self,
            id_type = None,
            release = config['release'],
            batch_size = config['batch_size'],
            submit_threshold = config['submit_threshold'],
            *args, **kwargs
        ):
        super().__init__(*args, **kwargs)
        self.id_type = id_type
        self.release = release
        self.batch_size = batch_size
        self.submit_threshold = submit_threshold
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
            conditions.append("id_type IN %(id_type)s")
            parameters['id_type'] = tuple(id_type)
        if self.release:
            conditions.append("release = %(release)s")
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
            conditions.append("release = %(release)s")
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

    def load(self, source, release=None, method='auto', chunksize=config['chunksize'], executable=config['executable']):
        """
        Load a copy of ``idmapping.dat`` into the table.

        Parameters
        ----------
        source : str or rotifer.db.uniprot.mirror.IdMappingCursor
            The flat file to load: either its path, the root of a
            local UniProt mirror, or a cursor already pointing at one.
        release : str, optional
            Value stored in the ``release`` column of every row
            loaded. Defaults to the cursor's ``release``.
        method : str, default 'auto'
            How to send the data:

            ``auto``
                Use ``client`` when the ClickHouse program is on the
                PATH, and ``python`` otherwise.
            ``client``
                Pipe the file through the ``clickhouse client``
                command line program. This is by far the fastest
                option and the one to use for a full release, but it
                requires the program to be installed and able to
                reach the server.
            ``python``
                Read the file in chunks with
                :meth:`rotifer.db.uniprot.mirror.IdMappingCursor.reader`
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
            If `method` is not ``auto``, ``client`` or ``python``.

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
        from rotifer.db.uniprot import mirror as rum

        if isinstance(release, types.NoneType):
            release = self.release
        if isinstance(source, rum.IdMappingCursor):
            reader = source
        else:
            reader = rum.IdMappingCursor(path=source, progress=self.progress)
        if isinstance(reader.datafile, types.NoneType):
            logger.error(f'No idmapping file found for {source}')
            return self.count()

        if not self.has_table():
            self.create(release=release)

        if method == 'auto':
            import shutil
            method = 'client' if shutil.which(executable) else 'python'
            if self.progress:
                logger.warn(f'Loading with method={method}')

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
            raise ValueError(f'Unknown load method {method}: use "auto", "client" or "python"')

        return self.count()

    def is_empty(self, release=None):
        """
        Find whether the table holds no rows for a release.

        Parameters
        ----------
        release : str, optional
            The release to look for. Defaults to the cursor's
            ``release``; when neither is set, the whole table is
            considered.

        Returns
        -------
        bool
            True when the table does not exist, or holds no matching
            row.
        """
        if not self.has_table():
            return True
        release = release if not isinstance(release, types.NoneType) else self.release
        if release:
            found = self.query(
                f'SELECT count() AS n FROM {self.qualified_name} WHERE release = {{release:String}}',
                parameters = {'release': release},
            )
            return not int(found.n.iloc[0])
        return not self.count()

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
        # release may have been cleared to None to widen queries, but
        # the column is a String and never takes None
        data['release'] = (release if not isinstance(release, types.NoneType) else self.release) or ''
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
            if len(targets) >= self.submit_threshold:
                # Too many to write into the statement: hand them over
                # as a table and let the query refer to it
                table = self.submit(targets)
                try:
                    parameters = {}
                    conditions = [f'{self.column} IN (SELECT id FROM {table})'] + self._filters(parameters)
                    stack.append(self.query(
                        f'SELECT accession, id_type, id FROM {self.qualified_name}'
                        f' WHERE {" AND ".join(conditions)}'
                        f' ORDER BY accession, id_type, id',
                        parameters = parameters,
                    ))
                finally:
                    self.cleanup()
            else:
                for batch in self._batches(targets, self.batch_size):
                    parameters = {'targets': tuple(batch)}
                    conditions = [f'{self.column} IN %(targets)s'] + self._filters(parameters)
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
    rotifer.db.uniprot.mirror.IdMappingCursor : same data, read from the flat file

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
        if len(targets) >= self.submit_threshold:
            # Too many to write into the statement: hand them over as a
            # table and let the join refer to it
            table = self.submit(targets)
            try:
                return [ self._mapping_query(f'f.{source} IN (SELECT id FROM {table})', {}, source, target) ]
            finally:
                self.cleanup()

        stack = []
        for batch in self._batches(targets, self.batch_size):
            parameters = {'targets': tuple(batch)}
            stack.append(self._mapping_query(f'f.{source} IN %(targets)s', parameters, source, target))
        return stack

    def _mapping_query(self, restriction, parameters, source, target):
        """
        Run the join for one set of queried identifiers.

        Parameters
        ----------
        restriction : str
            The SQL condition selecting them, either an inline list or
            a reference to the table filled by
            :meth:`~BaseClickHouseCursor.submit`.
        parameters : dict
            Query parameters, extended in place with the values bound
            by the conditions this method adds.
        source, target : str
            Names of the columns holding the queried and the returned
            identifiers, either ``accession`` or ``id``.

        Returns
        -------
        pandas.DataFrame
            Columns ``from``, ``accession`` and ``to``.
        """
        left = [restriction]
        if source == "id":
            left.append("f.id_type = %(from_type)s")
            parameters['from_type'] = self.from_type
        right = []
        if target == "id":
            right.append("t.id_type = %(to_type)s")
            parameters['to_type'] = self.to_type
        if self.release:
            left.append("f.release = %(release)s")
            right.append("t.release = %(release)s")
            parameters['release'] = self.release
        where = " AND ".join(left + right)
        return self.query(
            f'SELECT DISTINCT f.{source} AS `from`, f.accession AS accession, t.{target} AS `to`'
            f' FROM {self.qualified_name} AS f'
            f' INNER JOIN {self.qualified_name} AS t ON f.accession = t.accession'
            f' WHERE {where}'
            f' ORDER BY `from`, accession, `to`',
            parameters = parameters,
        )

if __name__ == '__main__':
    pass
