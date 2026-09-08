__doc__ = """
Shared machinery for cursors backed by ClickHouse.

:class:`BaseClickHouseCursor` is not meant to be used directly. It
holds everything a ClickHouse cursor needs that does not depend on
what the data means: the connection, the statement helpers, table
introspection and creation, bulk loading, and the temporary table used
to query against long lists of identifiers.

Concrete cursors combine it with one of the mixins in
:mod:`rotifer.db.methods`, which decide what the returned dataframe
looks like. :mod:`rotifer.db.uniprot.clickhouse` is the worked
example.
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
from rotifer.db.sql.clickhouse import config
from rotifer.core import functions as rcf
logger = rotifer.logging.getLogger(__name__)

class BaseClickHouseCursor(rotifer.db.core.BaseCursor):
    """
    Connection, statement and table helpers for ClickHouse cursors.

    This class is not meant to be used directly: it is the parent of
    every ClickHouse backed cursor. Subclasses add the queries that
    know what the rows mean, and may set :attr:`_schema_resource` so
    that :meth:`create` can build their table.

    Parameters
    ----------
    host : str, optional
        Host name of the ClickHouse server. Defaults to the ``host``
        entry of the :mod:`rotifer.db.sql.clickhouse` configuration.
    port : int, optional
        Port of the server's HTTP interface.
    user : str, optional
        User name.
    password : str, optional
        Password.
    dbname : str, optional
        Name of the database holding the cursor's table. This was
        called ``database`` until it was renamed for clarity; passing
        the old name now raises rather than being ignored.
    table : str, optional
        Name of the cursor's table.
    secure : bool, optional
        Whether to connect over HTTPS.
    batch_size : int, optional
        Number of identifiers written into one query by
        :meth:`_batches`.
    submit_threshold : int, optional
        Number of identifiers from which a query should send them to a
        temporary table with :meth:`submit` rather than list them in
        the statement.
    progress : bool, default False
        Whether to print progress messages.

    Attributes
    ----------
    client : clickhouse_connect.driver.client.Client
        The connection, opened on first access.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.BaseMappingCursor : a cursor built on this class
    """

    #: Resource name of the SQL file describing this cursor's table,
    #: as accepted by :func:`rotifer.core.functions.findDataFiles`.
    #: Subclasses that can create their own table set this.
    _schema_resource = None

    def __init__(
            self,
            host = config['host'],
            port = config['port'],
            user = config['user'],
            password = config['password'],
            dbname = config['dbname'],
            table = config['table'],
            secure = config['secure'],
            batch_size = config['batch_size'],
            submit_threshold = config['submit_threshold'],
            progress = False,
            *args, **kwargs
        ):
        # BaseCursor accepts and ignores unknown keywords, so a caller
        # still passing the old name would silently get the default
        # database instead of the one they asked for
        if 'database' in kwargs:
            raise TypeError(
                "the 'database' parameter is now called 'dbname'; "
                f"pass dbname={kwargs['database']!r} instead"
            )
        super().__init__(progress=progress, *args, **kwargs)
        self.host = host
        self.port = port
        self.user = user
        self.password = password
        self.dbname = dbname
        self.table = table
        self.secure = secure
        self.batch_size = batch_size
        self.submit_threshold = submit_threshold
        self._client = None
        # Each cursor owns a differently named temporary table, so that
        # cursors sharing a session cannot overwrite each other's query
        self._query_table = '_rotifer_query_' + uuid.uuid4().hex[:12]

    # Connection

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
            # The session is deliberately not bound to self.dbname:
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
        The cursor's table, qualified by its database.

        Returns
        -------
        str
            For example, ``rotifer.idmapping``.
        """
        return f'{self.dbname}.{self.table}'

    # Statements

    def query(self, sql, parameters=None):
        """
        Run a query and return its result as a dataframe.

        Parameters
        ----------
        sql : str
            A SELECT statement. Queries that carry a list of values
            must use client side binding, ``%(name)s``, so that the
            values travel in the request body: server side binding
            puts them in the URL, which the server rejects as "Field
            value too long" beyond roughly 6000 values. Longer lists
            belong in :meth:`submit`.
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

    # Introspection

    #: One registry per database, shared by every table in it, naming
    #: the file each table was loaded from. Not per table: provenance
    #: is the same kind of fact whatever was loaded, and a registry
    #: per table would multiply bookkeeping tables alongside the data.
    _sources_table = 'rotifer_sources'

    @property
    def sources_table(self):
        """
        Qualified name of this database's load registry.

        Returns
        -------
        str
        """
        return f'{self.dbname}.{self._sources_table}'

    @property
    def source_version(self):
        """
        Version label distinguishing loads into the same table.

        Tables that hold several versions at once, e.g. one partition
        per release, override this so each is recorded separately. The
        default is a single unversioned body of data.

        Returns
        -------
        str
        """
        return ''

    def create_sources(self):
        """
        Create this database's load registry, if it is missing.

        A row is written only by a load that ran to completion, which
        is what lets :meth:`content_id` distinguish a table holding a
        whole copy of a file from one left half filled by a load that
        was interrupted.
        """
        self.command(f"""
            CREATE TABLE IF NOT EXISTS {self.sources_table} (
                `table`   String,
                version   String,
                source    String,
                size      UInt64,
                mtime     Int64,
                rows      UInt64,
                checksum  String,
                loaded_at DateTime DEFAULT now()
            ) ENGINE = ReplacingMergeTree(loaded_at)
            ORDER BY (`table`, version, source)
        """)

    def record_source(self, datafile, rows, version=None, checksum=''):
        """
        Record that this table was loaded, whole, from a file.

        Call this only where a load is known to have finished: the
        presence of a row is the claim of completeness.

        Parameters
        ----------
        datafile : str
            Path of the file the rows came from. Recorded in full
            and with symbolic links resolved, so that the same file
            reached by two names is recorded as one file.
        rows : int
            Number of rows in the table afterwards.
        version : str, optional
            Version label. Defaults to :attr:`source_version`.
        checksum : str, optional
            Checksum of the file, when one was computed. Recorded for
            verification; routine comparisons use the cheap identity
            below, which costs a stat() rather than a full read.

        Returns
        -------
        bool
            Whether the record was written.
        """
        import os
        try:
            info = os.stat(datafile)
        except OSError:
            logger.warning(f'Cannot stat {datafile}: load not recorded')
            return False
        if isinstance(version, types.NoneType):
            version = self.source_version
        self.create_sources()
        self.client.insert(
            table = self._sources_table,
            database = self.dbname,
            column_names = ['table','version','source','size','mtime','rows','checksum'],
            data = [[str(self.table), str(version), os.path.realpath(datafile),
                     int(info.st_size), int(info.st_mtime), int(rows), str(checksum)]],
        )
        return True

    def forget_source(self, version=None, table=None):
        """
        Withdraw the record of a load.

        A row in the registry is the claim that this table holds a
        whole copy of a file. Whatever removes the rows has to remove
        the claim too, or the table goes on being taken for a copy of
        something it no longer has.

        Parameters
        ----------
        version : str, optional
            Version whose record to drop. Defaults to
            :attr:`source_version`; pass an empty string to drop the
            records of every version of this table.
        table : str, optional
            Table the record is about. Defaults to this cursor's.

        Returns
        -------
        bool
            Whether anything could be removed.
        """
        if isinstance(version, types.NoneType):
            version = self.source_version
        try:
            if not self.has_table(self._sources_table):
                return False
            sql = f'ALTER TABLE {self.sources_table} DELETE WHERE `table` = %(table)s'
            parameters = {'table': str(table or self.table)}
            if version:
                sql += ' AND version = %(version)s'
                parameters['version'] = str(version)
            # A mutation is asynchronous by default, so the statement
            # returns before the row is gone and whatever asked for the
            # claim to be withdrawn would still be shown it.
            sql += ' SETTINGS mutations_sync = 1'
            self.command(sql, parameters=parameters)
        except Exception:
            logger.debug('Could not remove the load record', exc_info=1)
            return False
        return True

    def content_id(self):
        """
        Identify the data this table holds.

        Answers only for a table loaded in full, and then with the
        identity of the file it was loaded from, so that a cursor
        reading that same file directly reports the same value and
        either can stand in for the other. A table never recorded,
        because its load was interrupted or predates this bookkeeping,
        yields None and nothing is skipped on its account.

        Returns
        -------
        str or None
            ``<path>:<size>:<mtime>``, or None.

        See Also
        --------
        rotifer.db.core.BaseCursor.content_id : what the value means
        record_source : how the value gets written
        """
        try:
            if not self.has_table(self._sources_table):
                return None
            sql = (f'SELECT source, size, mtime FROM {self.sources_table} '
                   f'WHERE `table` = %(table)s')
            parameters = {'table': self.table}
            version = self.source_version
            if version:
                sql += ' AND version = %(version)s'
                parameters['version'] = version
            frame = self.query(sql + ' ORDER BY loaded_at DESC', parameters=parameters)
        except Exception:
            logger.debug('Could not read the load registry', exc_info=1)
            return None
        if frame is None or frame.empty:
            return None
        # Without a version filter the cursor reads every version at
        # once, so it only matches a source that is the only one there.
        if not self.source_version and len(frame.drop_duplicates()) > 1:
            return None
        row = frame.iloc[0]
        return f'{row["source"]}:{int(row["size"])}:{int(row["mtime"])}'

    def has_table(self, name=None, dbname=None):
        """
        Find whether a table exists.

        Parameters
        ----------
        name : str, optional
            Table name. Defaults to the cursor's table.
        dbname : str, optional
            Database name. Defaults to the cursor's own.

        Returns
        -------
        bool
        """
        name = name or self.table
        dbname = dbname or self.dbname
        found = self.query(
            "SELECT count() AS n FROM system.tables WHERE database = {db:String} AND name = {tb:String}",
            parameters = {'db': dbname, 'tb': name},
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

    def count(self, where=None, parameters=None):
        """
        Count the rows of the cursor's table.

        Parameters
        ----------
        where : str, optional
            A condition restricting what is counted, without the
            ``WHERE`` keyword.
        parameters : dict, optional
            Values bound to the placeholders in `where`.

        Returns
        -------
        int
            Zero when the table does not exist.
        """
        if not self.has_table():
            return 0
        clause = f' WHERE {where}' if where else ''
        return int(self.query(
            f'SELECT count() AS n FROM {self.qualified_name}{clause}',
            parameters = parameters,
        ).n.iloc[0])

    def is_empty(self, where=None, parameters=None):
        """
        Find whether the table holds no rows.

        Parameters
        ----------
        where : str, optional
            A condition restricting what is looked for, without the
            ``WHERE`` keyword.
        parameters : dict, optional
            Values bound to the placeholders in `where`.

        Returns
        -------
        bool
            True when the table does not exist, or holds no matching
            row.
        """
        if not self.has_table():
            return True
        return not self.count(where=where, parameters=parameters)

    # Long queries

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
        >>> c = ruch.MappingCursor()  # doctest: +SKIP
        >>> table = c.submit(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP
        >>> c.query(f"SELECT count() FROM {table}")  # doctest: +SKIP
        """
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

    def _batches(self, targets, batch_size=None):
        """
        Split a set of identifiers into batches.

        Parameters
        ----------
        targets : iterable
            Identifiers to split.
        batch_size : int, optional
            Maximum number of identifiers per batch. Defaults to the
            cursor's ``batch_size``.

        Yields
        ------
        list
            One batch of identifiers.
        """
        batch_size = batch_size or self.batch_size
        targets = list(targets)
        for start in range(0, len(targets), batch_size):
            yield targets[start:start + batch_size]

    # Tables and data

    def create(self, replace=False, schema=None, **parameters):
        """
        Create the cursor's database and table.

        The statements are read from a SQL file located through
        :func:`rotifer.core.functions.findDataFiles`, so a copy under
        ``~/.rotifer/share`` takes precedence over the one shipped
        with rotifer. ``{dbname}`` and ``{table}`` in that file are
        filled in from the cursor; any other placeholder must be given
        in `parameters`.

        Parameters
        ----------
        replace : bool, default False
            If True, drop any existing table before creating it. Every
            row it holds is lost.
        schema : str, optional
            Resource name of the SQL file. Defaults to the subclass's
            :attr:`_schema_resource`.
        **parameters
            Further values substituted into the file.

        Returns
        -------
        bool
            Whether the table exists after the call.

        Raises
        ------
        ValueError
            If no schema resource is given and the subclass declares
            none.
        """
        schema = schema or self._schema_resource
        if isinstance(schema, types.NoneType):
            raise ValueError(f'{self.__name__} declares no _schema_resource and none was given')
        sqlfile = rcf.findDataFiles(schema)
        if not sqlfile:
            logger.error(f'Could not find the SQL file {schema}')
            return False
        sql = open(sqlfile, "rt").read().format(
            dbname = self.dbname,
            table = self.table,
            # accepted too, so that schema files written against the
            # older placeholder keep working
            database = self.dbname,
            **parameters,
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

    def insert(self, data):
        """
        Insert a dataframe into the cursor's table.

        Parameters
        ----------
        data : pandas.DataFrame
            Rows to insert. Its columns must match the table's.
        """
        if data.empty:
            return
        self.client.insert_df(table=self.table, df=data, database=self.dbname)

    def drop_partition(self, partition):
        """
        Delete every row of one partition.

        Dropping a partition is nearly instantaneous, which is what
        makes a partitioned table cheap to refresh: an obsolete slice
        goes in one statement instead of a mutation over every part.

        Parameters
        ----------
        partition : str
            The partition to remove.
        """
        self.command(
            f'ALTER TABLE {self.qualified_name} DROP PARTITION {{partition:String}}',
            parameters = {'partition': partition},
        )

    def load_file(self, path, select, columns, compressed=False, executable=None):
        """
        Stream a delimited file into the table with the ClickHouse client.

        The file is piped through the ``clickhouse client`` program,
        which is much faster than sending the rows through this
        driver, and is the only sane way to load a file of tens of
        gigabytes. It needs that program on the PATH and able to reach
        the server.

        Parameters
        ----------
        path : str
            The file to load.
        select : str
            The SELECT list turning the file's columns into the
            table's, e.g. ``"c1, c2, c3, '2026_01'"``.
        columns : str
            The file's columns, as accepted by ClickHouse's ``input``
            table function, e.g. ``"c1 String, c2 String, c3 String"``.
        compressed : bool, default False
            Whether `path` is gzip compressed.
        executable : str, optional
            Name or path of the ClickHouse program. Defaults to the
            ``executable`` entry of the
            :mod:`rotifer.db.sql.clickhouse` configuration.

        Returns
        -------
        bool
            Whether the load succeeded.
        """
        executable = executable or config['executable']
        insert = (
            f'INSERT INTO {self.qualified_name} '
            f'SELECT {select} '
            f"FROM input('{columns}') FORMAT TabSeparated"
        )
        command = [
            executable, "client",
            "--host", str(self.host),
            "--user", str(self.user),
            "--query", insert,
        ]
        if self.password:
            command += ["--password", str(self.password)]
        if self.progress:
            logger.warning(f'Loading {path} into {self.qualified_name}...')

        # No shell. The statement carries quotes of its own, around the
        # release and around the input() column list, so any attempt to
        # quote it into a command line ends with those closing the
        # quoting rather than being part of it, and the server is sent
        # something that will not parse.
        reader = None
        try:
            if compressed:
                reader = subprocess.Popen(['zcat', '-f', '--', path], stdout=subprocess.PIPE)
                stream = reader.stdout
            else:
                stream = open(path, 'rb')
            try:
                pipeline = subprocess.run(
                    command,
                    stdin = stream,
                    capture_output = True,
                    text = True,
                )
            finally:
                stream.close()
                if reader is not None:
                    reader.wait()
        except OSError as error:
            logger.error(f'Failed to load {path}: {error}')
            return False
        if reader is not None and reader.returncode:
            logger.error(f'Failed to read {path}: zcat exited {reader.returncode}')
            return False
        if pipeline.returncode != 0:
            logger.error(f'Failed to load {path}: {pipeline.stderr}')
            return False
        return True

if __name__ == '__main__':
    pass
