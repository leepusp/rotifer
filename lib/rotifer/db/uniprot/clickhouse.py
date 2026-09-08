__doc__ = """
Query UniProt identifier mappings stored in ClickHouse.

UniProt's ``idmapping.dat`` relates every UniProtKB accession to the
identifier of the same protein in each database UniProt
cross-references. The file is far too large to be searched
interactively (about 90 GB and a few billion rows in the 2026_01
release), so this module keeps it in a ClickHouse table and provides
cursors for the queries that table was designed to answer:

:class:`MappingCursor`
    Translate identifiers in any direction: from UniProtKB accessions
    to the databases they cross-reference, from those databases back
    to the accessions, or between two of them. Which direction is
    asked for is decided per query by its ``source`` and ``target``,
    not by picking a different class.

Everything here that is not about identifier mappings lives in
:mod:`rotifer.db.sql.clickhouse.core`, which any other ClickHouse backed
cursor can build on.

The table itself is created and populated through
:meth:`BaseMappingCursor.create` and
:meth:`BaseMappingCursor.load`, from a local copy of the flat file
read by :class:`rotifer.db.uniprot.mirror.MappingCursor`. Its schema
lives in ``share/rotifer/db/uniprot/clickhouse/idmapping.sql``.

Configuration
-------------
Connection parameters are read from ``~/.rotifer/etc/db/uniprot/clickhouse.yml``
and fall back to the shared :mod:`rotifer.db.sql.clickhouse` defaults, so
a server can be named once for every cursor or separately here.
"""

# Dependencies
import os
import types
import typing
import pandas as pd

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.sql.clickhouse.core
import rotifer.db.sql.mapping as sqlmap
import rotifer.db.sql.progress as sqlprog
from rotifer.db.sql.clickhouse import config as clickhouse_config
from rotifer.core import functions as rcf
logger = rotifer.logging.getLogger(__name__)

#: Kept so that ``from rotifer.db.uniprot.clickhouse import
#: BaseClickHouseCursor`` still works; the class itself now lives in
#: :mod:`rotifer.db.sql.clickhouse.core`.
BaseClickHouseCursor = rotifer.db.sql.clickhouse.core.BaseClickHouseCursor

# Defaults: the shared connection settings, with what UniProt adds
_defaults = dict(clickhouse_config)
_defaults.update({
    'dbname': 'rotifer',
    'release': '',
    'chunksize': 5000000,
})
config = rcf.loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

class BaseMappingCursor(rotifer.db.methods.MappingCursor, BaseClickHouseCursor):
    """
    Shared behaviour of the cursors reading the identifier mapping table.

    This class is not meant to be used directly: it builds, loads and
    describes the table that :class:`MappingCursor` queries.

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
    rotifer.db.uniprot.mirror.MappingCursor : the flat file this table is loaded from
    """

    #: Columns of the mapping table itself, which are not the columns
    #: this cursor returns: the table stores one row per
    #: cross-reference, while a mapping names both of its ends.
    _table_columns = ['accession','id_type','id']

    #: The one table these cursors read, and the role it plays.
    tables = {'mapping': 'idmapping'}

    #: Where :meth:`create` reads this table's definition from.
    _schema_resource = __name__ + ".idmapping.sql"

    def __init__(
            self,
            id_type = None,
            release = config['release'],
            *args, **kwargs
        ):
        # Connection settings default to this module's configuration
        # rather than the shared one, so that a UniProt server can be
        # named separately from every other ClickHouse table
        for key in ('host','port','user','password','dbname',
                    'secure','batch_size','submit_threshold'):
            kwargs.setdefault(key, config[key])
        super().__init__(*args, **kwargs)
        self.id_type = id_type
        self.release = release
        self.maxgetitem = 1000000
        self._databases = None

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

    def databases(self):
        """
        Name the cross-referenced databases this table holds.

        Read once and remembered: the query is a scan of a
        low cardinality column, which costs seconds on a table of a
        few billion rows and would otherwise be paid on every query.

        Returns
        -------
        set of str or None
            None when the table cannot be reached, so that an
            unreachable server narrows nothing.
        """
        if not isinstance(getattr(self, '_databases', None), types.NoneType):
            return self._databases
        try:
            frame = self.query(f'SELECT DISTINCT id_type FROM {self.qualified_name}')
        except Exception:
            logger.debug(f'Could not list the id_types of {self.qualified_name}', exc_info=1)
            return None
        self._databases = set(frame.id_type.dropna().astype(str))
        return self._databases

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
        >>> ruch.MappingCursor().id_types()  # doctest: +SKIP
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

        See Also
        --------
        rotifer.db.sql.clickhouse.core.BaseClickHouseCursor.create : the generic form

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> ic = ruch.MappingCursor(release='2026_01')  # doctest: +SKIP
        >>> ic.create()  # doctest: +SKIP
        """
        return super().create(
            replace = replace,
            release = release if not isinstance(release, types.NoneType) else self.release,
        )

    def load(self, mirror, release=None, method='auto', chunksize=config['chunksize'], executable=config['executable']):
        """
        Load a copy of ``idmapping.dat`` into the table.

        Parameters
        ----------
        mirror : str or rotifer.db.uniprot.mirror.MappingCursor
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
                :meth:`rotifer.db.uniprot.mirror.MappingCursor.reader`
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
        >>> ic = ruch.MappingCursor(release='2026_01')  # doctest: +SKIP
        >>> ic.create()  # doctest: +SKIP
        >>> ic.load("/scratch/global/databases/uniprot")  # doctest: +SKIP
        """
        from rotifer.db.uniprot import mirror as rum

        if isinstance(release, types.NoneType):
            release = self.release
        if isinstance(mirror, rum.MappingCursor):
            reader = mirror
        else:
            reader = rum.MappingCursor(path=mirror, progress=self.progress)
        if isinstance(reader.datafile, types.NoneType):
            logger.error(f'No idmapping file found for {mirror}')
            return self.count()

        if not self.has_table():
            self.create(release=release)

        if method == 'auto':
            import shutil
            method = 'client' if shutil.which(executable) else 'python'
            if self.progress:
                logger.warning(f'Loading with method={method}')

        estimate = sqlprog.estimate_rows(reader.datafile, compressed=reader.compressed)
        if method == 'client':
            # The rows go straight from the file into the server, so
            # there is nothing here to count: ask the table instead.
            with sqlprog.Watcher(lambda: self.count(), total=estimate,
                                 desc='loading', enabled=self.progress):
                self.load_file(
                    reader.datafile,
                    select = f"c1, c2, c3, '{release}'",
                    columns = 'c1 String, c2 String, c3 String',
                    compressed = reader.compressed,
                    executable = executable,
                )

        elif method == 'python':
            if self.progress:
                logger.warning(f'Loading {reader.datafile} into {self.qualified_name} in chunks of {chunksize} rows...')
            with sqlprog.Progress(total=estimate, desc='loading',
                                  enabled=self.progress) as bar:
                for chunk in reader.reader(chunksize=chunksize):
                    # insert() keeps only the table's own columns and
                    # stamps the release itself, so it has to be told
                    # which one: setting it on the frame would be
                    # dropped and the cursor's own release used.
                    self.insert(chunk, release=release)
                    bar.update(len(chunk))

        else:
            raise ValueError(f'Unknown load method {method}: use "auto", "client" or "python"')

        # Reached only when the load ran to the end, which is what
        # makes the record meaningful: an interrupted load raises or
        # dies before this and leaves the release unrecorded, so
        # nothing will take the half filled partition for a complete
        # copy of the file.
        rows = self.count()
        self.record_source(reader.datafile, rows, version=release)
        return rows

    @property
    def source_version(self):
        """
        Release distinguishing one load of this table from another.

        The table is partitioned by release, so each release is a
        separate body of data and is recorded separately.

        Returns
        -------
        str
        """
        return self.release or ''

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
        release = release if not isinstance(release, types.NoneType) else self.release
        if release:
            return super().is_empty(where="release = %(release)s", parameters={'release': release})
        return super().is_empty()

    def insert(self, data, release=None):
        """
        Insert identifier mappings into the table.

        Parameters
        ----------
        data : pandas.DataFrame
            Rows to insert. Must have the columns listed in
            :attr:`_table_columns`; a
            ``release`` column is added when missing.
        release : str, optional
            Value for the ``release`` column of rows that lack one.
            Defaults to the cursor's ``release``.
        """
        if data.empty:
            return
        data = data[self._table_columns].copy()
        # release may have been cleared to None to widen queries, but
        # the column is a String and never takes None
        data['release'] = (release if not isinstance(release, types.NoneType) else self.release) or ''
        super().insert(data)

    def drop_release(self, release):
        """
        Delete every row of one UniProt release.

        Because the table is partitioned by release, this drops a
        whole partition and is nearly instantaneous.

        The record of where that release was loaded from goes with
        it. It is what tells a delegator this table holds a copy of
        that file, and keeping it would have the mirror skipped in
        favour of a partition that is no longer there.

        Parameters
        ----------
        release : str
            The release to remove, e.g. ``2024_06``.
        """
        self.drop_partition(release)
        self.forget_source(version=release)

class MappingCursor(BaseMappingCursor):
    """
    Translate identifiers through UniProt's mapping table.

    One query answers what used to need three cursors. The
    identifiers given are matched in the databases named by
    ``source``, their UniProtKB accessions are found, and the
    identifiers those accessions carry in the databases named by
    ``target`` are returned. Pinning either end to
    :attr:`~rotifer.db.methods.MappingCursor.UNIPROTKB` gives the
    accession itself, so:

    ``source=UNIPROTKB``
        the cross-references of a UniProtKB accession
    ``target=UNIPROTKB``
        the accession an external identifier belongs to
    neither
        a translation between two external databases

    Both ends are named per call rather than per cursor, since they
    describe the question rather than the table.

    Parameters
    ----------
    release : str, optional
        Restrict results to one UniProt release.
    batch_size : int, optional
        Number of identifiers sent to the server per query.
    **kwargs
        Connection parameters, passed to
        :class:`~rotifer.db.sql.clickhouse.core.BaseClickHouseCursor`.

    Examples
    --------
    Every cross-reference of an accession:

    >>> from rotifer.db.uniprot import clickhouse as ruch
    >>> mc = ruch.MappingCursor()                                  # doctest: +SKIP
    >>> mc.fetchall(["Q6GZX4"], source=mc.UNIPROTKB)               # doctest: +SKIP

    Which accession a RefSeq protein belongs to:

    >>> mc.fetchall(["YP_031579.1"], target=mc.UNIPROTKB)          # doctest: +SKIP

    From GenBank CDS to RefSeq and KEGG at once:

    >>> mc.fetchall(["AAT09660.1"], source=['EMBL-CDS'],           # doctest: +SKIP
    ...             target=['RefSeq','KEGG'])
    """

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
            Identifiers with no translation are registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()
        source = self.parse_databases(source)
        target = self.parse_databases(target)

        stack = []
        try:
            if len(targets) >= self.submit_threshold:
                # Too many to write into the statement: hand them over
                # as a table and let every branch refer to it.
                table = self.submit(targets)
                try:
                    binder = sqlmap.NamedBinder()
                    restrict = lambda column, _t=table: f'f.{column} IN (SELECT id FROM {_t})'
                    sql = sqlmap.mapping_query(self.qualified(), restrict, source, target,
                                               binder, release=self.release)
                    stack.append(self.query(sql, parameters=binder.parameters))
                finally:
                    self.cleanup()
            else:
                batches = list(self._batches(targets, self.batch_size))
                with sqlprog.Progress(total=len(targets), unit='ids', desc='querying',
                                      enabled=self.progress,
                                      position=1, leave=False) as bar:
                  for batch in batches:
                    binder = sqlmap.NamedBinder()

                    def restrict(column, _batch=batch, _binder=binder):
                        return f'f.{column} IN {_binder.collection(_batch)}'

                    sql = sqlmap.mapping_query(self.qualified(), restrict, source, target,
                                               binder, release=self.release)
                    stack.append(self.query(sql, parameters=binder.parameters))
                    bar.update(len(batch))
        except Exception as error:
            # An unreachable or broken server must not abort the caller:
            # registering the query as missing lets a delegator hand it
            # to the next backend, and retry stays True because another
            # attempt may well succeed.
            logger.error(f'Query to {self.qualified()} at {self.host} failed: {error}')
            self.update_missing(targets, error=f'ClickHouse query failed: {error}', retry=True)
            return self.empty()

        df = pd.concat(stack, ignore_index=True) if stack else self.empty()
        if not df.empty:
            df = df[self.columns].drop_duplicates().reset_index(drop=True)

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'No mapping found in {self.qualified()}', retry=False)

        return df

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
            The mappings of one batch. Input order is not preserved.
        """
        targets = self.parse_ids(accessions)
        for batch in self._batches(targets, self.batch_size):
            found = self.__getitem__(batch, source=source, target=target)
            if not found.empty:
                yield found

if __name__ == '__main__':
    pass
