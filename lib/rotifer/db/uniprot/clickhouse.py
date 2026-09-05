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

Everything here that is not about identifier mappings lives in
:mod:`rotifer.db.clickhouse.core`, which any other ClickHouse backed
cursor can build on.

The table itself is created and populated through
:meth:`BaseIdMappingCursor.create` and
:meth:`BaseIdMappingCursor.load`, from a local copy of the flat file
read by :class:`rotifer.db.uniprot.mirror.IdMappingCursor`. Its schema
lives in ``share/rotifer/db/uniprot/clickhouse/idmapping.sql``.

Configuration
-------------
Connection parameters are read from ``~/.rotifer/etc/db/uniprot/clickhouse.yml``
and fall back to the shared :mod:`rotifer.db.clickhouse` defaults, so
a server can be named once for every cursor or separately here.
"""

# Dependencies
import types
import typing
import pandas as pd

# Rotifer
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.clickhouse.core
from rotifer.db.clickhouse import config as clickhouse_config
from rotifer.core import functions as rcf
logger = rotifer.logging.getLogger(__name__)

#: Kept so that ``from rotifer.db.uniprot.clickhouse import
#: BaseClickHouseCursor`` still works; the class itself now lives in
#: :mod:`rotifer.db.clickhouse.core`.
BaseClickHouseCursor = rotifer.db.clickhouse.core.BaseClickHouseCursor

# Defaults: the shared connection settings, with what UniProt adds
_defaults = dict(clickhouse_config)
_defaults.update({
    'dbname': 'uniprot',
    'table': 'idmapping',
    'release': '',
    'chunksize': 5000000,
})
config = rcf.loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)

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
        for key in ('host','port','user','password','dbname','table',
                    'secure','batch_size','submit_threshold'):
            kwargs.setdefault(key, config[key])
        super().__init__(*args, **kwargs)
        self.id_type = id_type
        self.release = release
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

        See Also
        --------
        rotifer.db.clickhouse.core.BaseClickHouseCursor.create : the generic form

        Examples
        --------
        >>> from rotifer.db.uniprot import clickhouse as ruch
        >>> ic = ruch.IdMappingCursor(release='2026_01')  # doctest: +SKIP
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
        mirror : str or rotifer.db.uniprot.mirror.IdMappingCursor
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
        if isinstance(mirror, rum.IdMappingCursor):
            reader = mirror
        else:
            reader = rum.IdMappingCursor(path=mirror, progress=self.progress)
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

        if method == 'client':
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
        super().insert(data)

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
        self.drop_partition(release)

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
    the identifiers of `source` are looked up in the table, and
    every identifier of `target` sharing their UniProtKB accession is
    returned. The two lookups are done in one server side join, so the
    intermediate accessions never travel over the network.

    Parameters
    ----------
    source : str
        Name of the database the queried identifiers belong to, as
        written in ``idmapping.dat``, e.g. ``EMBL-CDS``. Use
        ``UniProtKB-AC`` to start from UniProtKB accessions
        themselves.
    target : str
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
    >>> mc = ruch.MappingCursor(source='EMBL-CDS', target='RefSeq')  # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"])  # doctest: +SKIP
    """

    _columns = ['from','accession','to']
    column = 'from'

    def __init__(self, source, target, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.source = source
        self.target = target

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
            Columns ``from``, ``accession`` and ``to``: the queried
            identifier, the UniProtKB accession that links it to the
            result, and the identifier in the database named by
            ``target``. Identifiers with no translation are
            registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        targets = self.parse_ids(accessions)
        if not targets:
            return self.empty()

        # UniProtKB accessions are the join key, not rows of the table
        source_column = "accession" if self.source == "UniProtKB-AC" else "id"
        target_column = "accession" if self.target == "UniProtKB-AC" else "id"

        stack = []
        try:
            stack = self._mapping_batches(targets, source_column, target_column)
        except Exception as error:
            logger.error(f'Query to {self.qualified_name} at {self.host} failed: {error}')
            self.update_missing(targets, error=f'ClickHouse query failed: {error}', retry=True)
            return self.empty()

        df = pd.concat(stack, ignore_index=True) if stack else self.empty()

        missing = targets.difference(self.getids(df))
        if missing:
            self.update_missing(missing, error=f'No {self.target} identifier found for this {self.source} identifier', retry=False)

        return df

    def _mapping_batches(self, targets, source_column, target_column):
        """
        Run the join, one batch of identifiers at a time.

        Parameters
        ----------
        targets : set of str
            Identifiers to translate.
        source_column, target_column : str
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
                return [ self._mapping_query(f'f.{source_column} IN (SELECT id FROM {table})', {}, source_column, target_column) ]
            finally:
                self.cleanup()

        stack = []
        for batch in self._batches(targets, self.batch_size):
            parameters = {'targets': tuple(batch)}
            stack.append(self._mapping_query(f'f.{source_column} IN %(targets)s', parameters, source_column, target_column))
        return stack

    def _mapping_query(self, restriction, parameters, source_column, target_column):
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
        source_column, target_column : str
            Names of the columns holding the queried and the returned
            identifiers, either ``accession`` or ``id``.

        Returns
        -------
        pandas.DataFrame
            Columns ``from``, ``accession`` and ``to``.
        """
        left = [restriction]
        if source_column == "id":
            left.append("f.id_type = %(source)s")
            parameters['source'] = self.source
        right = []
        if target_column == "id":
            right.append("t.id_type = %(target)s")
            parameters['target'] = self.target
        if self.release:
            left.append("f.release = %(release)s")
            right.append("t.release = %(release)s")
            parameters['release'] = self.release
        where = " AND ".join(left + right)
        return self.query(
            f'SELECT DISTINCT f.{source_column} AS `from`, f.accession AS accession, t.{target_column} AS `to`'
            f' FROM {self.qualified_name} AS f'
            f' INNER JOIN {self.qualified_name} AS t ON f.accession = t.accession'
            f' WHERE {where}'
            f' ORDER BY `from`, accession, `to`',
            parameters = parameters,
        )

if __name__ == '__main__':
    pass
