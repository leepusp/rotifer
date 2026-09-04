"""
Access data published by UniProt.

This package is the main entry point for UniProt identifier mappings.
The cursors defined here are delegators: each one combines several
backends and tries them in order until every requested identifier is
resolved, so that a query is answered by the fastest source that
knows the answer.

Two backends are available, and the default order puts the fast one
first:

:mod:`rotifer.db.uniprot.clickhouse`
    An indexed copy of ``idmapping.dat`` in ClickHouse. Point
    lookups are answered in milliseconds. Tried first.
:mod:`rotifer.db.uniprot.mirror`
    The flat files of a local UniProt mirror. Every query scans the
    whole file, which costs minutes, so this backend only receives
    the identifiers ClickHouse could not resolve, and covers the
    cases where the table is out of date, incomplete or unreachable.

A third source, :mod:`rotifer.db.uniprot.webapi`, queries UniProt's
REST service. It is not wired in as a backend yet, because it is
written as functions rather than as cursor classes.

Configuration
-------------
The module level ``config`` dictionary is loaded from
``~/.rotifer/etc/db/uniprot.yml`` when that file exists. Its keys
include the path of the local UniProt mirror and the mapping of
backend names to reader and writer modules.

Examples
--------
Fetch every cross-reference of two UniProtKB accessions, from
ClickHouse if it has them and from the flat file otherwise:

>>> from rotifer.db import uniprot
>>> ic = uniprot.IdMappingCursor()  # doctest: +SKIP
>>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP

Ask which backend answered, and what is still missing:

>>> ic.missing  # doctest: +SKIP

Starting from a mirror and an empty ClickHouse database, fill the
table once and query it from then on:

>>> ic = uniprot.IdMappingCursor(  # doctest: +SKIP
...     local_database_path="/scratch/global/databases/uniprot",
...     dbname="uniprot", release="2026_01", initialize='load')

The same thing on demand, instead of in one sitting: every query
answered by the mirror is stored, so the second time it is answered by
the table.

>>> ic = uniprot.IdMappingCursor(  # doctest: +SKIP
...     local_database_path="/scratch/global/databases/uniprot",
...     dbname="uniprot", release="2026_01", cache=True)
>>> ic.fetchall(["Q6GZX4"])   # scans the file, then stores what it found
>>> ic.fetchall(["Q6GZX4"])   # answered by ClickHouse
"""

# Import external modules
import os
import types
import pandas as pd
from copy import deepcopy

# Import rotifer modules
import rotifer
import rotifer.db.core
import rotifer.db.methods
import rotifer.db.delegator
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Configuration
config = loadConfig(__name__.replace('rotifer.',':'), defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"uniprot"),
    'readers': {
        'clickhouse': 'rotifer.db.uniprot.clickhouse',
        'mirror': 'rotifer.db.uniprot.mirror',
    },
    'writers': {
        'clickhouse': 'rotifer.db.uniprot.clickhouse',
    },
})

# Classes

class BaseUniProtDelegatorCursor(rotifer.db.methods.IdMappingCursor, rotifer.db.delegator.SequentialDelegatorCursor):
    """
    Shared behaviour of the UniProt delegator cursors.

    This class is not meant to be used directly. It exists so that
    every cursor in this package returns one dataframe, whichever
    backend answered: the generic delegator collects one result per
    backend and per batch, and this class concatenates them.

    See Also
    --------
    rotifer.db.delegator.SequentialDelegatorCursor : the delegation logic
    """

    #: Name of the backend that stores data, used by ``cache``.
    _store_backend = 'clickhouse'

    #: Both are query filters, so None means "no filter" and has to
    #: reach the backends. Without this a filter set on the delegator
    #: could be changed but never cleared.
    _nullable_attributes = frozenset({'id_type','release'})

    def __init__(self, *args, **kwargs):
        """
        Build the delegator, rejecting the old name for ``dbname``.

        A delegator keeps its own keywords rather than handing them
        all to its backends, so a caller still passing ``database``
        would have it quietly dropped here instead of reaching the
        ClickHouse cursor that would have complained.
        """
        if 'database' in kwargs:
            raise TypeError(
                "the 'database' parameter is now called 'dbname'; "
                f"pass dbname={kwargs['database']!r} instead"
            )
        super().__init__(*args, **kwargs)

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Fetch identifier mappings, dictionary style.

        Equivalent to :meth:`fetchall`, so that dictionary style
        access caches its results like the other two access styles.

        Parameters
        ----------
        accessions : str or iterable of str
            Database identifiers.

        Returns
        -------
        pandas.DataFrame
            The rows found by every backend consulted, concatenated.
            Identifiers no backend could resolve are registered in
            :attr:`~rotifer.db.core.BaseCursor.missing`.
        """
        return self.fetchall(accessions, *args, **kwargs)

    def fetchone(self, accessions, *args, **kwargs):
        """
        Iterate over identifier mappings, trying each backend in turn.

        Backends listed in ``readers`` are consulted in order and each
        one receives only the identifiers its predecessors could not
        resolve. Rows are handed to the backends listed in
        ``writers``, which is how ``cache`` stores what a query just
        retrieved.

        Parameters
        ----------
        accessions : str or iterable of str
            Database identifiers.

        Yields
        ------
        pandas.DataFrame
            One block of rows, as produced by the backend that found
            them.

        Note
        ----
        This overrides
        :meth:`rotifer.db.delegator.SequentialDelegatorCursor.fetchone`
        for one reason: the generic version hands every result to
        every writer, including results a writer just returned as a
        reader. Here the same backend can be both, so rows are never
        written back to the backend they came from, which would
        duplicate them.
        """
        targets = self.parse_ids(accessions)
        todo = deepcopy(targets)
        for position, name in enumerate(self.readers):
            if not todo:
                break
            if name not in self.cursors:
                continue
            cursor = self.cursors[name]
            for result in cursor.fetchone(todo, *args, **kwargs):
                found = self.getids(result, *args, **kwargs)
                done = todo.intersection(found)
                for earlier in self.readers[:position+1]:
                    if earlier in self.cursors:
                        self.cursors[earlier].remove_missing(done)
                self.remove_missing(done)
                self.update_missing(data=cursor._missing)
                for writer in self.writers:
                    if writer == name or writer not in self.cursors:
                        continue
                    rows = self._rows_to_store(result, name)
                    if not rows.empty:
                        self.cursors[writer].insert(rows)
                todo = todo - done
                yield result

    def _rows_to_store(self, result, source):
        """
        Choose which rows to hand to the writers.

        Parameters
        ----------
        result : pandas.DataFrame
            Rows a reader just returned.
        source : str
            Name of the backend that produced them.

        Returns
        -------
        pandas.DataFrame
            Rows laid out like the storage table.
        """
        return result

    @property
    def store(self):
        """
        The backend cursor that stores data.

        Returns
        -------
        object or None
            The ClickHouse cursor, or None when it is not among this
            delegator's backends.
        """
        return self.cursors.get(self._store_backend)

    def create(self, replace=False):
        """
        Create the ClickHouse database and table.

        Parameters
        ----------
        replace : bool, default False
            If True, drop an existing table before creating it. Every
            row it holds is lost.

        Returns
        -------
        bool
            Whether the table exists after the call.

        Raises
        ------
        ValueError
            If the ClickHouse backend is not among this delegator's
            readers or writers.

        Examples
        --------
        >>> from rotifer.db import uniprot
        >>> ic = uniprot.IdMappingCursor(dbname='uniprot')  # doctest: +SKIP
        >>> ic.create()  # doctest: +SKIP
        """
        store = self.store
        if isinstance(store, types.NoneType):
            raise ValueError(f'No {self._store_backend} backend: add it to readers or writers')
        return store.create(replace=replace)

    def load(self, source=None, release=None, method='auto', **kwargs):
        """
        Load a whole release from the mirror into ClickHouse.

        The table is created when it does not exist yet, then every
        row of the mirror's ``idmapping.dat`` is inserted. This is the
        one call that turns an empty database into one worth querying.

        Parameters
        ----------
        source : str or cursor, optional
            Where to read the mappings from. Defaults to this
            delegator's own mirror backend, so that the path given at
            construction is used.
        release : str, optional
            Value stored in the ``release`` column of every row
            loaded. Defaults to the delegator's ``release``.
        method : str, default 'auto'
            How to send the data. See
            :meth:`rotifer.db.uniprot.clickhouse.BaseIdMappingCursor.load`.
        **kwargs
            Passed on to the ClickHouse backend's ``load``.

        Returns
        -------
        int
            Number of rows in the table after the load.

        Raises
        ------
        ValueError
            If the ClickHouse backend, or a source to read from, is
            missing.

        Note
        ----
        A full release is a few billion rows and takes about an hour.
        The table is partitioned by release, so an interrupted load is
        cleaned up with
        ``ALTER TABLE ... DROP PARTITION '<release>'`` before trying
        again.

        Examples
        --------
        Point a cursor at a mirror and an empty database, then fill it:

        >>> from rotifer.db import uniprot
        >>> ic = uniprot.IdMappingCursor(  # doctest: +SKIP
        ...     local_database_path="/scratch/global/databases/uniprot",
        ...     dbname="uniprot", release="2026_01")
        >>> ic.load()  # doctest: +SKIP
        2647104040
        """
        store = self.store
        if isinstance(store, types.NoneType):
            raise ValueError(f'No {self._store_backend} backend: add it to readers or writers')
        if isinstance(source, types.NoneType):
            source = self.cursors.get('mirror')
            if isinstance(source, types.NoneType):
                source = self.path
        if isinstance(source, types.NoneType):
            raise ValueError('No mirror backend and no source given: nothing to load from')
        if self.progress:
            logger.warn(f'Loading the whole mapping table into {store.qualified_name}. This takes about an hour.')
        return store.load(source, release=release, method=method, **kwargs)

    def _initialize(self, initialize, strict=True):
        """
        Prepare the storage backend at construction time.

        Parameters
        ----------
        initialize : bool or str
            One of False, ``create`` or ``load``. See the
            ``initialize`` parameter of the cursors in this module.

        Raises
        ------
        ValueError
            If `initialize` is not one of the accepted values.
        """
        if not initialize:
            return
        if initialize is True:
            initialize = 'create'
        if initialize not in ('create','load'):
            raise ValueError(f"Unknown initialize {initialize}: expected False, 'create' or 'load'")
        store = self.store
        if isinstance(store, types.NoneType):
            raise ValueError(f'No {self._store_backend} backend: add it to readers or writers')
        try:
            if not store.has_table():
                self.create()
            if initialize == 'load' and store.is_empty():
                self.load()
        except Exception as error:
            # An explicit request must fail loudly; the create implied
            # by cache must not stop a session that can still read the
            # mirror
            if strict:
                raise
            logger.error(f'Could not prepare {store.qualified_name}, caching is off: {error}')
            self.writers = [ x for x in self.writers if x != self._store_backend ]

class IdMappingCursor(BaseUniProtDelegatorCursor):
    """
    Fetch the cross-references of UniProtKB accessions.

    Backends are tried in the order given by `readers`, and each one
    receives only the accessions the previous backends could not
    resolve.

    Parameters
    ----------
    readers : list of str, default ``['clickhouse', 'mirror']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules. Setting this to ``['clickhouse']``
        stores rows recovered from the flat file into the ClickHouse
        table, which is useful to fill gaps in an incomplete load.
    id_type : str or list of str, optional
        Restrict results to these cross-referenced databases, e.g.
        ``RefSeq``. Shared with every backend.
    release : str, optional
        Restrict the ClickHouse backend to one UniProt release.
    local_database_path : str, optional
        Root directory of the local UniProt mirror, used by the
        ``mirror`` backend. Defaults to the ``local_database_path``
        configuration entry.
    engine : str, optional
        Matching engine of the ``mirror`` backend, one of ``auto``,
        ``arrow`` or ``python``.
    host, port, dbname, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration, which is
        also where credentials belong.
    initialize : bool or str, default False
        What to do about the ClickHouse table when the cursor is
        built:

        ``False``
            Nothing. The table is expected to exist.
        ``'create'`` or True
            Create the database and table when they are missing, so
            that a cursor can be pointed at an empty database.
        ``'load'``
            Create them, and when no row of the release is present,
            read the whole mirror into the table. This is the one
            call that turns an empty database into one worth
            querying, and it takes about an hour.

    cache : bool, default False
        Store rows retrieved from the mirror into ClickHouse as
        ``fetchall``, ``fetchone`` and item access return them, so
        that repeating a query is answered by the table. Implies
        ``initialize='create'``, and never writes rows back to the
        backend that produced them.
    progress : bool, default True
        Whether to print progress messages.
    batch_size : int, optional
        Number of identifiers per query, used by the ClickHouse
        backend.
    threads : int, optional
        Number of worker processes, used by the ``mirror`` backend.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.IdMappingCursor : the fast backend
    rotifer.db.uniprot.mirror.IdMappingCursor : the fallback backend
    CrossReferenceCursor : the reverse lookup

    Examples
    --------
    >>> from rotifer.db import uniprot
    >>> ic = uniprot.IdMappingCursor()  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4","Q6GZX3"])  # doctest: +SKIP

    Only their RefSeq proteins, straight from the flat file:

    >>> ic = uniprot.IdMappingCursor(readers=['mirror'], id_type='RefSeq')  # doctest: +SKIP
    >>> df = ic.fetchall(["Q6GZX4"])  # doctest: +SKIP
    """

    column = 'accession'

    def __init__(
            self,
            readers = ['clickhouse','mirror'],
            writers = [],
            id_type = None,
            release = None,
            local_database_path = config['local_database_path'],
            engine = None,
            host = None,
            port = None,
            dbname = None,
            table = None,
            initialize = False,
            cache = False,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        self._shared_attributes = ['progress','id_type','release','path','engine','host','port','dbname','table','batch_size','threads']
        self.id_type = id_type
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.dbname = dbname
        self.table = table
        writers = list(writers)
        if cache and self._store_backend not in writers:
            writers.append(self._store_backend)
        self.cache = cache
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        # Caching needs somewhere to write, so it implies a table
        self._initialize(initialize or (cache and 'create'), strict=bool(initialize))

class CrossReferenceCursor(BaseUniProtDelegatorCursor):
    """
    Fetch the UniProtKB accessions of identifiers from other databases.

    This is the reverse of :class:`IdMappingCursor`. The ClickHouse
    backend answers it from the table's ``by_id`` projection; the
    flat file backend scans the third column of ``idmapping.dat``.

    Parameters
    ----------
    readers : list of str, default ``['clickhouse', 'mirror']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules.
    id_type : str or list of str, optional
        Restrict the search to these cross-referenced databases.
        Naming the database makes the ClickHouse query cheaper.
    release : str, optional
        Restrict the ClickHouse backend to one UniProt release.
    local_database_path : str, optional
        Root directory of the local UniProt mirror.
    engine : str, optional
        Matching engine of the ``mirror`` backend.
    host, port, dbname, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration.
    initialize : bool or str, default False
        What to do about the ClickHouse table when the cursor is
        built:

        ``False``
            Nothing. The table is expected to exist.
        ``'create'`` or True
            Create the database and table when they are missing, so
            that a cursor can be pointed at an empty database.
        ``'load'``
            Create them, and when no row of the release is present,
            read the whole mirror into the table. This is the one
            call that turns an empty database into one worth
            querying, and it takes about an hour.

    cache : bool, default False
        Store rows retrieved from the mirror into ClickHouse as
        ``fetchall``, ``fetchone`` and item access return them, so
        that repeating a query is answered by the table. Implies
        ``initialize='create'``, and never writes rows back to the
        backend that produced them.
    progress : bool, default True
        Whether to print progress messages.
    batch_size : int, optional
        Number of identifiers per query.
    threads : int, optional
        Number of worker processes used by the ``mirror`` backend.

    See Also
    --------
    IdMappingCursor : the forward lookup
    rotifer.db.uniprot.clickhouse.CrossReferenceCursor : the fast backend

    Examples
    --------
    >>> from rotifer.db import uniprot
    >>> xc = uniprot.CrossReferenceCursor(id_type='RefSeq')  # doctest: +SKIP
    >>> xc.fetchall(["YP_031579.1"])  # doctest: +SKIP
    """

    column = 'id'

    def _rows_to_store(self, result, source):
        """
        Expand reverse lookup rows into whole accession groups.

        A reverse lookup returns only the rows whose identifier was
        queried, not every row of the accessions behind them. Storing
        those as they are would break the invariant the forward lookup
        depends on, that an accession present in the table is present
        in full: a later
        :class:`IdMappingCursor` query for one of these accessions
        would be answered from the table alone and silently return a
        fraction of its cross-references. So the mirror is asked for
        the complete rows of every accession found, which costs one
        further scan of the file.

        Parameters
        ----------
        result : pandas.DataFrame
            Rows a reader just returned.
        source : str
            Name of the backend that produced them.

        Returns
        -------
        pandas.DataFrame
            Every row of the accessions named in `result`.
        """
        if result.empty:
            return result
        from rotifer.db.uniprot import mirror as rum
        settings = { k: v for k, v in (
            ('path', self.path), ('threads', self.threads),
            ('engine', self.engine), ('progress', self.progress),
        ) if not isinstance(v, types.NoneType) }
        if self.progress:
            logger.warn('Caching a reverse lookup: scanning the mirror again to store whole accession groups')
        return rum.IdMappingCursor(**settings).fetchall(set(result.accession))

    def __init__(
            self,
            readers = ['clickhouse','mirror'],
            writers = [],
            id_type = None,
            release = None,
            local_database_path = config['local_database_path'],
            engine = None,
            host = None,
            port = None,
            dbname = None,
            table = None,
            initialize = False,
            cache = False,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        self._shared_attributes = ['progress','id_type','release','path','engine','host','port','dbname','table','batch_size','threads']
        self.id_type = id_type
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.dbname = dbname
        self.table = table
        writers = list(writers)
        if cache and self._store_backend not in writers:
            writers.append(self._store_backend)
        self.cache = cache
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        # Caching needs somewhere to write, so it implies a table
        self._initialize(initialize or (cache and 'create'), strict=bool(initialize))

class MappingCursor(BaseUniProtDelegatorCursor):
    """
    Translate identifiers from one database into another.

    The ClickHouse backend answers this with a single server side
    join. The flat file backend needs two passes over
    ``idmapping.dat``, so falling back to it is expensive.

    Parameters
    ----------
    source : str
        Name of the database the queried identifiers belong to, as
        written in ``idmapping.dat``, e.g. ``EMBL-CDS``. Use
        ``UniProtKB-AC`` to start from UniProtKB accessions.
    target : str
        Name of the database to translate into, e.g. ``RefSeq``. Use
        ``UniProtKB-AC`` to translate into UniProtKB accessions.
    readers : list of str, default ``['clickhouse', 'mirror']``
        Backend reader modules, tried in order.
    writers : list of str, default []
        Backend writer modules. Leave this empty: the rows this
        cursor returns do not match the layout of the ClickHouse
        table, so they cannot be stored there.
    release : str, optional
        Restrict the ClickHouse backend to one UniProt release.
    local_database_path : str, optional
        Root directory of the local UniProt mirror.
    engine : str, optional
        Matching engine of the ``mirror`` backend.
    host, port, dbname, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration.
    initialize : bool or str, default False
        What to do about the ClickHouse table when the cursor is
        built:

        ``False``
            Nothing. The table is expected to exist.
        ``'create'`` or True
            Create the database and table when they are missing, so
            that a cursor can be pointed at an empty database.
        ``'load'``
            Create them, and when no row of the release is present,
            read the whole mirror into the table. This is the one
            call that turns an empty database into one worth
            querying, and it takes about an hour.

    cache : bool, default False
        Store rows retrieved from the mirror into ClickHouse as
        ``fetchall``, ``fetchone`` and item access return them, so
        that repeating a query is answered by the table. Implies
        ``initialize='create'``, and never writes rows back to the
        backend that produced them.
    progress : bool, default True
        Whether to print progress messages.
    batch_size : int, optional
        Number of identifiers per query.
    threads : int, optional
        Number of worker processes used by the ``mirror`` backend.

    Attributes
    ----------
    columns : list of str
        ``['from', 'accession', 'to']``.

    See Also
    --------
    rotifer.db.uniprot.clickhouse.MappingCursor : the fast backend
    rotifer.db.uniprot.webapi.idmapping : the same query, run by UniProt

    Examples
    --------
    >>> from rotifer.db import uniprot
    >>> mc = uniprot.MappingCursor(source='EMBL-CDS', target='RefSeq')  # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"])  # doctest: +SKIP
    """

    _columns = ['from','accession','to']
    column = 'from'

    def __init__(
            self,
            source,
            target,
            readers = ['clickhouse','mirror'],
            writers = [],
            release = None,
            local_database_path = config['local_database_path'],
            engine = None,
            host = None,
            port = None,
            dbname = None,
            table = None,
            initialize = False,
            cache = False,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        if cache:
            raise ValueError(
                'MappingCursor cannot cache: its rows are from/accession/to, which do not '
                'match the layout of the storage table. Cache with IdMappingCursor, or load '
                'the whole release with load().'
            )
        self._shared_attributes = ['progress','source','target','release','path','engine','host','port','dbname','table','batch_size','threads']
        self.source = source
        self.target = target
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.dbname = dbname
        self.table = table
        writers = list(writers)
        if cache and self._store_backend not in writers:
            writers.append(self._store_backend)
        self.cache = cache
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)
        # Caching needs somewhere to write, so it implies a table
        self._initialize(initialize or (cache and 'create'), strict=bool(initialize))

if __name__ == '__main__':
    pass
