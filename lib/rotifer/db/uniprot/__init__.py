"""
Access data published by UniProt.

This package is the main entry point for UniProt identifier mappings.
The cursors defined here are delegators: each one combines several
backends and tries them in order until every requested identifier is
resolved, so that a query is answered by the fastest source that
knows the answer.

Three backends are available, and the default order tries them
cheapest first:

:mod:`rotifer.db.uniprot.clickhouse`
    An indexed copy of ``idmapping.dat``. Point lookups are answered
    in milliseconds, so it is asked first.
:mod:`rotifer.db.uniprot.webapi`
    UniProt's REST service. Always current and needs no local copy,
    and answers in seconds. It comes second because it is the only
    source for what the local table does not carry, and because a
    round trip is still far cheaper than the alternative.
:mod:`rotifer.db.uniprot.mirror`
    The flat files of a local UniProt mirror. Every query scans about
    90 GB, which costs a minute and a half whatever is asked, so it is
    asked last: it is the source that has everything, kept for what
    the other two could not answer and for when they are unreachable
    or out of date.

The three do not hold the same vocabulary. The table and the mapping
service each carry databases the other does not, and some, Pfam and GO
among them, are cross-references of an entry rather than identifier
mappings and are in neither. A backend is therefore asked only for the
databases it says it can map, and a database none of them supports is
reported in ``missing`` under its own name rather than silently
returning nothing.

The web backend also covers what the other two do not at all:
sequences, proteomes, taxonomy and UniProt's own search, through the
cursors documented in that module.

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
>>> ic = uniprot.MappingCursor()  # doctest: +SKIP
>>> df = ic.fetchall(["Q6GZX4","Q6GZX3"], source=ic.UNIPROTKB)  # doctest: +SKIP

Ask which backend answered, and what is still missing:

>>> ic.missing  # doctest: +SKIP

Starting from a mirror and an empty ClickHouse database, fill the
table once and query it from then on:

>>> ic = uniprot.MappingCursor(  # doctest: +SKIP
...     local_database_path="/scratch/global/databases/uniprot",
...     dbname="rotifer", release="2026_01", initialize='load')

The same thing on demand, instead of in one sitting: every query
answered by the mirror is stored, so the second time it is answered by
the table.

>>> ic = uniprot.MappingCursor(  # doctest: +SKIP
...     local_database_path="/scratch/global/databases/uniprot",
...     dbname="rotifer", release="2026_01", cache=True)
>>> ic.fetchall(["Q6GZX4"], source=ic.UNIPROTKB)   # scans the file, then stores what it found
>>> ic.fetchall(["Q6GZX4"], source=ic.UNIPROTKB)   # answered by ClickHouse
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
        # Registered but not enabled by default: every query is a
        # round trip to UniProt, so it is opted into per cursor.
        'webapi': 'rotifer.db.uniprot.webapi',
    },
    'writers': {
        'clickhouse': 'rotifer.db.uniprot.clickhouse',
    },
})

# Classes

class BaseUniProtDelegatorCursor(rotifer.db.methods.MappingCursor, rotifer.db.delegator.SequentialDelegatorCursor):
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
        consulted = dict()

        # A query names identifiers and, at both ends, databases. Any
        # of the three can go unanswered, and they fail differently: an
        # identifier absent from one copy of the data may be in the
        # next, while a database a backend does not carry is a gap in
        # that backend whatever it is asked. Both ends are therefore
        # tracked alongside the identifiers, so that what is handed on
        # is narrowed to what the next backend can actually answer.
        source = self.parse_databases(kwargs.pop('source', None))
        target = self.parse_databases(kwargs.pop('target', None))
        wanted_source = set(source or [])
        wanted_target = set(target or [])
        served_source, served_target = set(), set()

        for position, name in enumerate(self.readers):
            # What is still owed at each end. An end that is already
            # covered falls back to the whole request rather than to
            # nothing: it is a constraint on the query, not a thing to
            # be collected, so a query still owed at the other end
            # needs it stated in full.
            pending_source = sorted(wanted_source - served_source) if wanted_source else None
            pending_target = sorted(wanted_target - served_target) if wanted_target else None
            ask_source = pending_source or source
            ask_target = pending_target or target

            # Finding every identifier is not the same as answering
            # every question: a database no backend has looked at yet
            # is still owed, even when nothing is left to look up.
            if not todo and not pending_source and not pending_target:
                break
            if name not in self.cursors:
                continue
            cursor = self.cursors[name]
            served_by = self.redundant(cursor, consulted)
            if served_by:
                logger.info(f'Skipping backend {name}: same data as {served_by}')
                continue

            # Note what this backend holds before deciding whether to
            # ask it. A backend passed over for lacking a database has
            # still told us what its data contains, and a later backend
            # holding that same data lacks that database too: without
            # this, being skipped here would hide the very fact that
            # spares the next one a pointless scan.
            content = self.content_of(cursor)
            if not isinstance(content, types.NoneType):
                consulted.setdefault(content, name)

            # Ask each backend only for the databases it says it can
            # answer for, so that an unsupported one falls through to
            # the next backend instead of coming back empty and looking
            # like an absent identifier.
            here_source = cursor.supported(ask_source)
            here_target = cursor.supported(ask_target)
            if (ask_source and not here_source) or (ask_target and not here_target):
                logger.info(f'Skipping backend {name}: none of the databases still owed are available there')
                continue

            # Identifiers already found still have to be asked about
            # when a database remains uncovered: the rows owed are the
            # ones that database would have contributed.
            asking = todo if todo else deepcopy(targets)

            for result in cursor.fetchone(asking, source=here_source, target=here_target, *args, **kwargs):
                found = self.getids(result, *args, **kwargs)
                done = todo.intersection(found)
                for earlier in self.readers[:position+1]:
                    if earlier in self.cursors:
                        self.cursors[earlier].remove_missing(done)
                self.remove_missing(done)
                for writer in self.writers:
                    if writer == name or writer not in self.cursors:
                        continue
                    rows = self._rows_to_store(result, name)
                    if not rows.empty:
                        self.cursors[writer].insert(rows)
                todo = todo - done
                if isinstance(result, pd.DataFrame) and not result.empty:
                    served_source.update(result.source_type.dropna().astype(str))
                    served_target.update(result.target_type.dropna().astype(str))
                yield result

            # A backend that finds nothing yields nothing, so what it
            # could not do has to be collected once it is exhausted
            self.absorb_missing(cursor)

            # Entries some backend declared final will not be found by
            # any of the others either
            todo = todo - self.missing_ids(final=True)

            # Databases this backend could have answered for count as
            # covered even when they returned nothing: the answer is
            # then simply that there is no such mapping, which the next
            # backend would only repeat.
            served_source.update(x for x in (here_source or []) if x not in cursor.unsupported(here_source))
            served_target.update(x for x in (here_target or []) if x not in cursor.unsupported(here_target))
            if isinstance(source, types.NoneType):
                wanted_source.update(served_source)
            if isinstance(target, types.NoneType):
                wanted_target.update(served_target)

        self._register_database_failures('source', wanted_source - served_source)
        self._register_database_failures('target', wanted_target - served_target)

    def _register_database_failures(self, end, databases):
        """
        Record databases no backend was able to answer for.

        These are not missing identifiers but missing capabilities, so
        they are registered under the database name and said to be
        final: no backend carries them, and asking again would not
        change that.

        Parameters
        ----------
        end : str
            Which end of the mapping they were asked for, ``source``
            or ``target``.
        databases : iterable of str
            Names no backend could serve.
        """
        for database in sorted(databases):
            if database == self.UNIPROTKB:
                continue
            self.update_missing(
                [database],
                error = f'No backend can map {end} database {database!r}',
                retry = False,
                final = True,
            )

    def _rows_to_store(self, result, backend):
        """
        Choose which rows to hand to the writers.

        A mapping names both of its ends, while the storage table
        holds one row per cross-reference, so each end that is a real
        database contributes a row and the ends that are the
        accession itself contribute none: the accession is the table's
        key rather than a row of it.

        Parameters
        ----------
        result : pandas.DataFrame
            Rows a reader just returned, with this cursor's columns.
        backend : str
            Name of the backend that produced them.

        Returns
        -------
        pandas.DataFrame
            Rows laid out like the storage table, i.e. ``accession``,
            ``id_type`` and ``id``.
        """
        table_columns = ['accession','id_type','id']
        if not isinstance(result, pd.DataFrame) or result.empty:
            return pd.DataFrame([], columns=table_columns)
        sides = []
        for value, kind in (('source','source_type'), ('target','target_type')):
            if value not in result.columns or kind not in result.columns:
                continue
            side = result[result[kind] != self.UNIPROTKB]
            if side.empty:
                continue
            sides.append(side[['accession', kind, value]]
                         .rename(columns={kind:'id_type', value:'id'}))
        if not sides:
            return pd.DataFrame([], columns=table_columns)
        return pd.concat(sides, ignore_index=True)[table_columns].drop_duplicates()

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
        >>> ic = uniprot.MappingCursor(dbname='rotifer')  # doctest: +SKIP
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
            :meth:`rotifer.db.uniprot.clickhouse.BaseMappingCursor.load`.
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
        >>> ic = uniprot.MappingCursor(  # doctest: +SKIP
        ...     local_database_path="/scratch/global/databases/uniprot",
        ...     dbname="rotifer", release="2026_01")
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
            logger.warning(f'Loading the whole mapping table into {store.qualified_name}. This takes about an hour.')
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

class MappingCursor(BaseUniProtDelegatorCursor):
    """
    Translate identifiers through UniProt, from the fastest source.

    One cursor answers every direction of the question. The
    identifiers given are matched in the databases named by
    ``source``, their UniProtKB accessions are found, and the
    identifiers those accessions carry in the databases named by
    ``target`` are returned. Pinning either end to
    :attr:`~rotifer.db.methods.MappingCursor.UNIPROTKB` gives the
    accession itself:

    ``source=UNIPROTKB``
        the cross-references of a UniProtKB accession
    ``target=UNIPROTKB``
        the accession an external identifier belongs to
    neither
        a translation between two external databases

    Both ends are arguments of :meth:`fetchall` and :meth:`fetchone`
    rather than of the constructor, because they describe a question
    and not a data source: one cursor can be asked many of them.

    Parameters
    ----------
    readers : list of str, default ``['clickhouse', 'webapi', 'mirror']``
        Backend reader modules, tried in order, cheapest first. The
        web service comes before the mirror because it answers in
        seconds where a mirror scan costs about ninety, whatever is
        asked of it; the mirror stays last as the source that has
        everything, for whatever the other two could not map.
    writers : list of str, default []
        Backend writer modules.
    release : str, optional
        Restrict the ClickHouse backend to one UniProt release.
    local_database_path : str, optional
        Root directory of the local UniProt mirror.
    engine : str, optional
        Matching engine of the ``mirror`` backend.
    host, port, dbname, table : optional
        Where the ClickHouse backend should look. Each defaults to
        that backend's own configuration.
    initialize : bool or str, default False
        Create the storage table, and load it, before querying.
    cache : bool, default False
        Store what the slower backends return, so that the next
        query for the same identifiers is answered by the table.
    progress : bool, default True
        Whether backends report progress.

    Examples
    --------
    Every cross-reference of an accession:

    >>> from rotifer.db import uniprot
    >>> mc = uniprot.MappingCursor()                              # doctest: +SKIP
    >>> mc.fetchall(["Q6GZX4"], source=mc.UNIPROTKB)              # doctest: +SKIP

    Which accession a RefSeq protein belongs to:

    >>> mc.fetchall(["YP_031579.1"], target=mc.UNIPROTKB)         # doctest: +SKIP

    From GenBank CDS to RefSeq and KEGG at once:

    >>> mc.fetchall(["AAT09660.1"], source=['EMBL-CDS'],          # doctest: +SKIP
    ...             target=['RefSeq','KEGG'])
    """

    def __init__(
            self,
            readers = ['clickhouse','webapi','mirror'],
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
        self._shared_attributes = ['progress','release','path','engine','host','port','dbname','table','batch_size','threads']
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
