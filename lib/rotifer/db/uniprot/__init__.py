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
"""

# Import external modules
import os
import types
import pandas as pd

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

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Fetch identifier mappings, dictionary style.

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
        data = super().__getitem__(accessions, *args, **kwargs)
        if isinstance(data, pd.DataFrame):
            return data
        if isinstance(data, types.NoneType):
            return self.empty()
        if not isinstance(data, list):
            data = [data]
        data = [ x for x in data if isinstance(x, pd.DataFrame) and not x.empty ]
        return pd.concat(data, ignore_index=True) if data else self.empty()

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
    host, port, database, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration, which is
        also where credentials belong.
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
            database = None,
            table = None,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        self._shared_attributes = ['progress','id_type','release','path','engine','host','port','database','table','batch_size','threads']
        self.id_type = id_type
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.database = database
        self.table = table
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

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
    host, port, database, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration.
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
            database = None,
            table = None,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        self._shared_attributes = ['progress','id_type','release','path','engine','host','port','database','table','batch_size','threads']
        self.id_type = id_type
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.database = database
        self.table = table
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

class MappingCursor(BaseUniProtDelegatorCursor):
    """
    Translate identifiers from one database into another.

    The ClickHouse backend answers this with a single server side
    join. The flat file backend needs two passes over
    ``idmapping.dat``, so falling back to it is expensive.

    Parameters
    ----------
    from_type : str
        Name of the database the queried identifiers belong to, as
        written in ``idmapping.dat``, e.g. ``EMBL-CDS``. Use
        ``UniProtKB-AC`` to start from UniProtKB accessions.
    to_type : str
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
    host, port, database, table : optional
        Where the ClickHouse backend should look. Each defaults to
        the matching entry of the
        :mod:`rotifer.db.uniprot.clickhouse` configuration.
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
    >>> mc = uniprot.MappingCursor(from_type='EMBL-CDS', to_type='RefSeq')  # doctest: +SKIP
    >>> mc.fetchall(["AAT09660.1"])  # doctest: +SKIP
    """

    _columns = ['from','accession','to']
    column = 'from'

    def __init__(
            self,
            from_type,
            to_type,
            readers = ['clickhouse','mirror'],
            writers = [],
            release = None,
            local_database_path = config['local_database_path'],
            engine = None,
            host = None,
            port = None,
            database = None,
            table = None,
            progress = True,
            tries = None,
            batch_size = None,
            threads = None,
            *args, **kwargs
        ):
        self._shared_attributes = ['progress','from_type','to_type','release','path','engine','host','port','database','table','batch_size','threads']
        self.from_type = from_type
        self.to_type = to_type
        self.release = release
        self.path = local_database_path
        self.engine = engine
        self.host = host
        self.port = port
        self.database = database
        self.table = table
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

if __name__ == '__main__':
    pass
