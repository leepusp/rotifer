"""
Access data stored in ClickHouse.

This package holds what every ClickHouse backed cursor needs,
regardless of the data it carries: opening a connection, running
statements, describing and creating tables, moving large files in, and
submitting long lists of identifiers to query against.

:mod:`rotifer.db.clickhouse.core`
    :class:`~rotifer.db.clickhouse.core.BaseClickHouseCursor`, the
    parent class of the cursors in this and other packages.

The first cursors built on it live in
:mod:`rotifer.db.uniprot.clickhouse`, which adds only what is specific
to UniProt's identifier mappings.

Configuration
-------------
Connection defaults are read from ``~/.rotifer/etc/db/clickhouse.yml``
when that file exists. A package that keeps its data in its own
database is expected to layer its own configuration on top, as
:mod:`rotifer.db.uniprot.clickhouse` does, so that a server can be
named once for every cursor or separately for each.

Note
----
``port`` below is ClickHouse's HTTP port, which is not necessarily the
port a given server listens on; the native protocol usually sits on
9000 and is not what these cursors speak.
"""

import rotifer
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Defaults shared by every ClickHouse cursor
_defaults = {
    'host': 'localhost',
    'port': 8123,
    'user': 'default',
    'password': '',
    'dbname': 'default',
    'table': '',
    'secure': False,
    'batch_size': 5000,
    'submit_threshold': 1000,
    'executable': 'clickhouse',
}
config = loadConfig(__name__.replace('rotifer.',':'), defaults = _defaults)
