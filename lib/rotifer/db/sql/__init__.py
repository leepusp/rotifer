"""
Access data held in SQL databases.

This package gathers the backends that answer queries with SQL, so
that what they share is stated once and each module only adds what its
engine does differently:

:mod:`rotifer.db.sql.clickhouse`
    A ClickHouse server, for tables too large to scan, such as
    UniProt's identifier mappings.
:mod:`rotifer.db.sql.sqlite3`
    A local SQLite3 file, for data that travels with the analysis
    rather than living on a server.

Provenance
----------
Both record where their data came from, in a table named
``rotifer_sources``, with the same columns and through the same
methods: :meth:`record_source` writes a row when a load completes, and
:meth:`content_id` reports the identity of the file that was loaded.
There is one such registry per database rather than one per data
table, because provenance is the same kind of fact whatever was
loaded.

That record is what allows a delegator to notice that a cursor reading
a file directly holds nothing a table loaded from that same file does
not, and to skip the slower of the two. A table whose load was never
recorded reports nothing, and is never taken for a complete copy.

See Also
--------
rotifer.db.core.BaseCursor.content_id : what the recorded identity means
rotifer.db.delegator : where it is used to skip a redundant backend
"""

# base SQL module
