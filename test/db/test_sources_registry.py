#!/usr/bin/env python3
"""Tests for the SQL backends' load registry, `rotifer_sources`.

A SQL backend records where its data came from so a delegator can tell that
some other cursor reading that same file has nothing to add. The registry is
one table per database rather than one per data table, because provenance is
the same kind of fact whatever was loaded.

What has to hold, and why:

  * a table with no record yields content_id() None. That is the completeness
    gate: a load that was interrupted, or one that predates this bookkeeping,
    leaves no row, and None is never equal to anything, so nothing is skipped
    on its account
  * the identity is (name, size, mtime) of the source file, not a hash. The
    file this was built for is 90 GB, and hashing it costs about a minute,
    which is the cost the comparison exists to avoid
  * the value must match, byte for byte, what a cursor reading the file
    directly reports, or the two can never be recognised as the same data
  * a database holding several versions only answers when told which one to
    look at, since otherwise it holds more than the file does

Only SQLite3 is exercised here: it needs nothing but a temporary file, whereas
the ClickHouse half of this needs a server. Both implement the same names and
the same shape, and the last test pins that.

Run under pytest, or standalone: python test/db/test_sources_registry.py
"""

import os
import tempfile

import pytest

from rotifer.db.sql import sqlite3 as rdss


def cursor(tmpdir):
    return rdss.BaseSQLite3Cursor(path=os.path.join(tmpdir, 'test.sqlite3'))


def datafile(tmpdir, name='idmapping.dat', content=b'x' * 1024):
    path = os.path.join(tmpdir, name)
    with open(path, 'wb') as fh:
        fh.write(content)
    return path


def test_an_unrecorded_database_says_nothing(tmp_path):
    """The completeness gate: no record means no claim, so no skipping."""
    assert cursor(str(tmp_path)).content_id() is None


def test_recording_a_load_makes_it_identifiable(tmp_path):
    cur = cursor(str(tmp_path))
    path = datafile(str(tmp_path))
    assert cur.record_source(path, rows=7) is True
    info = os.stat(path)
    assert cur.content_id() == f'idmapping.dat:{info.st_size}:{int(info.st_mtime)}'


def test_the_identity_is_the_files_not_the_copys(tmp_path):
    """Two backends holding one file must produce the same string, or neither
    can be recognised as a substitute for the other. This is the exact format
    rotifer.db.uniprot.mirror derives from stat()."""
    cur = cursor(str(tmp_path))
    path = datafile(str(tmp_path))
    cur.record_source(path, rows=7)
    info = os.stat(path)
    mirror_style = f'{os.path.basename(path)}:{info.st_size}:{int(info.st_mtime)}'
    assert cur.content_id() == mirror_style


def test_a_missing_source_file_is_not_recorded(tmp_path):
    """Recording is a claim about a file; it cannot be made about one that is
    not there."""
    cur = cursor(str(tmp_path))
    assert cur.record_source(os.path.join(str(tmp_path), 'gone.dat'), rows=1) is False
    assert cur.content_id() is None


def test_records_are_kept_per_table(tmp_path):
    """One registry serves the whole database, so it has to say which table
    each row is about."""
    cur = cursor(str(tmp_path))
    cur.record_source(datafile(str(tmp_path), 'a.dat'), rows=1, table='features')
    cur.record_source(datafile(str(tmp_path), 'b.dat'), rows=1, table='other')
    assert cur.content_id().startswith('a.dat:')      # cursor.table is 'features'
    cur.table = 'other'
    assert cur.content_id().startswith('b.dat:')


def test_several_versions_are_ambiguous_without_a_filter(tmp_path):
    """Holding two releases at once is not the same as holding either one, so
    an unversioned cursor must not claim to be a substitute for a single file."""
    cur = cursor(str(tmp_path))
    cur.record_source(datafile(str(tmp_path), 'r1.dat'), rows=1, version='2026_01')
    cur.record_source(datafile(str(tmp_path), 'r2.dat'), rows=1, version='2026_02')
    assert cur.content_id() is None


def test_a_version_filter_resolves_the_ambiguity(tmp_path):
    class Versioned(rdss.BaseSQLite3Cursor):
        @property
        def source_version(self):
            return '2026_02'

    cur = Versioned(path=os.path.join(str(tmp_path), 'test.sqlite3'))
    cur.record_source(datafile(str(tmp_path), 'r1.dat'), rows=1, version='2026_01')
    cur.record_source(datafile(str(tmp_path), 'r2.dat'), rows=1, version='2026_02')
    assert cur.content_id().startswith('r2.dat:')


def test_reloading_replaces_rather_than_accumulates(tmp_path):
    """Loading the same table again describes the same slot, so the record is
    replaced; otherwise the table would look ambiguous to itself."""
    cur = cursor(str(tmp_path))
    path = datafile(str(tmp_path))
    cur.record_source(path, rows=1)
    cur.record_source(path, rows=2)
    rows = cur._dbconn.execute(f'SELECT rows FROM {cur.sources_table}').fetchall()
    assert rows == [(2,)]
    assert cur.content_id() is not None


def test_the_registry_is_named_the_same_on_every_sql_backend():
    """The point of the shared name: provenance is looked up the same way
    whichever SQL backend is holding the data."""
    import rotifer.db.sql.clickhouse.core as ch
    assert rdss.BaseSQLite3Cursor._sources_table == 'rotifer_sources'
    assert ch.BaseClickHouseCursor._sources_table == 'rotifer_sources'
    for name in ('sources_table', 'source_version', 'create_sources',
                 'record_source', 'content_id'):
        assert hasattr(rdss.BaseSQLite3Cursor, name), name
        assert hasattr(ch.BaseClickHouseCursor, name), name


if __name__ == '__main__':
    import sys
    failures = 0
    for name, test in sorted(globals().items()):
        if not name.startswith('test_') or not callable(test):
            continue
        try:
            if 'tmp_path' in test.__code__.co_varnames[:test.__code__.co_argcount]:
                with tempfile.TemporaryDirectory() as tmp:
                    test(tmp)
            else:
                test()
            ok, detail = True, ''
        except Exception as exc:
            ok, detail = False, f'{type(exc).__name__}: {exc}'
        failures += not ok
        print('  %-56s %s' % (name, 'ok' if ok else 'FAILED'))
        if detail:
            print('      %s' % detail)
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
