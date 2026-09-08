#!/usr/bin/env python3
"""Tests for loading idmapping.dat into ClickHouse.

`load()` is the path that puts the data there, and it was the one piece of
this package with no coverage at all: a rename once left it calling a class
that no longer existed, and nothing failed until someone tried to load a
release. It needs a server, so these tests skip when there is none rather
than fail.

They never touch the real table. Each builds a mirror of its own -- a
directory laid out the way UniProt's FTP site is, holding a handful of rows --
loads it into a table named after the test, and drops that table afterwards.
The point is the shape of the exchange rather than its volume: a real release
is 90 GB and takes hours, and nothing in that time is exercised that ten rows
do not exercise.

What is checked:

  * both load methods. `python` sends rows through the driver, `client` pipes
    the file through the clickhouse program, and they are different code
  * the rows arrive with the release stamped on them, and can be read back
    through the cursor that will be used to query them
  * a completed load records where it came from, so `content_id()` matches
    what the mirror reports and a delegator can tell the two hold one file
  * a table created but never loaded records nothing, which is the gate that
    keeps a half filled table from being taken for a complete one

Run under pytest, or standalone: python test/db/uniprot/test_clickhouse_load.py
"""

import gzip
import os
import shutil
import uuid

import pytest

from rotifer.db.uniprot import clickhouse as ruch
from rotifer.db.uniprot import mirror as rum

#: Rows written into every synthetic mirror, as idmapping.dat has them.
ROWS = [
    ('TEST0001', 'RefSeq', 'NP_TEST0001.1'),
    ('TEST0001', 'KEGG', 'vg:0001'),
    ('TEST0001', 'EMBL-CDS', 'AAT00001.1'),
    ('TEST0002', 'RefSeq', 'NP_TEST0002.1'),
    ('TEST0002', 'KEGG', 'vg:0002'),
    ('TEST0003', 'RefSeq', 'NP_TEST0003.1'),
]

RELEASE = 'test_release'


def server_available():
    """
    Find whether a ClickHouse server can be reached.

    Returns
    -------
    bool
    """
    try:
        cursor = ruch.MappingCursor(progress=False)
        cursor.query('SELECT 1')
        return True
    except Exception:
        return False


needs_server = pytest.mark.skipif(
    not server_available(),
    reason='needs a reachable ClickHouse server',
)


def build_mirror(root, compressed=False):
    """
    Write a directory laid out like a UniProt mirror.

    Parameters
    ----------
    root : str
        Directory to build it under.
    compressed : bool, default False
        Write the gzip compressed copy instead, which UniProt also
        publishes and a cursor finds just as well.

    Returns
    -------
    str
        The mirror's root, as a cursor expects to be given.
    """
    path = os.path.join(root, 'knowledgebase', 'idmapping')
    os.makedirs(path, exist_ok=True)
    text = "".join("\t".join(row) + "\n" for row in ROWS)
    if compressed:
        with gzip.open(os.path.join(path, 'idmapping.dat.gz'), 'wt') as fh:
            fh.write(text)
    else:
        with open(os.path.join(path, 'idmapping.dat'), 'wt') as fh:
            fh.write(text)
    return root


@pytest.fixture
def loader(tmp_path):
    """A cursor on a table of its own, dropped when the test ends."""
    table = 'test_load_' + uuid.uuid4().hex[:8]
    cursor = ruch.MappingCursor(table=table, release=RELEASE, progress=False)
    try:
        yield cursor, build_mirror(str(tmp_path))
    finally:
        try:
            cursor.command(f'DROP TABLE IF EXISTS {cursor.qualified_name}')
            if cursor.has_table(cursor._sources_table):
                cursor.command(
                    f'ALTER TABLE {cursor.sources_table} DELETE WHERE `table` = %(t)s',
                    parameters={'t': table},
                )
        except Exception:
            pass


@needs_server
@pytest.mark.parametrize('method', ['python', 'client'])
def test_load_puts_every_row_in_the_table(loader, method):
    """Both methods must land the same rows: one sends them through the
    driver, the other pipes the file through the clickhouse program."""
    if method == 'client' and not shutil.which(ruch.config['executable']):
        pytest.skip('needs the clickhouse client program')
    cursor, mirror = loader
    assert cursor.create() is True
    assert cursor.load(mirror, release=RELEASE, method=method) == len(ROWS)
    assert cursor.count() == len(ROWS)


@needs_server
def test_loaded_rows_are_readable_through_the_cursor(loader):
    """A load nobody can query is not a load: read it back the way the
    package will, rather than by counting rows."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    frame = cursor.fetchall(['TEST0001'], source=cursor.UNIPROTKB)
    assert sorted(zip(frame.target_type, frame.target)) == [
        ('EMBL-CDS', 'AAT00001.1'),
        ('KEGG', 'vg:0001'),
        ('RefSeq', 'NP_TEST0001.1'),
    ]


@needs_server
def test_the_release_is_stamped_on_every_row(loader):
    """The table is partitioned by release, so a row without one would be
    unreachable to a cursor filtering by it and undeletable as a partition."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    releases = cursor.query(f'SELECT DISTINCT release FROM {cursor.qualified_name}')
    assert releases.release.tolist() == [RELEASE]


@needs_server
def test_a_completed_load_records_where_it_came_from(loader):
    """The record is what lets a delegator see that the table and the mirror
    hold one file, and skip the ninety second scan."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    reader = rum.MappingCursor(path=mirror, progress=False)
    assert cursor.content_id() is not None
    assert cursor.content_id() == reader.content_id()


@needs_server
def test_a_table_that_was_never_loaded_records_nothing(loader):
    """The completeness gate: creating the table is not loading it, and an
    empty table must not be mistaken for a copy of the file."""
    cursor, mirror = loader
    cursor.create()
    assert cursor.count() == 0
    assert cursor.content_id() is None


@needs_server
def test_load_accepts_a_mirror_cursor_as_well_as_a_path(loader):
    """The signature takes either, and the branch that unwraps a cursor is
    exactly the one a rename once broke."""
    cursor, mirror = loader
    cursor.create()
    reader = rum.MappingCursor(path=mirror, progress=False)
    assert cursor.load(reader, release=RELEASE, method='python') == len(ROWS)


@needs_server
def test_loading_a_mirror_with_no_file_changes_nothing(loader, tmp_path):
    """An empty directory is not a mirror. Saying so by returning the row
    count leaves the table as it was, rather than half filled."""
    cursor, mirror = loader
    cursor.create()
    empty = str(tmp_path / 'not-a-mirror')
    os.makedirs(empty, exist_ok=True)
    assert cursor.load(empty, release=RELEASE, method='python') == 0
    assert cursor.content_id() is None


@needs_server
@pytest.mark.parametrize('method', ['python', 'client'])
def test_a_compressed_mirror_loads_the_same_rows(loader, tmp_path, method):
    """UniProt publishes the file gzipped, and a mirror may hold it that way.
    The client path decompresses it in a separate process, which is the half
    of that code a plain file never reaches."""
    if method == 'client' and not shutil.which(ruch.config['executable']):
        pytest.skip('needs the clickhouse client program')
    cursor, _ = loader
    gzipped = build_mirror(str(tmp_path / 'gz'), compressed=True)
    cursor.create()
    assert cursor.load(gzipped, release=RELEASE, method=method) == len(ROWS)
    frame = cursor.fetchall(['TEST0002'], source=cursor.UNIPROTKB)
    assert sorted(frame.target_type) == ['KEGG', 'RefSeq']


@needs_server
@pytest.mark.parametrize('method', ['python', 'client'])
def test_load_stamps_the_release_it_was_given(loader, method):
    """Regression: the two methods disagreed. The client path put the release
    it was told into the statement, while the python path let insert() drop
    the column and stamp the cursor's own release instead, so loading a new
    release into an existing cursor silently filled the wrong partition."""
    if method == 'client' and not shutil.which(ruch.config['executable']):
        pytest.skip('needs the clickhouse client program')
    cursor, mirror = loader          # the cursor's own release is RELEASE
    cursor.create()
    cursor.load(mirror, release='another_release', method=method)
    releases = cursor.query(f'SELECT DISTINCT release FROM {cursor.qualified_name}')
    assert releases.release.tolist() == ['another_release']


# --------------------------------------------------------- table lifecycle

@needs_server
def test_create_is_idempotent(loader):
    """Creating a table that is already there must not empty it: load() calls
    create() whenever the table is missing, and a caller may call it too."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    assert cursor.create() is True
    assert cursor.count() == len(ROWS)


@needs_server
def test_create_replace_discards_the_data(loader):
    """replace=True drops the table first, so every row goes. That is the
    point of it, and the reason it is not the default."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    assert cursor.count() == len(ROWS)
    assert cursor.create(replace=True) is True
    assert cursor.count() == 0


@needs_server
def test_create_applies_the_whole_schema(loader):
    """A bare table would answer forward lookups and be hopeless at reverse
    ones: the by_id projection is a second copy of the data sorted the other
    way, and is what makes them cheap. Its absence would not fail a query,
    only make it scan, so nothing else would notice."""
    cursor, mirror = loader
    cursor.create()
    schema = cursor.schema
    assert 'PROJECTION by_id' in schema
    assert 'PARTITION BY release' in schema


@needs_server
def test_create_survives_semicolons_inside_comments(loader):
    """The schema is split on semicolons, and its comments contain statements
    that end in one. Splitting first would cut a comment in half and send the
    remainder to the server as SQL."""
    cursor, mirror = loader
    assert cursor.create() is True
    assert cursor.has_table() is True


@needs_server
def test_drop_release_removes_only_that_release(loader):
    """The table holds several releases at once, one partition each, so
    retiring an old one must leave the current one untouched."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release='old_release', method='python')
    cursor.load(mirror, release=RELEASE, method='python')
    assert cursor.count() == 2 * len(ROWS)

    cursor.drop_release('old_release')
    assert cursor.count() == len(ROWS)
    releases = cursor.query(f'SELECT DISTINCT release FROM {cursor.qualified_name}')
    assert releases.release.tolist() == [RELEASE]


@needs_server
def test_dropping_a_release_that_is_not_there_changes_nothing(loader):
    """Retiring a release the table never held is not an error: it is already
    in the state the caller asked for."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    cursor.drop_release('no_such_release')
    assert cursor.count() == len(ROWS)


@needs_server
def test_a_dropped_release_no_longer_claims_to_hold_the_file(loader):
    """The load registry says this table holds a copy of that file, which is
    what lets a delegator skip the mirror. Dropping the rows must withdraw the
    claim, or the mirror would be skipped in favour of a partition that is no
    longer there."""
    cursor, mirror = loader
    cursor.create()
    cursor.load(mirror, release=RELEASE, method='python')
    assert cursor.content_id() is not None

    cursor.drop_release(RELEASE)
    assert cursor.count() == 0
    assert cursor.content_id() is None


if __name__ == '__main__':
    import sys
    if not server_available():
        print('no ClickHouse server reachable: nothing to do')
        sys.exit(0)
    sys.exit(pytest.main([__file__, '-q']))
