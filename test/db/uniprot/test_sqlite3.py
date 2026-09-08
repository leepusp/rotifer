#!/usr/bin/env python3
"""Tests for the SQLite3 UniProt mapping backend.

This backend answers the same query as the ClickHouse one, from a file rather
than a server, so it needs nothing to run: every test here builds its own
database in a temporary directory.

That matters more than convenience. The SQL both backends run is built once,
in rotifer.db.sql.mapping, and only ClickHouse is exercised against real data;
if the shared builder drifts, this is where it shows. In particular SQLite3
binds by position, so a fragment built out of order pairs the wrong value with
the wrong placeholder here while ClickHouse, binding by name, goes on working.

What is checked:

  * every query shape, including the two that must not join
  * the asymmetry between the ends -- an omitted source matches an accession
    too, an omitted target does not return one
  * batching, since SQLite3 limits how many values a statement may bind and a
    long query is therefore asked in pieces
  * loading from a mirror, and the rule that a filtered load records no
    provenance: a file holding part of a release is not a copy of it, and a
    delegator must not skip the mirror on its word

Run under pytest, or standalone: python test/db/uniprot/test_sqlite3.py
"""

import os

import pandas as pd
import pytest

from rotifer.db.uniprot import mirror as rum
from rotifer.db.uniprot import sqlite3 as rus

U = rus.MappingCursor.UNIPROTKB

#: The rows every fixture starts from, as idmapping.dat has them.
ROWS = [
    ('Q6GZX4', 'RefSeq', 'YP_031579.1'),
    ('Q6GZX4', 'KEGG', 'vg:2947773'),
    ('Q6GZX4', 'EMBL-CDS', 'AAT09660.1'),
    ('Q6GZX4', 'UniRef100', 'UniRef100_Q6GZX4'),
    ('P00750', 'RefSeq', 'NP_000921.1'),
    ('P00750', 'KEGG', 'hsa:5327'),
]


def cursor(tmp_path, rows=ROWS, release=None, **kwargs):
    """A cursor on a database of its own, holding `rows`."""
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'),
                           release=release, progress=False, **kwargs)
    mc.create()
    if rows:
        mc.insert(pd.DataFrame(list(rows), columns=['accession', 'id_type', 'id']),
                  release=release)
    return mc


def build_mirror(root, rows=ROWS):
    """A directory laid out like a UniProt mirror."""
    path = os.path.join(root, 'knowledgebase', 'idmapping')
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, 'idmapping.dat'), 'wt') as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return root


def pairs(frame):
    return sorted(set(zip(frame.target_type, frame.target)))


# ----------------------------------------------------------- the table

def test_create_makes_the_table_and_both_indexes(tmp_path):
    """Two indexes, because the two directions read different columns. The
    one on id is what a ClickHouse projection is for: without it a reverse
    lookup reads the whole table, which fails nothing and is merely slow."""
    mc = cursor(tmp_path, rows=[])
    names = {r[0] for r in mc._dbconn.execute(
        "SELECT name FROM sqlite_master WHERE type='index'").fetchall()}
    assert 'idmapping_by_accession' in names
    assert 'idmapping_by_id' in names


def test_create_is_idempotent(tmp_path):
    mc = cursor(tmp_path)
    assert mc.count() == len(ROWS)
    assert mc.create() is True
    assert mc.count() == len(ROWS)


def test_create_replace_empties_the_table(tmp_path):
    mc = cursor(tmp_path)
    assert mc.create(replace=True) is True
    assert mc.count() == 0


# ----------------------------------------------------- the query shapes

def test_an_accession_returns_its_cross_references(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4'], source=U)
    assert pairs(frame) == [
        ('EMBL-CDS', 'AAT09660.1'), ('KEGG', 'vg:2947773'),
        ('RefSeq', 'YP_031579.1'), ('UniRef100', 'UniRef100_Q6GZX4')]


def test_a_named_target_narrows_the_answer(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4'], source=U, target=['RefSeq', 'KEGG'])
    assert pairs(frame) == [('KEGG', 'vg:2947773'), ('RefSeq', 'YP_031579.1')]


def test_an_identifier_finds_its_accession(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['YP_031579.1'], source=['RefSeq'], target=U)
    assert frame.to_dict('records') == [{
        'source': 'YP_031579.1', 'source_type': 'RefSeq', 'accession': 'Q6GZX4',
        'target': 'Q6GZX4', 'target_type': 'UniProtKB-AC'}]


def test_one_database_translates_into_another(tmp_path):
    """The shape that needs a join, and needs the far side restricted to the
    accessions the near side found."""
    mc = cursor(tmp_path)
    frame = mc.fetchall(['AAT09660.1'], source=['EMBL-CDS'], target=['RefSeq'])
    assert pairs(frame) == [('RefSeq', 'YP_031579.1')]


def test_both_ends_open(tmp_path):
    """With neither end named the identifier may be anything and every
    database is wanted, which is the broadest question this answers."""
    mc = cursor(tmp_path)
    frame = mc.fetchall(['YP_031579.1'])
    assert set(frame.source_type) == {'RefSeq'}
    assert pairs(frame) == [
        ('EMBL-CDS', 'AAT09660.1'), ('KEGG', 'vg:2947773'),
        ('RefSeq', 'YP_031579.1'), ('UniRef100', 'UniRef100_Q6GZX4')]


def test_an_open_source_also_matches_an_accession(tmp_path):
    """Nothing said what the identifier is, so it may be an accession."""
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4'], target=['RefSeq'])
    assert (frame.source_type == U).any()


def test_an_open_target_does_not_return_the_accession(tmp_path):
    """'Every database' means the ones cross-referenced. The accession is the
    table's key rather than one of them."""
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4'], source=U)
    assert U not in set(frame.target_type)


def test_naming_the_accession_as_a_target_returns_it(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4'], source=U, target=['RefSeq', U])
    assert (U, 'Q6GZX4') in pairs(frame)


def test_several_identifiers_at_once(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['Q6GZX4', 'P00750'], source=U, target=['RefSeq'])
    assert sorted(set(frame.source)) == ['P00750', 'Q6GZX4']


def test_an_unknown_identifier_is_registered(tmp_path):
    mc = cursor(tmp_path)
    frame = mc.fetchall(['NO_SUCH_ID'], source=U)
    assert frame.empty
    assert 'NO_SUCH_ID' in mc.missing_ids()


def test_a_missing_table_is_reported_rather_than_raised(tmp_path):
    """A database that was never built is a query this backend cannot answer,
    which a delegator should hear about and pass on."""
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'empty.sqlite3'), progress=False)
    frame = mc.fetchall(['Q6GZX4'], source=U)
    assert frame.empty and 'Q6GZX4' in mc.missing_ids()


# --------------------------------------------------------- batching

def test_a_query_longer_than_a_batch_still_answers(tmp_path):
    """SQLite3 limits how many values one statement may bind, so a long query
    is asked in pieces and the pieces are put back together."""
    mc = cursor(tmp_path, batch_size=2)
    frame = mc.fetchall(['Q6GZX4', 'P00750', 'NO_SUCH_ID'], source=U, target=['RefSeq'])
    assert sorted(set(frame.source)) == ['P00750', 'Q6GZX4']
    assert 'NO_SUCH_ID' in mc.missing_ids()


def test_batching_does_not_repeat_rows(tmp_path):
    mc = cursor(tmp_path, batch_size=1)
    frame = mc.fetchall(['Q6GZX4', 'P00750'], source=U)
    assert frame.duplicated().sum() == 0


# --------------------------------------------------------- releases

def test_a_release_filter_selects_its_own_rows(tmp_path):
    """A file may hold more than one release, and a cursor asking for one must
    not be answered from another."""
    mc = cursor(tmp_path, rows=[], release='2026_01')
    mc.insert(pd.DataFrame([('Q6GZX4', 'RefSeq', 'OLD.1')],
                           columns=['accession', 'id_type', 'id']), release='2025_01')
    mc.insert(pd.DataFrame([('Q6GZX4', 'RefSeq', 'NEW.1')],
                           columns=['accession', 'id_type', 'id']), release='2026_01')
    frame = mc.fetchall(['Q6GZX4'], source=U)
    assert pairs(frame) == [('RefSeq', 'NEW.1')]


# ------------------------------------------------------------ loading

def test_load_fills_the_table_from_a_mirror(tmp_path):
    mirror = build_mirror(str(tmp_path / 'mirror'))
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'),
                           release='2026_01', progress=False)
    assert mc.load(mirror, release='2026_01') == len(ROWS)
    assert pairs(mc.fetchall(['Q6GZX4'], source=U, target=['RefSeq'])) == \
        [('RefSeq', 'YP_031579.1')]


def test_a_whole_load_records_where_it_came_from(tmp_path):
    """The record is what lets a delegator see that this file and the mirror
    hold the same thing."""
    mirror = build_mirror(str(tmp_path / 'mirror'))
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'),
                           release='2026_01', progress=False)
    mc.load(mirror, release='2026_01')
    reader = rum.MappingCursor(path=mirror, progress=False)
    assert mc.content_id() == reader.content_id()


def test_a_filtered_load_records_nothing(tmp_path):
    """A file holding part of a release is not a copy of it. Recording it
    would have a delegator skip the mirror for the databases this file was
    never given."""
    mirror = build_mirror(str(tmp_path / 'mirror'))
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'),
                           release='2026_01', progress=False)
    rows = mc.load(mirror, release='2026_01', id_type=['RefSeq'])
    assert rows == 2                      # only the RefSeq rows
    assert mc.content_id() is None
    assert pairs(mc.fetchall(['Q6GZX4'], source=U)) == [('RefSeq', 'YP_031579.1')]


def test_loading_a_directory_that_is_not_a_mirror_changes_nothing(tmp_path):
    empty = str(tmp_path / 'not-a-mirror')
    os.makedirs(empty, exist_ok=True)
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'), progress=False)
    mc.create()
    assert mc.load(empty) == 0
    assert mc.content_id() is None


def test_load_accepts_a_mirror_cursor(tmp_path):
    mirror = build_mirror(str(tmp_path / 'mirror'))
    reader = rum.MappingCursor(path=mirror, progress=False)
    mc = rus.MappingCursor(os.path.join(str(tmp_path), 'uniprot.sqlite3'),
                           release='2026_01', progress=False)
    assert mc.load(reader, release='2026_01') == len(ROWS)


# -------------------------------------------------- the shared contract

def test_the_columns_are_the_ones_every_backend_returns(tmp_path):
    """All three backends return this, which is what lets a delegator stack
    them."""
    mc = cursor(tmp_path)
    assert mc.fetchall(['Q6GZX4'], source=U).columns.tolist() == [
        'source', 'source_type', 'accession', 'target', 'target_type']
    assert mc.empty().columns.tolist() == mc.columns


def test_tables_are_declared_by_the_class(tmp_path):
    """As on the server backend: which tables a cursor reads is a property of
    its SQL, not something a caller chooses."""
    mc = cursor(tmp_path)
    assert mc.tables == {'mapping': 'idmapping'}
    assert mc.qualified() == 'idmapping'      # a file needs no qualifying


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
