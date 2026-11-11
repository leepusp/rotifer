#!/usr/bin/env python3
"""Tests for rotifer.genome.data.NeighborhoodDF.compact_neighborhood (alias compact).

compact() groups a neighborhood table by block_id (or by strand_group when
strand=True) and renders each group as a string of annotations and arrows,
oriented so the query always points forward. The groups must keep their key
column: pandas 3 stopped passing it to groupby().apply() callbacks, which broke
compact() with KeyError: 'block_id' (#129, #140).

Fixture: two blocks.
  block 1: A(+), *B(+) query, ?(-) unannotated   -> strand change before the last locus
  block 2: *C(-) query, D(-)                       -> reversed so the query points forward

Run under pytest, or standalone: python test/genome/test_compact.py
"""

import pandas as pd

from rotifer.genome.data import NeighborhoodDF

EXTRA = ['lineage', 'organism', 'seq_type', 'classification']
COMPACT_COLUMNS = ['compact', 'lineage', 'organism', 'seq_type', 'classification', 'pid_compact']


def neighborhood(rows=None):
    if rows is None:
        rows = [
            dict(block_id=1, pid='p1', cluster='A',  strand=1,  query=0),
            dict(block_id=1, pid='p2', cluster='B',  strand=1,  query=1),
            dict(block_id=1, pid='p3', cluster=None, strand=-1, query=0),
            dict(block_id=2, pid='p4', cluster='C',  strand=-1, query=1),
            dict(block_id=2, pid='p5', cluster='D',  strand=-1, query=0),
        ]
    df = pd.DataFrame(rows, columns=['block_id', 'pid', 'cluster', 'strand', 'query'])
    for c in EXTRA:
        df[c] = 'x'
    df['type'] = 'CDS'
    return NeighborhoodDF(df)


def compact(df, **kwargs):
    return df.compact(annotation='cluster', **kwargs)


def test_whole_blocks():
    r = compact(neighborhood())
    assert r.index.name == 'block_id'
    assert list(r.columns) == COMPACT_COLUMNS
    assert r.compact.to_dict() == {1: 'A->*B->|<-?', 2: 'D->*C->'}
    assert r.pid_compact.to_dict() == {1: 'p1->*p2->|<-p3', 2: 'p5->*p4->'}


def test_same_strand_only():
    r = compact(neighborhood(), strand=True)
    assert r.index.name == 'block_id'
    assert list(r.columns) == COMPACT_COLUMNS
    assert r.compact.to_dict() == {1: 'A->*B->', 2: 'D->*C->'}
    assert r.pid_compact.to_dict() == {1: 'p1->*p2->', 2: 'p5->*p4->'}


def test_fill_unknown_from_column():
    r = compact(neighborhood(), fill_unk='pid')
    assert r.compact.to_dict() == {1: 'A->*B->|<-p3', 2: 'D->*C->'}


def test_extra_columns_carried_over():
    r = compact(neighborhood())
    for c in EXTRA:
        assert r[c].to_list() == ['x', 'x']


def test_empty_input_matches_normal_schema():
    normal = compact(neighborhood())
    for strand in (False, True):
        r = compact(neighborhood(rows=[]), strand=strand)
        assert len(r) == 0, strand
        assert type(r) is type(normal), strand
        assert r.index.name == normal.index.name, strand
        assert list(r.columns) == list(normal.columns), strand


if __name__ == '__main__':
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith('test_') and callable(f)]
    failed = 0
    for name, test in tests:
        try:
            test()
            print(f'ok    {name}')
        except Exception as e:
            failed += 1
            print(f'FAIL  {name}: {type(e).__name__}: {e}')
    print(f'{len(tests) - failed}/{len(tests)} passed')
    raise SystemExit(1 if failed else 0)
