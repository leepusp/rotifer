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


def test_only_query_1_is_a_query():
    # B (query=2) and C (query=NaN) are neither marked nor used for orientation,
    # and with strand=True C's run has no query, so it is left out.
    df = neighborhood([
        dict(block_id=1, pid='p1', cluster='A', strand=1,  query=1),
        dict(block_id=1, pid='p2', cluster='B', strand=1,  query=2),
        dict(block_id=1, pid='p3', cluster='C', strand=-1, query=None),
    ])
    assert compact(df).compact.to_dict() == {1: '*A->B->|<-C'}
    assert compact(df, strand=True).compact.to_dict() == {1: '*A->B->'}


def test_missing_query_names_the_blocks():
    rows = [dict(block_id=b, pid=f'p{b}', cluster='A', strand=1, query=int(b == 1))
            for b in range(1, 14)]
    try:
        compact(neighborhood(rows))
    except NameError as e:
        assert str(e) == ('Query_Missing: no query in block_id(s) '
                          '2, 3, 4, 5, 6, 7, 8, 9, 10, 11 and 2 more'), str(e)
    else:
        raise AssertionError('expected NameError')
    try:
        compact(neighborhood(rows[:3]))
    except NameError as e:
        assert str(e) == 'Query_Missing: no query in block_id(s) 2, 3', str(e)
    else:
        raise AssertionError('expected NameError')


def test_split_block_gives_the_same_lines():
    # Block 1's rows are split by block 2's. A and B are on the same strand,
    # so there must be no '|' between them, and strand=True keeps them in
    # one row.
    block1 = [dict(block_id=1, pid='p1', cluster='A', strand=1,  query=1),
              dict(block_id=1, pid='p2', cluster='B', strand=1,  query=0),
              dict(block_id=1, pid='p3', cluster='C', strand=-1, query=0)]
    block2 = [dict(block_id=2, pid='p4', cluster='D', strand=-1, query=1)]
    contiguous = neighborhood(block1 + block2)
    split = neighborhood(block1[:1] + block2 + block1[1:])
    for strand in (False, True):
        expected = compact(contiguous, strand=strand)
        actual = compact(split, strand=strand)
        assert actual.sort_index().equals(expected.sort_index()), (strand, actual, expected)
    assert compact(split).compact.to_dict() == {1: '*A->B->|<-C', 2: '*D->'}
    assert compact(split, strand=True).compact.to_dict() == {1: '*A->B->', 2: '*D->'}


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
