#!/usr/bin/env python3
"""Tests for rotifer.graph.network.adjacency_edges.

What has to hold, and why:

  * each label is joined to the next one in the same block, never to a
    label of another block, and the input must already be in genome order
  * a multi-domain protein 'A+B' gives a fusion edge A->B, and its last
    domain a neighbour edge to the first label of the next protein; the
    rule for a fusion is that both labels come from rows with the same pid
  * target may be a list: the first non-missing column is the label, and
    '?' when all are missing
  * ignore drops labels (the ones around them become adjacent); replace
    renames them, and replace=None is a no-op, not an error
  * rows whose type is not in ftype are dropped before edges are made
  * orient=True writes neighbour pairs smaller label first and leaves
    fusions in N to C order
  * interaction=False merges fusion and neighbour counts; stats replaces
    the default weight/blocks aggregation
  * a missing required column raises KeyError, the input is not modified,
    and no edges give an empty frame with the right columns

The data are tiny hand-built tables; nothing touches the network.

Run under pytest: PYTHONPATH=lib pytest test/graph/test_adjacency_edges.py
"""

import pandas as pd
import pytest

pytest.importorskip('networkx')
pytest.importorskip('scipy')

from rotifer.graph.network import adjacency_edges


def frame(labels, blocks=None, pids=None, types=None, col='pfam'):
    n = len(labels)
    return pd.DataFrame({
        'block_id': blocks or ['b1'] * n,
        'pid': pids or [f'p{i}' for i in range(n)],
        'type': types or ['CDS'] * n,
        col: labels})


def rows(edges):
    return [tuple(r) for r in edges.itertuples(index=False)]


def test_single_block_chain():
    e = adjacency_edges(frame(['A', 'B', 'C']))
    assert list(e.columns) == ['source', 'target', 'interaction', 'weight', 'blocks']
    assert rows(e) == [('A', 'B', 'neighbor', 1, 1), ('B', 'C', 'neighbor', 1, 1)]


def test_no_edge_across_blocks():
    e = adjacency_edges(frame(['A', 'B', 'C', 'D'], blocks=['b1', 'b1', 'b2', 'b2']))
    assert rows(e) == [('A', 'B', 'neighbor', 1, 1), ('C', 'D', 'neighbor', 1, 1)]


def test_counts_rows_and_blocks():
    e = adjacency_edges(frame(['A', 'B', 'A', 'B', 'A', 'B'],
                              blocks=['b1'] * 4 + ['b2'] * 2))
    assert rows(e) == [('A', 'B', 'neighbor', 3, 2), ('B', 'A', 'neighbor', 1, 1)]


def test_multidomain_protein_fusion_then_neighbor():
    e = adjacency_edges(frame(['A+B', 'C']))
    assert rows(e) == [('A', 'B', 'fusion', 1, 1), ('B', 'C', 'neighbor', 1, 1)]


def test_fusion_needs_same_pid():
    e = adjacency_edges(frame(['A', 'B'], pids=['p1', 'p1']))
    assert rows(e) == [('A', 'B', 'fusion', 1, 1)]


def test_separator_none_keeps_labels_whole():
    e = adjacency_edges(frame(['A+B', 'C']), separator=None)
    assert rows(e) == [('A+B', 'C', 'neighbor', 1, 1)]


def test_separator_regex():
    e = adjacency_edges(frame(['A|B', 'C']), separator='[|]', separator_regex=True)
    assert rows(e) == [('A', 'B', 'fusion', 1, 1), ('B', 'C', 'neighbor', 1, 1)]


def test_target_list_falls_back_to_second_column():
    df = frame(['A', None, None])
    df['cdd'] = [None, 'X', None]
    e = adjacency_edges(df, target=['pfam', 'cdd'])
    assert rows(e) == [('A', 'X', 'neighbor', 1, 1), ('X', '?', 'neighbor', 1, 1)]


def test_missing_label_is_question_mark():
    e = adjacency_edges(frame(['A', None]))
    assert rows(e) == [('A', '?', 'neighbor', 1, 1)]


def test_ignore_makes_neighbors_adjacent():
    e = adjacency_edges(frame(['A', 'X', 'B']), ignore=['X'])
    assert rows(e) == [('A', 'B', 'neighbor', 1, 1)]


def test_ignore_string_is_one_label():
    e = adjacency_edges(frame(['A', 'N', 'NA', 'B']), ignore='NA')
    assert rows(e) == [('A', 'N', 'neighbor', 1, 1), ('N', 'B', 'neighbor', 1, 1)]


def test_integer_labels_keep_dtype():
    e = adjacency_edges(frame([1, 2, 3]), separator=None)
    assert e['target'].dtype == e['source'].dtype
    assert all(isinstance(v, int) for v in e['target'])
    assert rows(e) == [(1, 2, 'neighbor', 1, 1), (2, 3, 'neighbor', 1, 1)]


def test_replace_mapping_then_ignore():
    e = adjacency_edges(frame(['A', 'B', 'C']), replace={'B': 'Z'}, ignore=['Z'])
    assert rows(e) == [('A', 'C', 'neighbor', 1, 1)]


def test_replace_regex():
    e = adjacency_edges(frame(['A1', 'A2']), replace={r'\d': ''}, replace_regex=True)
    assert rows(e) == [('A', 'A', 'neighbor', 1, 1)]


def test_replace_none_is_noop():
    df = frame(['A', 'B'])
    assert rows(adjacency_edges(df, replace=None)) == rows(adjacency_edges(df))


def test_non_cds_rows_filtered():
    df = frame(['A', 'R', 'B'], types=['CDS', 'rRNA', 'CDS'])
    assert rows(adjacency_edges(df)) == [('A', 'B', 'neighbor', 1, 1)]
    assert rows(adjacency_edges(df, ftype=['CDS', 'rRNA'])) == [
        ('A', 'R', 'neighbor', 1, 1), ('R', 'B', 'neighbor', 1, 1)]


def test_orient_reorders_neighbors_only():
    df = frame(['B', 'A', 'D+C'], pids=['p1', 'p2', 'p3'])
    assert rows(adjacency_edges(df)) == [
        ('A', 'D', 'neighbor', 1, 1),
        ('B', 'A', 'neighbor', 1, 1),
        ('D', 'C', 'fusion', 1, 1)]
    assert rows(adjacency_edges(df, orient=True)) == [
        ('A', 'B', 'neighbor', 1, 1),
        ('A', 'D', 'neighbor', 1, 1),
        ('D', 'C', 'fusion', 1, 1)]


def test_orient_merges_opposite_neighbors():
    e = adjacency_edges(frame(['A', 'B', 'A'], blocks=['b1', 'b1', 'b1']), orient=True)
    assert rows(e) == [('A', 'B', 'neighbor', 2, 1)]


def test_interaction_false():
    e = adjacency_edges(frame(['A+B', 'A', 'B'], blocks=['b1', 'b1', 'b1']),
                        interaction=False)
    assert list(e.columns) == ['source', 'target', 'weight', 'blocks']
    assert rows(e) == [('A', 'B', 2, 1), ('B', 'A', 1, 1)]


def test_custom_stats():
    e = adjacency_edges(frame(['A', 'B']), stats={'n': ('pid', 'nunique')})
    assert list(e.columns) == ['source', 'target', 'interaction', 'n']
    assert rows(e) == [('A', 'B', 'neighbor', 1)]


def test_missing_column_raises_keyerror():
    with pytest.raises(KeyError, match='pid'):
        adjacency_edges(frame(['A', 'B']).drop(columns='pid'))
    with pytest.raises(KeyError, match='nope'):
        adjacency_edges(frame(['A', 'B']), target=['pfam', 'nope'])


def test_input_not_mutated():
    df = frame(['A+B', None, 'C'])
    before = df.copy()
    adjacency_edges(df, ignore=['C'], replace={'A': 'Z'}, orient=True)
    pd.testing.assert_frame_equal(df, before)


def test_empty_result_has_columns():
    e = adjacency_edges(frame(['A']))
    assert e.empty
    assert list(e.columns) == ['source', 'target', 'interaction', 'weight', 'blocks']
    e = adjacency_edges(frame([]).iloc[0:0])
    assert e.empty
    assert list(e.columns) == ['source', 'target', 'interaction', 'weight', 'blocks']
