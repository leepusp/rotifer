#!/usr/bin/env python3
"""Tests for rotifer.genome.data.NeighborhoodDF.to_network.

to_network() lists the edges between adjacent features of each block. With
strand=True the blocks are first cut to the query strand and read in the
query's direction; with strand=False the table is used as it is and neighbour
edges are undirected (orient=True in rotifer.graph.network.adjacency_edges).

Fixture: two blocks, one nucleotide each.
  block 1, plus strand : X, A+B (query, a fused protein), C
  block 2, minus strand: C, A (query), X in genome order, so it reads X, A, C
                         in the query's direction

Run under pytest: PYTHONPATH=lib pytest test/genome/test_neighborhood_network.py
"""

import sys

import pandas as pd
import pytest

from rotifer.genome.data import NeighborhoodDF

pytest.importorskip('networkx')


def neighborhood():
    rows = [
        (1, 'p1', 1, 1, 0, 'X'), (1, 'p2', 2, 1, 1, 'A+B'), (1, 'p3', 3, 1, 0, 'C'),
        (2, 'p4', 1, -1, 0, 'C'), (2, 'p5', 2, -1, 1, 'A'), (2, 'p6', 3, -1, 0, 'X'),
    ]
    df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'strand', 'query', 'pfam'])
    df['nucleotide'] = ['n1'] * 3 + ['n2'] * 3
    df['end'] = df.start + 1
    df['type'] = 'CDS'
    df['assembly'] = 'a'
    df['internal_id'] = range(6)
    df['feature_order'] = [1, 2, 3, 1, 2, 3]
    df['topology'] = 'linear'
    return NeighborhoodDF(df)


def edges(rows):
    return pd.DataFrame(rows, columns=['source', 'target', 'interaction', 'weight', 'blocks'])


def test_strand_reads_blocks_in_query_direction():
    got = neighborhood().to_network()
    expected = edges([
        ('A', 'B', 'fusion', 1, 1),
        ('A', 'C', 'neighbor', 1, 1),
        ('B', 'C', 'neighbor', 1, 1),
        ('X', 'A', 'neighbor', 2, 2),
    ])
    pd.testing.assert_frame_equal(got, expected, check_dtype=False)


def test_no_strand_orients_neighbor_pairs():
    got = neighborhood().to_network(strand=False)
    expected = edges([
        ('A', 'B', 'fusion', 1, 1),
        ('A', 'C', 'neighbor', 1, 1),
        ('A', 'X', 'neighbor', 2, 2),
        ('B', 'C', 'neighbor', 1, 1),
    ])
    pd.testing.assert_frame_equal(got, expected, check_dtype=False)


def test_does_not_modify_the_table():
    df = neighborhood()
    before = df.copy()
    df.to_network()
    df.to_network(strand=False)
    pd.testing.assert_frame_equal(pd.DataFrame(df), pd.DataFrame(before))


def test_without_the_graph_module_the_error_names_the_extra(monkeypatch):
    monkeypatch.setitem(sys.modules, 'rotifer.graph.network', None)
    with pytest.raises(ImportError, match=r'pip install rotifer\[graph\]') as info:
        neighborhood().to_network()
    assert isinstance(info.value.__cause__, ImportError)
