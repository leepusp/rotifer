#!/usr/bin/env python3
"""Tests for NeighborhoodDF.select_neighbors (step 5: each row once, any index)."""

import pandas as pd
import pytest

from rotifer.genome.data import NeighborhoodDF


def make_table(queries=(2, 3, 7), index=None):
    """Ten CDS rows in two blocks: p0..p5 in block A, p6..p9 in block B.

    Rows listed in `queries` get query == 1. Strand alternates + + - - + +
    for block A and + - - + for block B. `index` replaces the default
    0..9 index labels.
    """
    df = NeighborhoodDF({
        'block_id': ['A'] * 6 + ['B'] * 4,
        'pid': [f'p{i}' for i in range(10)],
        'type': ['CDS'] * 10,
        'strand': [1, 1, -1, -1, 1, 1, 1, -1, -1, 1],
        'query': [1 if i in queries else 0 for i in range(10)],
    })
    if index is not None:
        df.index = index
    return df


def pids(df):
    """The pid column of `df` as a list."""
    return df['pid'].tolist()


def test_overlapping_windows_give_each_row_once():
    """Adjacent targets 2 and 3, one row each side: p1..p4 once, not twice."""
    out = make_table().select_neighbors(before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']
    assert out['new_query'].tolist() == [0, 1, 1, 0, 0, 1, 0]


def test_wide_overlap_is_not_duplicated():
    """Two rows each side of 2 and 3 cover p0..p5 once (10 rows before the fix)."""
    out = make_table(queries=(2, 3)).select_neighbors(before=2, after=2)
    assert pids(out) == ['p0', 'p1', 'p2', 'p3', 'p4', 'p5']


def test_shifted_index_keeps_labels():
    """A non-RangeIndex works and its labels are kept."""
    labels = [50, 40, 30, 20, 10, 0, 90, 80, 70, 60]
    out = make_table(index=labels).select_neighbors(before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']
    assert out.index.tolist() == [40, 30, 20, 10, 90, 80, 70]


def test_windows_are_clipped_at_block_edges():
    """Targets on the first row of A and last row of B stay inside their block."""
    out = make_table(queries=(0, 9)).select_neighbors(before=2, after=2)
    assert pids(out) == ['p0', 'p1', 'p2', 'p7', 'p8', 'p9']


def test_two_blocks_are_independent():
    """A window never reaches into the next block: p5 and p6 are neighbours of nobody."""
    out = make_table(queries=(5, 6)).select_neighbors(before=1, after=1)
    assert pids(out) == ['p4', 'p5', 'p6', 'p7']


def test_before_counts_rows_below_and_after_rows_above():
    """Old convention, kept: before=2, after=0 on row 3 selects p3..p5."""
    out = make_table(queries=(3,)).select_neighbors(before=2, after=0)
    assert pids(out) == ['p3', 'p4', 'p5']


def test_query_column_is_not_changed_and_new_query_is_added():
    """The query column is left as it was; new_query marks the targets."""
    out = make_table(queries=(3,)).select_neighbors(before=1, after=1)
    assert out['query'].tolist() == [0, 1, 0]
    assert out['new_query'].tolist() == [0, 1, 0]


def test_value_counts_full_names_the_count_columns():
    """Columns are <column>, block_id_count, pid_count, relative_block_id, relative_pid.

    Group x has 6 CDS in both blocks, group y has 4 CDS in block A only; the
    table has 2 blocks and 10 distinct pids.
    """
    df = make_table()
    df['group'] = ['x', 'x', 'y', 'y', 'y', 'y', 'x', 'x', 'x', 'x']
    out = df.value_counts_full('group')
    expected = pd.DataFrame({
        'group': ['x', 'y'],
        'block_id_count': [2, 1],
        'pid_count': [6, 4],
        'relative_block_id': [100.0, 50.0],
        'relative_pid': [60.0, 40.0],
    })
    pd.testing.assert_frame_equal(out, expected, check_dtype=False)
