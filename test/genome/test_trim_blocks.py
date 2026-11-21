#!/usr/bin/env python3
"""Tests for NeighborhoodDF.trim_blocks, its shim select_neighbors and value_counts_full.

What has to hold, and why:

  * trim_blocks keeps, for each target, `before` rows above and `after` rows
    below it inside the same block_id, counting every row, and drops blocks
    without a target; a window never crosses into another block
  * overlapping windows give each row once, in the original order with the
    original index labels. The old select_neighbors repeated rows (two
    adjacent targets, two rows each side: 10 rows for 6 distinct ones)
  * any index works, because windows are built from row positions, not from
    index labels (the old code raised IndexError on a shifted index)
  * query is 1 for target rows and 0 for the others (int64); strand='same'
    keeps only the runs of same-strand rows that hold a target
  * no target gives an empty frame with the same columns; a negative window,
    an unknown strand, a rule that cannot be evaluated and a missing block_id
    raise instead of returning something misleading
  * select_neighbors warns (FutureWarning, pointing at the caller) and
    gives the old output: query untouched, a new_query column, strand_group
    when strand=True, stats= as a value_counts_full table, unknown keywords
    as TypeError. Its before counts rows below and after rows above, as it
    always did; trim_blocks counts them the other way round
  * value_counts_full names its columns block_id_count and pid_count; it
    used to discard the rename

The data are one hand-built table of ten CDS rows in two blocks (A: p0 to p5,
strands + + - - + +; B: p6 to p9, strands + - - +); nothing touches the
network. Where windows overlap the expected rows are worked out by hand in
each test.

Run under pytest: PYTHONPATH=lib pytest test/genome/test_trim_blocks.py
"""

import numpy as np
import pandas as pd
import pytest

from rotifer.genome.data import NeighborhoodDF


def make_table(queries=(2, 3, 7), index=None):
    """Ten CDS rows in two blocks: p0..p5 in block A, p6..p9 in block B.

    Rows listed in `queries` get query == 1. Strand is + + - - + + for
    block A and + - - + for block B. `index` replaces the default 0..9
    index labels.
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


def old_select(df, **kwargs):
    """Call the deprecated select_neighbors, checking that it warns."""
    with pytest.warns(FutureWarning, match="trim_blocks"):
        return df.select_neighbors(**kwargs)


# --- trim_blocks ----------------------------------------------------------------

def test_overlapping_windows_give_each_row_once():
    """Adjacent targets 2 and 3, one row each side: p1..p4 once, not twice."""
    out = make_table().trim_blocks(before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']
    assert out['query'].tolist() == [0, 1, 1, 0, 0, 1, 0]


def test_wide_overlap_is_not_duplicated():
    """Two rows each side of 2 and 3 cover p0..p5 once (10 rows with the old code)."""
    out = make_table(queries=(2, 3)).trim_blocks(before=2, after=2)
    assert pids(out) == ['p0', 'p1', 'p2', 'p3', 'p4', 'p5']


def test_shifted_index_keeps_labels():
    """A non-RangeIndex works and its labels are kept."""
    labels = [50, 40, 30, 20, 10, 0, 90, 80, 70, 60]
    out = make_table(index=labels).trim_blocks(before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']
    assert out.index.tolist() == [40, 30, 20, 10, 90, 80, 70]


def test_duplicated_index_labels_work():
    """Repeated index labels do not mix rows up."""
    out = make_table(index=[0] * 10).trim_blocks(before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']


def test_windows_are_clipped_at_block_edges():
    """Targets on the first row of A and last row of B stay inside their block."""
    out = make_table(queries=(0, 9)).trim_blocks(before=2, after=2)
    assert pids(out) == ['p0', 'p1', 'p2', 'p7', 'p8', 'p9']


def test_two_blocks_are_independent():
    """A window never reaches into the next block: p5 and p6 are targets, each keeps one row."""
    out = make_table(queries=(5, 6)).trim_blocks(before=1, after=1)
    assert pids(out) == ['p4', 'p5', 'p6', 'p7']


def test_blocks_without_targets_are_dropped():
    """Block B has no target, so none of its rows is kept."""
    out = make_table(queries=(2,)).trim_blocks(before=5, after=5)
    assert pids(out) == ['p0', 'p1', 'p2', 'p3', 'p4', 'p5']


def test_before_is_above_and_after_is_below():
    """before=2, after=0 on row 3 keeps p1..p3: two rows above, none below."""
    out = make_table(queries=(3,)).trim_blocks(before=2, after=0)
    assert pids(out) == ['p1', 'p2', 'p3']


def test_zero_window_keeps_only_targets():
    """before=0 and after=0 keep exactly the target rows."""
    out = make_table().trim_blocks(before=0, after=0)
    assert pids(out) == ['p2', 'p3', 'p7']
    assert out['query'].tolist() == [1, 1, 1]


def test_query_is_int64_and_marks_the_targets():
    """query is 1 for the rows a rule selects, whatever the old query column said."""
    out = make_table().trim_blocks("pid == 'p5'", before=1, after=0)
    assert pids(out) == ['p4', 'p5']
    assert out['query'].tolist() == [0, 1]
    assert out['query'].dtype == np.int64


def test_rules_combine_with_and():
    """A list of rules selects the rows where every rule holds: query 1 and strand -1."""
    out = make_table(queries=(1, 2, 3)).trim_blocks(['query == 1', 'strand == -1'], before=0, after=0)
    assert pids(out) == ['p2', 'p3']


def test_series_target_is_aligned_on_the_index():
    """A boolean Series lacking some labels counts them as False."""
    table = make_table()
    target = pd.Series(True, index=[3, 4])
    out = table.trim_blocks(target, before=0, after=0)
    assert pids(out) == ['p3', 'p4']


def test_strand_same_keeps_runs_with_a_target():
    """Targets p3 and p7, one row each side: only the same-strand runs holding them stay.

    Windows are p2 p3 p4 (strands - - +) and p6 p7 p8 (+ - -). The runs
    p2 p3 and p7 p8 hold a target; p4 and p6 do not.
    """
    out = make_table(queries=(3, 7)).trim_blocks(before=1, after=1, strand='same')
    assert pids(out) == ['p2', 'p3', 'p7', 'p8']
    assert out['query'].tolist() == [0, 1, 1, 0]


def test_strand_same_with_one_strand_keeps_the_window():
    """When the whole window has the target's strand, nothing is cut."""
    out = make_table(queries=(0,)).trim_blocks(before=0, after=1, strand='same')
    assert pids(out) == ['p0', 'p1']


def test_no_target_gives_an_empty_frame_with_the_same_columns():
    """No matching row: empty, same columns, query int64."""
    table = make_table()
    out = table.trim_blocks("pid == 'nothing'")
    assert len(out) == 0
    assert list(out.columns) == list(table.columns)
    assert out['query'].dtype == np.int64
    assert isinstance(out, NeighborhoodDF)


def test_input_is_not_modified():
    """The table keeps its rows and its query column."""
    table = make_table()
    before = table.copy()
    table.trim_blocks("pid == 'p5'", before=1, after=1)
    pd.testing.assert_frame_equal(pd.DataFrame(table), pd.DataFrame(before))


@pytest.mark.parametrize("kwargs", [dict(before=-1), dict(after=-1), dict(before=-2, after=-2)])
def test_negative_window_raises(kwargs):
    """A negative number of rows is an error."""
    with pytest.raises(ValueError, match="must not be negative"):
        make_table().trim_blocks(**kwargs)


@pytest.mark.parametrize("strand", ["+", "-", None, True, "different"])
def test_bad_strand_raises(strand):
    """strand is 'any' or 'same'; nothing else."""
    with pytest.raises(ValueError, match=r"strand=.* not supported.*'any'.*'same'"):
        make_table().trim_blocks(strand=strand)


def test_bad_rule_raises_value_error_naming_the_rule():
    """A rule that cannot be evaluated raises ValueError, chained from the cause."""
    with pytest.raises(ValueError, match="nonexistent == 1") as err:
        make_table().trim_blocks("nonexistent == 1")
    assert err.value.__cause__ is not None


def test_target_of_wrong_type_raises():
    """A target that is neither a string nor a Series raises TypeError."""
    with pytest.raises(TypeError, match="targets must be"):
        make_table().trim_blocks(3)


def test_missing_block_id_raises_key_error():
    """Without a block_id column there are no blocks to cut."""
    with pytest.raises(KeyError, match="block_id"):
        make_table().drop(columns='block_id').trim_blocks()


# --- select_neighbors (deprecated shim) -----------------------------------------

def test_shim_warns_where_it_is_called():
    """The FutureWarning names trim_blocks, the removal date, and points at the caller."""
    with pytest.warns(FutureWarning) as record:
        make_table().select_neighbors()
    assert "trim_blocks" in str(record[0].message)
    assert "2026-11-30" in str(record[0].message)
    assert record[0].filename == __file__


def test_shim_overlapping_windows_give_each_row_once():
    """Same rows as trim_blocks for symmetric windows: p1..p4 and p6..p8, once each."""
    out = old_select(make_table(), before=1, after=1)
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']
    assert out['new_query'].tolist() == [0, 1, 1, 0, 0, 1, 0]


def test_shim_shifted_index_keeps_labels():
    """A non-RangeIndex works (IndexError before) and its labels are kept."""
    labels = [50, 40, 30, 20, 10, 0, 90, 80, 70, 60]
    out = old_select(make_table(index=labels), before=1, after=1)
    assert out.index.tolist() == [40, 30, 20, 10, 90, 80, 70]


def test_shim_windows_are_clipped_at_block_edges():
    """Targets at the two ends of the table stay inside their block."""
    out = old_select(make_table(queries=(0, 9)), before=2, after=2)
    assert pids(out) == ['p0', 'p1', 'p2', 'p7', 'p8', 'p9']


def test_shim_two_blocks_are_independent():
    """p5 and p6 are targets in different blocks and keep one row each side."""
    out = old_select(make_table(queries=(5, 6)), before=1, after=1)
    assert pids(out) == ['p4', 'p5', 'p6', 'p7']


def test_shim_before_counts_rows_below():
    """Old convention: before=2, after=0 on row 3 keeps p3..p5."""
    out = old_select(make_table(queries=(3,)), before=2, after=0)
    assert pids(out) == ['p3', 'p4', 'p5']


def test_shim_leaves_query_and_adds_new_query():
    """query keeps the table's values; new_query marks the rows that match column == value."""
    out = old_select(make_table(queries=(3,)), before=1, after=1)
    assert out['query'].tolist() == [0, 1, 0]
    assert out['new_query'].tolist() == [0, 1, 0]
    assert list(out.columns) == list(make_table().columns) + ['new_query']


def test_shim_column_and_value_pick_the_targets():
    """column='pid', value='p8' targets that row; query is untouched."""
    out = old_select(make_table(), column='pid', value='p8', before=1, after=1)
    assert pids(out) == ['p7', 'p8', 'p9']
    assert out['new_query'].tolist() == [0, 1, 0]
    assert out['query'].tolist() == [1, 0, 0]


def test_shim_new_query_false_drops_the_column():
    """new_query=False returns the rows without the new_query column."""
    out = old_select(make_table(), before=1, after=1, new_query=False)
    assert 'new_query' not in out.columns
    assert pids(out) == ['p1', 'p2', 'p3', 'p4', 'p6', 'p7', 'p8']


def test_shim_strand_adds_strand_group():
    """strand=True keeps same-strand runs with a target and keeps their strand_group numbers.

    The selected rows p1 p2 p3 p4 p6 p7 p8 form runs 1, 2 2, 3, 4, 5 5; only
    runs 2 and 5 hold a target.
    """
    out = old_select(make_table(), before=1, after=1, strand=True)
    assert pids(out) == ['p2', 'p3', 'p7', 'p8']
    assert out['strand_group'].tolist() == [2, 2, 5, 5]
    assert out['new_query'].tolist() == [1, 1, 1, 0]


def test_shim_stats_gives_value_counts_full():
    """stats='strand' returns the frequency table of the (strand-filtered) rows.

    After strand=True the rows p2 p3 p7 p8 all have strand -1, in 2 blocks.
    """
    out = old_select(make_table(), before=1, after=1, strand=True, stats='strand')
    expected = pd.DataFrame({
        'strand': [-1], 'block_id_count': [2], 'pid_count': [4],
        'relative_block_id': [100.0], 'relative_pid': [100.0],
    })
    pd.testing.assert_frame_equal(out, expected, check_dtype=False)


def test_shim_unknown_keyword_raises_type_error():
    """Keywords other than before, after, strand and stats are an error."""
    with pytest.raises(TypeError, match="sideways"):
        make_table().select_neighbors(sideways=1)


def test_shim_no_match_gives_an_empty_frame():
    """No row with column == value: an empty frame with the new_query column."""
    out = old_select(make_table(), column='pid', value='nothing')
    assert len(out) == 0
    assert 'new_query' in out.columns


# --- value_counts_full ------------------------------------------------------------

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
