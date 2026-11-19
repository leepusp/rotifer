"""
Tests for rotifer.interval.utils.complement.

What has to hold, and why:

- Gaps are listed before the first region (from position 1), between
  regions, and after the last one up to the sequence length; each gap is
  inclusive and 1-based, like the regions.
- Overlapping and nested regions are merged first, so a gap never has a
  negative length and no position is in two gaps (the old code returned a
  gap 91..19 of length -71 for regions 10-90 and 20-40 on 100 residues).
- A region starting at 1 has no gap before it, and one reaching the end of
  the sequence has no gap after it.
- Without a length column there is no trailing gap and no length column.
- ``reference``, ``start``, ``end`` and ``length`` may be named anything,
  and ``reference`` may be several columns.
- ``include``, ``exclude`` and ``add`` shape the columns of the result;
  ``add`` with a dict used to raise NameError, and ``complete.py`` calls
  ``complement(df, include, exclude, values)`` positionally.
- Empty input gives an empty table, the input is not modified, ``start``
  and ``end`` are int64 and the index is 0..n-1.

The fixtures are hand-built tables small enough to check by eye: a
sequence of length 100 with the regions of each test, and a second
sequence that is fully covered and has no gap.

Run under pytest: PYTHONPATH=lib pytest test/interval/test_complement.py
"""

import pandas as pd
import pytest

from rotifer.interval.utils import complement

def _table(regions, length=100, ids=None):
    """
    Build a table with one row per (start, end) pair of ``regions``, all on
    sequence s1 (or on ``ids``, one per row) of the given length.
    """
    ids = ids or ['s1'] * len(regions)
    return pd.DataFrame({'ID': ids, 'start': [r[0] for r in regions],
                         'end': [r[1] for r in regions], 'length': length})

def _gaps(out):
    """The (start, end) pairs of a result, as a list of tuples."""
    return list(zip(out['start'].tolist(), out['end'].tolist()))

def test_before_between_after():
    """Gaps before, between and after the regions, with their lengths."""
    out = complement(_table([(10, 30), (50, 80)]))
    assert out.to_dict('list') == {'ID': ['s1'] * 3, 'start': [1, 31, 81],
                                   'end': [9, 49, 100], 'length': [9, 19, 20]}

def test_fully_covered_sequence_has_no_gap():
    """A sequence covered from 1 to its length has no row."""
    out = complement(_table([(1, 40)], length=40, ids=['s2']))
    assert out.empty

def test_adjacent_regions_leave_no_gap():
    """Regions that touch (30 then 31) leave no gap between them."""
    assert _gaps(complement(_table([(1, 30), (31, 60)], length=60))) == []

def test_region_starting_at_one():
    """No gap before a region that starts at 1."""
    assert _gaps(complement(_table([(1, 20)]))) == [(21, 100)]

def test_region_reaching_the_end():
    """No gap after a region that ends at the sequence length."""
    assert _gaps(complement(_table([(30, 100)]))) == [(1, 29)]

def test_nested_region():
    """
    A region inside another changes nothing. The old code returned a gap
    91..19 of length -71 for 10-90 and 20-40.
    """
    out = complement(_table([(10, 90), (20, 40)]))
    assert _gaps(out) == [(1, 9), (91, 100)]
    assert (out['length'] > 0).all()

def test_nested_region_listed_first_and_later():
    """The merge does not depend on row order."""
    a = complement(_table([(20, 40), (10, 90)]))
    b = complement(_table([(10, 90), (20, 40)]))
    pd.testing.assert_frame_equal(a, b)

def test_overlapping_regions():
    """Overlapping regions are merged: 10-50 and 40-70 cover 10-70."""
    assert _gaps(complement(_table([(10, 50), (40, 70)]))) == [(1, 9), (71, 100)]

def test_region_after_nested_one_keeps_running_end():
    """
    The end that counts is the largest so far, not the previous row's:
    after 10-90 and 20-40, a region at 50-60 leaves no gap before it.
    """
    assert _gaps(complement(_table([(10, 90), (20, 40), (50, 60)]))) == [(1, 9), (91, 100)]

def test_no_length_column():
    """Without length there is no trailing gap and no length column."""
    regions = _table([(10, 30), (50, 80)]).drop(columns='length')
    out = complement(regions)
    assert out.columns.tolist() == ['ID', 'start', 'end']
    assert _gaps(out) == [(1, 9), (31, 49)]

def test_custom_names():
    """Reference, start, end and length columns can be named anything."""
    regions = pd.DataFrame({'seq': ['a', 'a'], 'b': [10, 50], 'e': [30, 80], 'n': [100, 100]})
    out = complement(regions, reference='seq', start='b', end='e', length='n')
    assert out.to_dict('list') == {'seq': ['a'] * 3, 'b': [1, 31, 81], 'e': [9, 49, 100], 'n': [9, 19, 20]}

def test_several_reference_columns():
    """Sequences are told apart by all reference columns, and sorted by them."""
    regions = pd.DataFrame({'db': ['x', 'x', 'y'], 'ID': ['s', 's', 's'],
                            'start': [10, 50, 10], 'end': [30, 80, 100], 'length': [100, 100, 100]})
    out = complement(regions, reference=['db', 'ID'])
    assert out['db'].tolist() == ['x', 'x', 'x', 'y']
    assert _gaps(out) == [(1, 9), (31, 49), (81, 100), (1, 9)]

def test_sorted_by_reference_and_start():
    """Result is sorted by reference then start, whatever the input order."""
    regions = _table([(50, 80), (10, 30), (5, 20)], ids=['b', 'a', 'b'])
    out = complement(regions)
    assert out['ID'].tolist() == ['a', 'a', 'b', 'b', 'b']
    assert out['start'].tolist() == [1, 31, 1, 21, 81]

def test_include_true_keeps_input_columns():
    """include=True gives the columns of df in its order; others are missing values."""
    regions = _table([(10, 30)]).assign(domain='A').loc[:, ['domain', 'length', 'ID', 'end', 'start']]
    out = complement(regions)
    assert out.columns.tolist() == ['domain', 'length', 'ID', 'end', 'start']
    assert out['domain'].isna().all()

def test_include_false_gives_computed_columns():
    """include=False gives reference, start, end and length only."""
    regions = _table([(10, 30)]).assign(domain='A')
    out = complement(regions, include=False)
    assert out.columns.tolist() == ['ID', 'start', 'end', 'length']

def test_include_list():
    """include as a list gives these columns in this order."""
    regions = _table([(10, 30)]).assign(domain='A')
    out = complement(regions, include=['end', 'ID', 'domain'])
    assert out.columns.tolist() == ['end', 'ID', 'domain']
    assert out['end'].tolist() == [9, 100]

def test_exclude():
    """exclude drops columns; names not in the result are ignored."""
    out = complement(_table([(10, 30)]), exclude=['length', 'nothing'])
    assert out.columns.tolist() == ['ID', 'start', 'end']

def test_add_dict():
    """add sets a column on every gap (the old code raised NameError)."""
    out = complement(_table([(10, 30)]).assign(domain='A'), add={'domain': 'unk'})
    assert out['domain'].tolist() == ['unk', 'unk']

def test_positional_call_as_complete_does():
    """complete.py calls complement(df, include, exclude, values) positionally."""
    regions = _table([(10, 30)]).assign(domain='A')
    out = complement(regions, True, ['length'], {'domain': 'unk'})
    assert out.columns.tolist() == ['ID', 'start', 'end', 'domain']
    assert out['domain'].tolist() == ['unk', 'unk']

def test_empty_input():
    """An empty table gives an empty table with the same columns."""
    regions = _table([(10, 30)]).iloc[:0]
    out = complement(regions)
    assert out.empty
    assert out.columns.tolist() == regions.columns.tolist()

def test_input_not_modified():
    """The input table keeps its rows, columns and order."""
    regions = _table([(50, 80), (10, 30)])
    expected = regions.copy()
    complement(regions)
    pd.testing.assert_frame_equal(regions, expected)

def test_dtypes_and_index():
    """start and end are int64 (even from float input) and the index is 0..n-1."""
    regions = _table([(50, 80), (10, 30)], ids=['b', 'a']).astype({'start': float, 'end': float})
    out = complement(regions)
    assert out['start'].dtype == 'int64' and out['end'].dtype == 'int64'
    assert out.index.tolist() == list(range(len(out)))

def test_invalid_rows_cover_nothing():
    """A row with a missing end, or with end < start, covers nothing."""
    regions = _table([(10, 30), (50, None), (70, 60)])
    assert _gaps(complement(regions)) == [(1, 9), (31, 100)]
