"""
Tests for rotifer.interval.utils: filter_nonoverlapping_regions,
architecture and complement.

What has to hold, and why:

- With the ``hmmer`` preset, the best region of each sequence (lowest
  E-value, then longest) is kept, and every other region is kept only if
  it does not overlap a kept region by more than the allowed overlap.
  Devel code (``devel.beta.hmmer.add_arch_to_df``) builds Pfam
  architectures this way.
- A fractional ``maximum_overlap`` (below 1) is a fraction of the shorter
  region; a value of 1 or more is a number of residues.
- ``architecture`` joins the kept regions of each sequence with ``+`` in
  coordinate order.
- The default columns are the HH-suite ones (query, hit, qstart, qend).
- ``complement`` returns the gaps between regions, including the gap
  after the last region when a ``length`` column is given.

The fixtures are hand-built tables small enough to check by eye: two
proteins with five HMMER-like hits, one HH-suite query with three hits,
and two sequences with regions and a length column.

Run under pytest: PYTHONPATH=lib pytest test/interval/test_utils.py
"""

import warnings

import pandas as pd
import pytest

import rotifer.interval.utils as riu

# The old behaviour is pinned below; only the tests that check the
# deprecation look at the warning.
pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")

def _hmmer():
    """
    Five hmmscan-like hits on two proteins, columns as hmmer2table writes
    them (sequence, model, estart, eend, evalue).

    On p1, A (1-60) is best; B (50-100) overlaps A by 11 residues, below
    40% of B's 51; C (90-200) overlaps B by 11, below 40% of B; D (10-30)
    lies inside A. p2 has a single hit, E.
    """
    return pd.DataFrame({
        'sequence': ['p1', 'p1', 'p1', 'p1', 'p2'],
        'model':    ['A', 'B', 'C', 'D', 'E'],
        'estart':   [1, 50, 90, 10, 5],
        'eend':     [60, 100, 200, 30, 80],
        'evalue':   [1e-10, 1e-5, 1e-3, 1e-2, 1e-4],
    })

def _hhsuite():
    """
    One HH-suite query with three hits. Y has the largest sum_probs and is
    kept; X overlaps Y by 71 residues and is dropped; Z (120-200) overlaps
    Y by 31 residues, below 40% of Z's 81.
    """
    return pd.DataFrame({
        'query':       ['q', 'q', 'q'],
        'hit':         ['X', 'Y', 'Z'],
        'qstart':      [1, 30, 120],
        'qend':        [100, 150, 200],
        'probability': [99.0, 95.0, 90.0],
        'sum_probs':   [80.0, 90.0, 10.0],
        'score':       [50.0, 60.0, 5.0],
    })

def test_filter_hmmer_preset():
    """The hmmer preset keeps A, B and C on p1 and E on p2; D lies inside A."""
    out = riu.filter_nonoverlapping_regions(_hmmer(), **riu.config['hmmer'])
    assert out['model'].tolist() == ['A', 'B', 'C', 'E']
    assert out.columns.tolist() == _hmmer().columns.tolist()

def test_filter_residue_overlap():
    """With maximum_overlap=5 residues, B's 11-residue overlap with A is too much."""
    out = riu.filter_nonoverlapping_regions(_hmmer(), **{**riu.config['hmmer'], 'maximum_overlap': 5})
    assert out['model'].tolist() == ['A', 'C', 'E']

def test_filter_does_not_modify_input():
    """The input table keeps its columns and rows."""
    df = _hmmer()
    riu.filter_nonoverlapping_regions(df, **riu.config['hmmer'])
    pd.testing.assert_frame_equal(df, _hmmer())

def test_filter_index():
    """
    The result keeps the input index labels: Y and Z are input rows 1 and 2.

    Changed: the old code returned a fresh index from an internal sort
    (labels 0 and 2, which matched no input row).
    """
    out = riu.filter_nonoverlapping_regions(_hhsuite())
    assert out['hit'].tolist() == ['Y', 'Z']
    assert out.index.tolist() == [1, 2]

def test_architecture_hmmer_preset():
    """Kept regions are joined with '+' in coordinate order, one row per sequence."""
    out = riu.architecture(_hmmer(), **riu.config['hmmer'])
    assert out.to_dict('list') == {'sequence': ['p1', 'p2'], 'architecture': ['A+B+C', 'E']}

def test_architecture_default_columns():
    """The defaults read HH-suite columns and rank by sum_probs first."""
    out = riu.architecture(_hhsuite())
    assert out.to_dict('list') == {'query': ['q'], 'architecture': ['Y+Z']}

def test_architecture_reference_list():
    """Several reference columns group together and all appear in the result."""
    out = riu.architecture(_hmmer().assign(g='x'), **{**riu.config['hmmer'], 'reference': ['g', 'sequence']})
    assert out.to_dict('list') == {'g': ['x', 'x'], 'sequence': ['p1', 'p2'], 'architecture': ['A+B+C', 'E']}

def test_architecture_empty():
    """An empty table gives an empty architecture table with the same key columns."""
    out = riu.architecture(_hmmer().iloc[:0], **riu.config['hmmer'])
    assert out.empty
    assert out.columns.tolist() == ['sequence', 'architecture']

def test_complement_with_length():
    """
    Gaps before, between and after the regions of s1 (length 100); s2 is
    covered from 1 to its length, so it has no gap.
    """
    regions = pd.DataFrame({'ID': ['s1', 's1', 's2'], 'start': [10, 50, 1],
                            'end': [30, 80, 40], 'length': [100, 100, 40]})
    out = riu.complement(regions)
    assert out.to_dict('list') == {'ID': ['s1', 's1', 's1'], 'start': [1, 31, 81],
                                   'end': [9, 49, 100], 'length': [9, 19, 20]}

def _caught(func, *args, **kwargs):
    """Call ``func`` and return its FutureWarnings (it is called from this file)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        func(*args, **kwargs)
    return [w for w in caught if issubclass(w.category, FutureWarning)]

def test_filter_warns_at_caller():
    """One FutureWarning, naming the replacement and pointing at this file."""
    caught = _caught(riu.filter_nonoverlapping_regions, _hmmer(), **riu.config['hmmer'])
    assert len(caught) == 1
    assert 'select_regions' in str(caught[0].message)
    assert caught[0].filename == __file__

def test_architecture_warns_once_at_caller():
    """architecture warns once (not also for the filter) and points at this file."""
    caught = _caught(riu.architecture, _hmmer(), **riu.config['hmmer'])
    assert len(caught) == 1
    assert 'domain2architecture' in str(caught[0].message)
    assert caught[0].filename == __file__

def test_fractional_overlap_without_region_length():
    """
    A fractional overlap works when criteria has no 'region_length'.

    Changed: the old code raised AttributeError, because it read the
    region_length column that it only created when ranking by it.
    """
    out = riu.filter_nonoverlapping_regions(_hmmer(), **{**riu.config['hmmer'], 'criteria': {'evalue': True}})
    assert out['model'].tolist() == ['A', 'B', 'C', 'E']

@pytest.mark.parametrize('order', [[0, 1], [1, 0]])
def test_tie_broken_by_smallest_start(order):
    """
    Two equally good overlapping regions: the one with the smaller start
    wins, whatever the row order.

    Changed: the old code kept the first row in input order.
    """
    df = pd.DataFrame({'sequence': ['p', 'p'], 'model': ['late', 'early'],
                       'estart': [20, 1], 'eend': [60, 41], 'evalue': [1e-5, 1e-5]}).iloc[order]
    out = riu.filter_nonoverlapping_regions(df, **{**riu.config['hmmer'], 'maximum_overlap': 0})
    assert out['model'].tolist() == ['early']

def test_overlap_one_and_fraction_are_different_units():
    """
    maximum_overlap=1 is one residue (A and B share 11), not 100%.
    """
    out = riu.filter_nonoverlapping_regions(_hmmer(), **{**riu.config['hmmer'], 'maximum_overlap': 1})
    assert out['model'].tolist() == ['A', 'C', 'E']

def test_overlap_non_integer_is_rounded_down():
    """5.7 residues means 5: an overlap of 11 is too much, one of 5 is not."""
    df = pd.DataFrame({'sequence': ['p', 'p'], 'model': ['A', 'B'],
                       'estart': [1, 6], 'eend': [10, 40], 'evalue': [1e-9, 1e-5]})
    assert riu.filter_nonoverlapping_regions(df, **{**riu.config['hmmer'], 'maximum_overlap': 5.7})['model'].tolist() == ['A', 'B']
    assert riu.filter_nonoverlapping_regions(df, **{**riu.config['hmmer'], 'maximum_overlap': 4.9})['model'].tolist() == ['A']

def test_default_criteria_not_mutated():
    """Calling with the default criteria leaves config['default']['criteria'] as it was."""
    before = dict(riu.config['default']['criteria'])
    riu.filter_nonoverlapping_regions(_hhsuite())
    riu.architecture(_hhsuite())
    assert dict(riu.config['default']['criteria']) == before
    assert list(riu.config['default']['criteria']) == list(before)

def test_preset_with_collapse_key_unpacks():
    """A whole preset, including its 'collapse' key, unpacks into the filter."""
    for name in ('default', 'archtsv', 'hhsuite', 'hmmer'):
        preset = riu.config[name]
        df = pd.DataFrame({k: [1] for k in preset['reference'] + [preset['start'], preset['end'], preset['collapse']]
                                         + list(preset['criteria']) if k != 'region_length'})
        out = riu.filter_nonoverlapping_regions(df, **preset)
        assert len(out) == 1
