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

import pandas as pd
import pytest

import rotifer.interval.utils as riu

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
    The result has a fresh 0..n-1 index from an internal sort, not the
    input index: Y and Z are input rows 1 and 2.
    """
    out = riu.filter_nonoverlapping_regions(_hhsuite())
    assert out['hit'].tolist() == ['Y', 'Z']
    assert out.index.tolist() == [0, 2]

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
