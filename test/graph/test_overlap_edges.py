#!/usr/bin/env python3
"""Tests for rotifer.graph.network.overlap_edges.

What has to hold, and why:

  * the edges equal a brute-force set computation, for all three metrics,
    so the sparse product and the weight formulas are right
  * duplicate (group, member) rows and missing values in either column
    are ignored, as the docstring promises
  * min_members drops small groups before pairs are made, and min_weight
    is inclusive: a pair exactly at the threshold stays
  * each pair appears once, source sorts before target, no group is paired
    with itself, and the rows are sorted by source and target
  * _overlap_pairs gives the same pairs whatever the block size
  * integer labels stay integers and string labels work
  * empty input, or groups that share nothing, give an empty frame with
    the four columns
  * bad arguments raise KeyError or ValueError, and the input frame is
    not modified
  * a table of 20,000 groups finishes and gives a valid frame

The data are small seeded random tables; nothing touches the network.

Run under pytest: PYTHONPATH=lib pytest test/graph/test_overlap_edges.py
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip('networkx')
pytest.importorskip('scipy')

import scipy.sparse

from rotifer.graph.network import _overlap_pairs, overlap_edges

COLUMNS = ['source', 'target', 'shared', 'weight']


def random_table(seed, n_groups=40, n_members=25, n_rows=200):
    """
    Build a random long table of (group, member) assignments.

    Parameters
    ----------
    seed : int
        Seed of the random generator.
    n_groups, n_members : int
        Number of distinct group and member labels to draw from.
    n_rows : int
        Number of rows, so duplicates are likely.

    Returns
    -------
    pandas.DataFrame
        Columns ``g`` (int labels) and ``m`` (str labels).
    """
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        'g': rng.integers(0, n_groups, n_rows),
        'm': [f'm{i}' for i in rng.integers(0, n_members, n_rows)]})


def brute_force(df, metric, min_members=1, min_weight=0):
    """
    Compute the expected edges with plain Python sets.

    Parameters
    ----------
    df : pandas.DataFrame
        Table with columns ``g`` and ``m``, without missing values.
    metric : {'count', 'jaccard', 'overlap'}
        Score of a pair.
    min_members : int, default 1
        Groups with fewer distinct members are left out.
    min_weight : float, default 0
        Pairs scoring less are left out.

    Returns
    -------
    pandas.DataFrame
        Edges with the columns of ``overlap_edges``.
    """
    sets = {}
    for g, m in zip(df['g'], df['m']):
        sets.setdefault(g, set()).add(m)
    sets = {g: s for g, s in sets.items() if len(s) >= min_members}
    out = []
    for a in sorted(sets):
        for b in sorted(sets):
            if a >= b:
                continue
            shared = len(sets[a] & sets[b])
            if not shared:
                continue
            if metric == 'count':
                w = shared
            elif metric == 'jaccard':
                w = shared / len(sets[a] | sets[b])
            else:
                w = shared / min(len(sets[a]), len(sets[b]))
            if w >= min_weight:
                out.append((a, b, shared, w))
    expected = pd.DataFrame(out, columns=COLUMNS)
    expected['shared'] = expected['shared'].astype('int64')
    expected['weight'] = expected['weight'].astype(
        'int64' if metric == 'count' else 'float64')
    return expected


def assert_same(got, expected):
    """
    Compare two edge frames, with a tolerance on float weights.

    Parameters
    ----------
    got, expected : pandas.DataFrame
        Frames to compare.

    Returns
    -------
    None
    """
    got = got.copy()
    expected = expected.copy()
    for frame in (got, expected):
        frame['source'] = frame['source'].astype('int64')
        frame['target'] = frame['target'].astype('int64')
    pd.testing.assert_frame_equal(got, expected, rtol=1e-12)


@pytest.mark.parametrize('metric', ['count', 'jaccard', 'overlap'])
@pytest.mark.parametrize('seed', [0, 1, 2, 3])
def test_matches_brute_force(seed, metric):
    """Edges equal the set-based computation for every metric."""
    df = random_table(seed)
    assert_same(overlap_edges(df, 'g', 'm', metric=metric),
                brute_force(df, metric))


@pytest.mark.parametrize('metric', ['count', 'jaccard', 'overlap'])
def test_filters_match_brute_force(metric):
    """min_members and min_weight together match the set computation."""
    df = random_table(7)
    got = overlap_edges(df, 'g', 'm', metric=metric,
                        min_members=4, min_weight=2 if metric == 'count' else 0.2)
    expected = brute_force(df, metric, min_members=4,
                           min_weight=2 if metric == 'count' else 0.2)
    assert len(expected) > 0
    assert_same(got, expected)


def test_duplicates_and_missing_ignored():
    """Repeated rows and NaN in either column change nothing."""
    df = random_table(5)
    noisy = pd.concat([df, df.iloc[:30],
                       pd.DataFrame({'g': [1.0, np.nan], 'm': [np.nan, 'm1']})],
                      ignore_index=True)
    for metric in ('count', 'jaccard'):
        a = overlap_edges(df, 'g', 'm', metric=metric)
        b = overlap_edges(noisy, 'g', 'm', metric=metric)
        pd.testing.assert_frame_equal(a.astype({'source': 'float64', 'target': 'float64'}),
                                      b, check_dtype=True)


def test_min_members_drops_small_groups():
    """A group with too few distinct members is not paired."""
    df = pd.DataFrame({'g': [1, 1, 1, 2, 2, 3], 'm': list('abcabc')[:5] + ['a']})
    got = overlap_edges(df, 'g', 'm', metric='count', min_members=2)
    assert got[['source', 'target']].values.tolist() == [[1, 2]]
    assert 3 not in set(got['source']) | set(got['target'])


def test_min_weight_is_inclusive():
    """A pair whose score equals min_weight is kept."""
    df = pd.DataFrame({'g': [1, 1, 2, 2], 'm': list('abbc')})
    assert len(overlap_edges(df, 'g', 'm', min_weight=1 / 3)) == 1
    assert len(overlap_edges(df, 'g', 'm', min_weight=0.34)) == 0
    assert len(overlap_edges(df, 'g', 'm', metric='count', min_weight=1)) == 1
    assert len(overlap_edges(df, 'g', 'm', metric='count', min_weight=2)) == 0


def test_pair_structure_and_order():
    """Pairs are unique, ordered source < target, and sorted."""
    got = overlap_edges(random_table(11), 'g', 'm')
    assert (got['source'] < got['target']).all()
    assert not got.duplicated(['source', 'target']).any()
    keys = list(zip(got['source'], got['target']))
    assert keys == sorted(keys)
    assert isinstance(got.index, pd.RangeIndex)
    assert list(got.columns) == COLUMNS


@pytest.mark.parametrize('block', [1, 2, 1000])
@pytest.mark.parametrize('metric', ['count', 'jaccard', 'overlap'])
def test_block_size_does_not_matter(metric, block):
    """_overlap_pairs returns the same pairs for every block size."""
    df = random_table(13, n_groups=30)
    df = df.drop_duplicates()
    g, _ = pd.factorize(df['g'], sort=True)
    m, _ = pd.factorize(df['m'], sort=True)
    inc = scipy.sparse.csr_matrix((np.ones(len(g), dtype=np.int64), (g, m)))
    sizes = np.asarray(inc.sum(axis=1)).ravel().astype(np.int64)
    ref = _overlap_pairs(inc, sizes, metric, 0.1, block=1000)
    got = _overlap_pairs(inc, sizes, metric, 0.1, block=block)
    for a, b in zip(ref, got):
        assert a.dtype == b.dtype
        np.testing.assert_allclose(a, b)
    assert len(ref[0]) > 0


def test_integer_and_string_labels():
    """Integer group labels stay integers; strings work."""
    df = pd.DataFrame({'g': [1, 1, 2, 2], 'm': list('abbc')})
    assert pd.api.types.is_integer_dtype(overlap_edges(df, 'g', 'm')['source'])
    s = df.assign(g=df['g'].map({1: 'x', 2: 'y'}))
    got = overlap_edges(s, 'g', 'm')
    assert got[['source', 'target']].values.tolist() == [['x', 'y']]


def test_empty_input_and_no_overlap():
    """Empty tables and disjoint groups give an empty four-column frame."""
    empty = pd.DataFrame({'g': pd.Series([], dtype=object),
                          'm': pd.Series([], dtype=object)})
    disjoint = pd.DataFrame({'g': [1, 2], 'm': ['a', 'b']})
    for df in (empty, disjoint):
        for metric in ('count', 'jaccard'):
            got = overlap_edges(df, 'g', 'm', metric=metric)
            assert list(got.columns) == COLUMNS
            assert len(got) == 0
            assert got['shared'].dtype == np.int64
            assert got['weight'].dtype == (
                np.int64 if metric == 'count' else np.float64)


def test_bad_arguments():
    """Missing columns raise KeyError; bad values raise ValueError."""
    df = pd.DataFrame({'g': [1], 'm': ['a']})
    with pytest.raises(KeyError, match='nope'):
        overlap_edges(df, 'g', 'nope')
    with pytest.raises(ValueError):
        overlap_edges(df, 'g', 'm', metric='cosine')
    with pytest.raises(ValueError):
        overlap_edges(df, 'g', 'm', min_members=0)
    with pytest.raises(ValueError):
        overlap_edges(df, 'g', 'm', min_weight=-1)


def test_input_not_modified():
    """The input frame is the same after the call."""
    df = random_table(3)
    df.loc[0, 'm'] = np.nan
    before = df.copy()
    overlap_edges(df, 'g', 'm')
    pd.testing.assert_frame_equal(df, before)


def test_scale_smoke():
    """20,000 groups of 10 members out of 50,000 finish with a valid frame."""
    rng = np.random.default_rng(0)
    n = 20000
    members = np.array([rng.choice(50000, 10, replace=False) for _ in range(n)])
    df = pd.DataFrame({'g': np.repeat(np.arange(n), 10),
                       'm': members.ravel()})
    df.loc[df['g'] % 100 == 1, 'm'] = df['m'] % 12  # plant overlapping groups
    got = overlap_edges(df, 'g', 'm', min_weight=0.3)
    assert len(got) > 0
    assert list(got.columns) == COLUMNS
    assert (got['source'] < got['target']).all()
    assert (got['weight'] >= 0.3).all()
