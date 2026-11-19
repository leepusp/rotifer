#!/usr/bin/env python3
"""Tests for rotifer.genome.data.NeighborhoodDF.jaccard.

jaccard() describes each block by the set of distinct values of a column
(c80e3 by default), joins blocks whose sets have a Jaccard index of at least
min_jaccard, and numbers the Louvain communities and connected components of
that graph, 0 being the largest. Blocks with fewer than min_members values,
or with no partner, are in no community and no component (<NA>).

Fixture: seven blocks.
  g1, g2, g3 : {A, B, C}, {A, B, C}, {A, B, D}  (one tight group)
  h1, h2     : {X, Y, Z}, {X, Y, W}              (a smaller group)
  s1         : {A}                               (below min_members)
  lone       : {Q, R, S}                         (no partner)

Run under pytest: PYTHONPATH=lib pytest test/genome/test_neighborhood_jaccard.py
"""

import pandas as pd
import pytest

from rotifer.genome.data import NeighborhoodDF

pytest.importorskip('networkx')


def make(blocks, column='c80e3'):
    """
    Build a NeighborhoodDF from (block, values) pairs.

    Parameters
    ----------
    blocks : list of tuple
        ``(block_id, values)`` pairs, kept in the given order.
    column : str, default 'c80e3'
        Name of the column that holds the values.

    Returns
    -------
    NeighborhoodDF
        One row per value.
    """
    rows = [(b, v) for b, vals in blocks for v in vals]
    return NeighborhoodDF(pd.DataFrame(rows, columns=['block_id', column]))


def fixture():
    """
    Build the seven-block fixture described in the module docstring.

    Returns
    -------
    NeighborhoodDF
        The fixture table.
    """
    return make([
        ('g1', 'ABC'), ('g2', 'ABC'), ('g3', 'ABD'),
        ('h1', 'XYZ'), ('h2', 'XYW'), ('s1', 'A'), ('lone', 'QRS'),
    ])


def by_block(res):
    """
    Index a jaccard() result by block.

    Parameters
    ----------
    res : pandas.DataFrame
        Result of ``jaccard()``.

    Returns
    -------
    pandas.DataFrame
        The same rows indexed by ``block_id``.
    """
    return res.set_index('block_id')


def test_identical_sets_share_community_and_component():
    """Blocks with the same member set land together."""
    r = by_block(fixture().jaccard())
    assert r.loc['g1', 'community'] == r.loc['g2', 'community']
    assert r.loc['g1', 'component'] == r.loc['g2', 'component']
    assert r.loc['g1', 'component'] == r.loc['g3', 'component']


def test_small_and_unpartnered_blocks_are_na():
    """Blocks below min_members or without a partner get <NA>."""
    r = by_block(fixture().jaccard())
    for b in ('s1', 'lone'):
        assert pd.isna(r.loc[b, 'community'])
        assert pd.isna(r.loc[b, 'component'])


def test_threshold_is_inclusive():
    """A pair with Jaccard exactly equal to min_jaccard is joined."""
    ndf = make([('a', 'ABC'), ('b', 'BCD')])  # 2 shared / 4 = 0.5
    r = ndf.jaccard(min_jaccard=0.5)
    assert r.component.notna().all()
    r = ndf.jaccard(min_jaccard=0.51)
    assert r.component.isna().all()


def test_same_seed_same_result():
    """Two runs with the same seed agree."""
    pd.testing.assert_frame_equal(fixture().jaccard(seed=3), fixture().jaccard(seed=3))


def test_groups_numbered_by_size():
    """Separate groups get different labels and the bigger one is 0."""
    r = by_block(fixture().jaccard())
    assert r.loc[['g1', 'g2', 'g3'], 'community'].tolist() == [0, 0, 0]
    assert r.loc[['g1', 'g2', 'g3'], 'component'].tolist() == [0, 0, 0]
    assert r.loc[['h1', 'h2'], 'community'].tolist() == [1, 1]
    assert r.loc[['h1', 'h2'], 'component'].tolist() == [1, 1]


def test_members_sorted_tuple_without_missing():
    """members is a sorted tuple, missing values dropped, all-missing gives ()."""
    ndf = NeighborhoodDF(pd.DataFrame({
        'block_id': ['a', 'a', 'a', 'a', 'z', 'z'],
        'c80e3': ['C', None, 'A', 'A', None, None]}))
    r = by_block(ndf.jaccard(min_members=1))
    assert r.loc['a', 'members'] == ('A', 'C')
    assert r.loc['a', 'n_members'] == 2
    assert r.loc['z', 'members'] == ()
    assert r.loc['z', 'n_members'] == 0


def test_column_argument():
    """column= reads another column."""
    ndf = make([('a', 'ABC'), ('b', 'ABC')], column='pfam')
    r = ndf.jaccard('pfam')
    assert r.component.tolist() == [0, 0]
    assert r.members.tolist() == [('A', 'B', 'C')] * 2


def test_output_shape_and_dtypes():
    """Columns, order, index and dtypes are fixed."""
    r = fixture().jaccard()
    assert list(r.columns) == ['block_id', 'members', 'n_members', 'community', 'component']
    assert isinstance(r.index, pd.RangeIndex)
    assert r.n_members.dtype == 'int64'
    assert r.community.dtype == 'Int64'
    assert r.component.dtype == 'Int64'


def test_rows_in_order_of_first_appearance():
    """Interleaved blocks keep the order of their first row."""
    ndf = NeighborhoodDF(pd.DataFrame({
        'block_id': ['b', 'a', 'b', 'a', 'c'],
        'c80e3': ['X', 'Y', 'Z', 'Z', 'Q']}))
    assert ndf.jaccard(min_members=1).block_id.tolist() == ['b', 'a', 'c']


def test_empty_table():
    """An empty table gives an empty frame with the five columns."""
    ndf = NeighborhoodDF(pd.DataFrame({'block_id': [], 'c80e3': []}))
    r = ndf.jaccard()
    assert r.empty
    assert list(r.columns) == ['block_id', 'members', 'n_members', 'community', 'component']


def test_missing_column_raises_keyerror():
    """A missing column is a KeyError that names it."""
    ndf = NeighborhoodDF(pd.DataFrame({'block_id': [1], 'pfam': ['A']}))
    with pytest.raises(KeyError, match='c80e3'):
        ndf.jaccard()


def test_bad_parameters_raise_valueerror():
    """min_members below 1 and min_jaccard outside [0, 1] are rejected."""
    with pytest.raises(ValueError):
        fixture().jaccard(min_members=0)
    with pytest.raises(ValueError):
        fixture().jaccard(min_jaccard=1.5)


def test_self_unchanged():
    """The table is not modified."""
    ndf = fixture()
    before = ndf.copy()
    ndf.jaccard()
    pd.testing.assert_frame_equal(pd.DataFrame(ndf), pd.DataFrame(before))


def test_identical_blocks_alone_are_joined():
    """Two identical blocks with no other partner still form a group."""
    r = make([('a', 'ABC'), ('b', 'ABC'), ('c', 'XYZ')]).jaccard()
    assert r.community.tolist() == [0, 0, pd.NA]
    assert r.component.tolist() == [0, 0, pd.NA]
