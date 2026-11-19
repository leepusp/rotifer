#!/usr/bin/env python3
"""Tests for the group key rendered in rotifer.pandas block headers.

Every block header reads `# <groupby>=<key>; <stat>=<value>; ...`. The key is the
part under test here, observed only through the public to_blocks and
default_header output. It follows what pandas hands to GroupBy.apply as `.name`,
including the shapes that are easy to get wrong:

  * a single-element list groupby, which pandas unwraps to a scalar, is rendered
    as a 1-tuple so that label and key keep matching arity
  * a list of labels renders a tuple key, scalars stay bare
  * a list groupby is unhashable and must not reach a column-membership test

Outside GroupBy.apply (default_header called directly) there is no `.name`, so the
key is derived from the first row; on an empty frame, or an unknown column, it
renders as None.

Run under pytest, or standalone: python test/pandas/test_group_key.py
"""

import numpy as np
import pandas as pd
import pytest

import rotifer.pandas as rp


def frames():
    flat = pd.DataFrame({'g': ['A', 'A', 'B'], 'h': ['x', 'x', 'y'],
                         'n': [1, 2, 3], 'v': [9, 8, 7]})
    mi = pd.DataFrame({('meta', 'g'): ['A', 'A', 'B'], ('meta', 'h'): ['x', 'x', 'y'],
                       ('data', 'v'): [1, 2, 3]})
    mi.columns = pd.MultiIndex.from_tuples(mi.columns)
    mi_int = pd.DataFrame({('m', 'g'): [1, 1, 2], ('d', 'v'): [1, 2, 3]})
    mi_int.columns = pd.MultiIndex.from_tuples(mi_int.columns)
    return flat, mi, mi_int


STATS = {'flat': {'v': 'nunique'}, 'mi': {('data', 'v'): 'nunique'},
         'mi_int': {('d', 'v'): 'nunique'}}


def keys_in_headers(frame_name, groupby):
    """The `<key>` text of every block header, in block order."""
    df = dict(zip(('flat', 'mi', 'mi_int'), frames()))[frame_name]
    text = rp.to_blocks(df, groupby=groupby, stats=STATS[frame_name], columns=[], sortrows={})
    heads = [ln for ln in text.splitlines() if ln.startswith('# ')]
    return [h[2:].split('; ')[0].split('=', 1)[1] for h in heads]


@pytest.mark.parametrize('frame_name,groupby,expected', [
    ('flat', 'g', ['A', 'B']),
    ('flat', ['g'], ["('A',)", "('B',)"]),              # unwrapped by pandas, padded back
    ('flat', ['g', 'h'], ["('A', 'x')", "('B', 'y')"]),
    ('flat', ['g', 'h', 'n'], ["('A', 'x', 1)", "('A', 'x', 2)", "('B', 'y', 3)"]),
    ('flat', 'n', ['1', '2', '3']),                      # numpy ints render bare
    # One multi-level column label: within to_blocks the grouping column is gone
    # from the frame, so the key is padded to a 1-tuple (default_header called
    # directly on a frame that still holds the column keeps it scalar; see
    # test_to_blocks.test_multilevel_column_label_keeps_a_scalar_key).
    ('mi', ('meta', 'g'), ["('A',)", "('B',)"]),
    ('mi', [('meta', 'g')], ["('A',)", "('B',)"]),
    ('mi', [('meta', 'g'), ('meta', 'h')], ["('A', 'x')", "('B', 'y')"]),
    ('mi_int', ('m', 'g'), ['(1,)', '(2,)']),
    ('mi_int', [('m', 'g')], ['(1,)', '(2,)']),
])
def test_key_rendered_in_block_headers(frame_name, groupby, expected):
    assert keys_in_headers(frame_name, groupby) == expected


@pytest.mark.parametrize('groupby,expected', [
    ('g', '#\n# g=A; v=24'),                    # scalar label -> scalar key
    (['g'], "#\n# ('g',)=('A',); v=24"),         # pandas unwraps; padded back to match
    (['g', 'h'], "#\n# ('g', 'h')=('A', 'x'); v=24"),
    ('nope', '#\n# nope=None; v=24'),            # unknown column
])
def test_standalone_header_derives_key_from_first_row(groupby, expected):
    """Outside GroupBy.apply there is no `.name`, so the key comes from row 0."""
    flat = frames()[0]
    assert rp.default_header(flat, groupby=groupby, stats={'v': 'sum'}) == expected


def test_empty_frame_renders_a_none_key():
    """Empty groups do reach apply (categorical grouper with observed=False), so
    deriving the key must not index into row 0."""
    flat = frames()[0]
    assert rp.default_header(flat.iloc[:0], groupby='g', stats={'v': 'sum'}) == '#\n# g=None; v=0'
    assert rp.default_header(flat.iloc[:0], groupby=['g'], stats={'v': 'sum'}) \
        == "#\n# ('g',)=(None,); v=0"


def test_list_groupby_is_not_hashed():
    """A list is unhashable; `groupby in df.columns` on one raises TypeError."""
    flat = frames()[0]
    assert rp.default_header(flat, groupby=['g', 'h'], stats={'v': 'sum'}).endswith("=('A', 'x'); v=24")


def test_numpy_scalars_render_as_python_literals():
    """repr(np.int64(1)) is not a literal under numpy 2, so keys inside a tuple are
    demoted to builtins before rendering."""
    flat = frames()[0]
    head = rp.default_header(flat.iloc[:1], groupby=['n', 'g'], stats={'v': 'sum'})
    assert head == "#\n# ('n', 'g')=(1, 'A'); v=9"
    assert 'np.' not in head and 'numpy' not in head
    assert rp.default_header(flat, groupby='g', stats={'n': np.int64(7)}) == '#\n# g=A; n=7'


def test_name_wins_over_a_column_called_name():
    """pandas sets .name via object.__setattr__, so it shadows column attribute
    access -- a column literally named 'name' must not be mistaken for the key."""
    df = pd.DataFrame({'g': ['A', 'A'], 'name': ['zzz', 'yyy']})
    text = rp.to_blocks(df, groupby='g', stats={'name': 'nunique'}, columns=['name'], sortrows={})
    assert text == '#\n# g=A; name=2\n  : name\n0 : zzz \n1 : yyy '


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
