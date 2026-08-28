#!/usr/bin/env python3
"""Tests for rotifer.pandas._group_key.

Ground truth is pandas itself: inside GroupBy.apply, pandas sets `.name` on each
group frame, and _group_key must reproduce that value exactly -- including the
shapes that are easy to get wrong:

  * a single-element list groupby, which pandas unwraps to a scalar (not a 1-tuple)
  * a tuple that is one multi-level column label, whose key stays scalar, versus a
    list of labels, whose key is a tuple
  * a list groupby, which is unhashable and must not reach an `in df.columns` test

It must also work standalone (no `.name`), and on empty frames -- which do reach
apply, e.g. a categorical grouper with observed=False.

Run under pytest, or standalone: python test/pandas/test_group_key.py
"""

import pandas as pd
import pytest

from rotifer.pandas import _group_key


def frames():
    flat = pd.DataFrame({'g': ['A', 'A', 'B'], 'h': ['x', 'x', 'y'],
                         'n': [1, 2, 3], 'v': [9, 8, 7]})
    mi = pd.DataFrame({('meta', 'g'): ['A', 'A', 'B'], ('meta', 'h'): ['x', 'x', 'y'],
                       ('data', 'v'): [1, 2, 3]})
    mi.columns = pd.MultiIndex.from_tuples(mi.columns)
    mi_int = pd.DataFrame({('m', 'g'): [1, 1, 2], ('d', 'v'): [1, 2, 3]})
    mi_int.columns = pd.MultiIndex.from_tuples(mi_int.columns)
    return flat, mi, mi_int


def groupby_specs():
    """(frame_name, groupby) pairs whose keys are checked against pandas' .name."""
    return [
        ('flat', 'g'), ('flat', ['g']), ('flat', ['g', 'h']),
        ('flat', ['g', 'h', 'n']), ('flat', 'n'),
        ('mi', ('meta', 'g')), ('mi', [('meta', 'g')]),
        ('mi', [('meta', 'g'), ('meta', 'h')]),
        ('mi_int', ('m', 'g')), ('mi_int', [('m', 'g')]),
    ]


@pytest.mark.parametrize('frame_name,groupby', groupby_specs())
def test_matches_pandas_name(frame_name, groupby):
    df = dict(zip(('flat', 'mi', 'mi_int'), frames()))[frame_name]
    seen = []
    df.groupby(groupby, sort=False).apply(
        lambda g: seen.append((getattr(g, 'name'), _group_key(g, groupby))))
    assert seen, 'no groups were produced'
    for name, got in seen:
        assert got == name
        # Tuple-vs-scalar is the shape that regresses (a 1-element list groupby
        # must not yield a 1-tuple). Exact dtype is not part of the contract:
        # pandas hands back a Python int where deriving from the frame yields
        # np.int64, and _fmt_label demotes numpy scalars before rendering anyway.
        assert isinstance(got, tuple) == isinstance(name, tuple)


@pytest.mark.parametrize('groupby,expected', [
    ('g', 'A'),                 # scalar label -> scalar key
    (['g'], 'A'),               # pandas unwraps single-column groupings
    (['g', 'h'], ('A', 'x')),   # sequence of labels -> tuple key
    ('nope', None),             # unknown column
])
def test_standalone_without_name(groupby, expected):
    """Outside GroupBy.apply there is no `.name`, so the key is derived from row 0."""
    flat = frames()[0]
    assert _group_key(flat, groupby) == expected


@pytest.mark.parametrize('groupby', ['g', ['g'], ['g', 'h']])
def test_empty_frame_returns_none(groupby):
    """Empty groups do reach apply (categorical grouper with observed=False), so
    deriving the key must not index into row 0."""
    flat = frames()[0]
    assert _group_key(flat.iloc[:0], groupby) is None


def test_list_groupby_is_not_hashed():
    """A list is unhashable; `groupby in df.columns` on one raises TypeError."""
    flat = frames()[0]
    assert _group_key(flat, ['g', 'h']) == ('A', 'x')


def test_name_wins_over_a_column_called_name():
    """pandas sets .name via object.__setattr__, so it shadows column attribute
    access -- a column literally named 'name' must not be mistaken for the key."""
    df = pd.DataFrame({'g': ['A', 'A'], 'name': ['zzz', 'yyy']})
    seen = []
    df.groupby('g', sort=False).apply(lambda d: seen.append(_group_key(d, 'g')))
    assert seen == ['A']


if __name__ == '__main__':
    import sys
    failures = 0
    named = dict(zip(('flat', 'mi', 'mi_int'), frames()))
    for frame_name, groupby in groupby_specs():
        df = named[frame_name]
        seen = []
        df.groupby(groupby, sort=False).apply(
            lambda g, gb=groupby: seen.append((getattr(g, 'name'), _group_key(g, gb))))
        for name, got in seen:
            ok = got == name and isinstance(got, tuple) == isinstance(name, tuple)
            failures += not ok
            print('  %-8s %-30s name=%-16r got=%-16r %s'
                  % (frame_name, str(groupby), name, got, '' if ok else '<-- MISMATCH'))
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
