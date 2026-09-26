#!/usr/bin/env python3
"""Tests for rotifer.pandas.to_blocks headers and rotifer.pandas._resolve_stats.

Block headers have the form

    # <groupby>=<key>; <stat>=<value>; ...

A scalar groupby stays bare (`# c80e3=A`) so the common case reads cleanly, while
array-like labels are emitted as Python literals (`# ('c80e3', 'pid')=('A', 'p1')`)
so a downstream parser can recover them with ast.literal_eval. Label and key
always have matching arity, so the two tuples can be zipped.

Run under pytest, or standalone: python test/pandas/test_to_blocks.py
"""

import ast

import pandas as pd
import pytest

import rotifer.pandas as rp
from rotifer.pandas import _resolve_stats


def frame():
    return pd.DataFrame({'c80e3': ['A', 'A', 'B'], 'pid': ['p1', 'p2', 'p3'],
                         'lineage': ['L1', 'L1', 'L2'], 'block_id': [1, 1, 2]})


BLOCK_KWARGS = dict(columns=['pid', 'lineage'], stats={'lineage': 'nunique'}, sortrows={})


def parse_label(text):
    """Recover a label emitted by _fmt_label: tuples arrive as Python literals."""
    text = text.strip()
    if text.startswith('(') and text.endswith(')'):
        return ast.literal_eval(text)
    return text


def parse_header(line):
    """'# gb=key; s1=v1; s2=v2' -> (groupby, key, {stat: value})"""
    fields = line.lstrip('#').strip().split('; ')
    groupby, key = fields[0].split('=', 1)
    stats = dict(f.split('=', 1) for f in fields[1:] if '=' in f)
    return parse_label(groupby), parse_label(key), {parse_label(k): v for k, v in stats.items()}


def first_header(df, **kwargs):
    return rp.to_blocks(df, **kwargs).splitlines()[1]


@pytest.mark.parametrize('groupby,label,key', [
    ('c80e3',            'c80e3',            'A'),
    (['c80e3'],          ('c80e3',),         ('A',)),
    (['c80e3', 'pid'],   ('c80e3', 'pid'),   ('A', 'p1')),
    (['block_id'],       ('block_id',),      (1,)),
])
def test_header_round_trips(groupby, label, key):
    """Both value and type must survive: an int key has to come back an int, which
    is why numpy scalars are demoted before repr -- repr(np.int64(1)) is not a
    literal under numpy 2 and would not parse."""
    got_label, got_key, _ = parse_header(first_header(frame(), groupby=groupby, **BLOCK_KWARGS))
    assert got_label == label and type(got_label) is type(label)
    assert got_key == key and type(got_key) is type(key)


def test_multilevel_column_label_keeps_a_scalar_key():
    """('m','g') is one column label here, not two labels, so its key stays scalar
    and must not be padded to a 1-tuple."""
    mi = pd.DataFrame({('m', 'g'): ['A', 'A', 'B'], ('m', 'n'): [1, 1, 2], ('d', 'v'): [1, 2, 3]})
    mi.columns = pd.MultiIndex.from_tuples(mi.columns)
    header = rp.default_header(mi[mi[('m', 'g')] == 'A'],
                               groupby=('m', 'g'), stats={('m', 'n'): 'nunique'})
    label, key, stats = parse_header(header.splitlines()[1])
    assert label == ('m', 'g')
    assert key == 'A'
    assert ('m', 'n') in stats


def test_scalar_groupby_header_is_unquoted():
    assert first_header(frame(), groupby='c80e3', **BLOCK_KWARGS) == '# c80e3=A; lineage=1'


def test_block_body_is_independent_of_groupby_shape():
    """The groupby shape changes only the header, never the rendered rows."""
    scalar = rp.to_blocks(frame(), groupby='c80e3', **BLOCK_KWARGS).splitlines()[2:4]
    listed = rp.to_blocks(frame(), groupby=['c80e3'], **BLOCK_KWARGS).splitlines()[2:4]
    assert scalar == listed


def test_missing_stats_column_raises():
    """A column that is not in the frame is unambiguous, so it fails loudly rather
    than emitting the aggregation name in place of a value."""
    with pytest.raises(KeyError):
        rp.to_blocks(frame(), groupby='c80e3', columns=['pid'], sortrows={},
                     stats={('lineage',): 'nunique'})


def test_resolve_stats_aggregations():
    df = pd.DataFrame({'organism': ['E. coli', 'X'], 'lineage': ['L1', 'L2'], 'n': [1, 2]})
    assert _resolve_stats(df, {'lineage': 'nunique'}) == {'lineage': 2}
    assert _resolve_stats(df, {'n': 'size'}) == {'n': 2}     # a property, not a method
    assert _resolve_stats(df, {'n': 7}) == {'n': 7}          # value supplied directly


def test_resolve_stats_is_idempotent_for_string_valued_aggregations():
    """_resolve_stats runs twice per block (default_block_rows, then default_header).
    An aggregation may itself return a string, which the second pass must not try
    to re-resolve."""
    df = pd.DataFrame({'organism': ['E. coli', 'X']})
    once = _resolve_stats(df, {'organism': 'max'})
    assert once == {'organism': 'X'}
    assert _resolve_stats(df, once) == once


if __name__ == '__main__':
    import sys
    failures = 0
    for groupby, label, key in [('c80e3', 'c80e3', 'A'), (['c80e3'], ('c80e3',), ('A',)),
                                (['c80e3', 'pid'], ('c80e3', 'pid'), ('A', 'p1')),
                                (['block_id'], ('block_id',), (1,))]:
        line = first_header(frame(), groupby=groupby, **BLOCK_KWARGS)
        got_label, got_key, _ = parse_header(line)
        ok = (got_label == label and type(got_label) is type(label)
              and got_key == key and type(got_key) is type(key))
        failures += not ok
        print('  %-20s %-34s -> label=%-18r key=%-14r %s'
              % (str(groupby), line, got_label, got_key, '' if ok else '<-- MISMATCH'))
    try:
        rp.to_blocks(frame(), groupby='c80e3', columns=['pid'], sortrows={},
                     stats={('lineage',): 'nunique'})
        failures += 1
        print('  bad stats key        no KeyError raised <-- MISMATCH')
    except KeyError as e:
        print('  bad stats key        KeyError: %s' % e)
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
