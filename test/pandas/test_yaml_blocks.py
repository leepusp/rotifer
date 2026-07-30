#!/usr/bin/env python3
"""Tests for rotifer.pandas.yaml_blocks.

The contract is that a document written by to_yaml is (a) valid YAML that
yaml.safe_load parses without help, (b) readable and editable by hand, and
(c) restorable to the original frame by read_yaml.

Run under pytest, or standalone: python test/pandas/test_yaml_blocks.py
"""

import io

import pandas as pd
import pytest
import yaml

from rotifer.pandas import yaml_blocks as yb


def frame():
    return pd.DataFrame({'c80e3': ['A', 'A', 'B'], 'pid': ['p1', 'p2', 'p3'],
                         'plen': [10, 200, 3000],
                         'product': ['prod one', 'has | a pipe', 'p3'],
                         'lineage': ['L1', 'L1', 'L2']})


DOC_KWARGS = dict(columns=['pid', 'plen', 'product'],
                  stats={'lineage': 'nunique', 'pid': 'nunique'}, sortrows={})


def document(df=None, groupby='c80e3', **kwargs):
    return yb.to_yaml(df if df is not None else frame(), groupby=groupby,
                      **{**DOC_KWARGS, **kwargs})


# ------------------------------------------------------------------- structure
def test_document_is_valid_yaml():
    doc = yaml.safe_load(document())
    assert isinstance(doc, list)
    assert doc[0] == {'parameters': {'groupby': ['c80e3']}}


def test_groupby_is_stored_once():
    """The grouping columns live in the parameters element, not on every block."""
    text = document()
    assert text.count('groupby') == 1
    assert [b['key'] for b in yaml.safe_load(text)[1:]] == [['A'], ['B']]


def test_blocks_carry_key_and_tables():
    for block in yaml.safe_load(document())[1:]:
        assert 'key' in block
        assert {'stats', 'rows'} <= set(block)
        assert all(isinstance(v, str) for k, v in block.items() if k != 'key')


def test_separator_and_blank_lines():
    """One separator line per block, each preceded by a blank line; no blank lines
    inside a block."""
    lines = document().split('\n')
    seps = [i for i, l in enumerate(lines) if l.startswith('# ---')]
    assert len(seps) == 2
    for i in seps:
        assert lines[i - 1] == ''
    body = '\n'.join(lines[seps[0]:seps[1]])
    assert '\n\n' not in body.rstrip('\n')


# ---------------------------------------------------------------- round-tripping
def test_round_trip_restores_the_frame():
    got = yb.read_yaml(document(), merge=True, drop=['stats'])
    want = frame()[['c80e3', 'pid', 'plen', 'product']]
    assert got[want.columns.tolist()].equals(want)


def test_round_trip_multicolumn_groupby():
    df = pd.DataFrame({'g': ['A', 'A', 'B'], 'h': ['x', 'y', 'x'],
                       'pid': ['p1', 'p2', 'p3'], 'plen': [10, 200, 3000]})
    text = yb.to_yaml(df, groupby=['g', 'h'], columns=['pid', 'plen'],
                      stats={'pid': 'nunique'}, sortrows={})
    assert yaml.safe_load(text)[0] == {'parameters': {'groupby': ['g', 'h']}}
    got = yb.read_yaml(text, merge=True, drop=['stats'])
    want = df[['g', 'h', 'pid', 'plen']]
    assert got[want.columns.tolist()].equals(want)


def test_escaped_pipe_survives():
    rows = yb.read_yaml(document())['rows']
    assert 'has | a pipe' in rows['product'].tolist()


def test_newline_in_a_cell_is_flattened():
    """An embedded newline would terminate the YAML block scalar, so it is
    replaced on write rather than corrupting the document."""
    df = frame()
    df.loc[0, 'product'] = 'two\nlines'
    text = document(df)
    yaml.safe_load(text)                       # must still parse
    assert 'two lines' in yb.read_yaml(text)['rows']['product'].tolist()


def test_dtypes_are_inferred():
    rows = yb.read_yaml(document())['rows']
    assert rows['plen'].dtype.kind == 'i'
    assert rows['pid'].dtype == object


# ------------------------------------------------------------- combining tables
def test_merge_false_returns_one_frame_per_table():
    frames = yb.read_yaml(document())
    assert set(frames) == {'stats', 'rows'}
    assert len(frames['stats']) == 2         # one row per block, concatenated
    assert len(frames['rows']) == 3


def test_merge_broadcasts_block_values_onto_rows():
    merged = yb.read_yaml(document(), merge=True)
    assert len(merged) == 3
    assert merged.loc[merged.c80e3 == 'A', 'pid_stats'].tolist() == [2, 2]


def test_merge_suffix_disambiguates_colliding_names():
    """Stats are conventionally named after the columns they summarise, so `pid`
    means the identifier in one table and its distinct count in the other."""
    merged = yb.read_yaml(document(), merge=True)
    assert 'pid' in merged and 'pid_stats' in merged
    assert merged['pid'].dtype == object          # identifiers keep the plain name
    assert merged['pid_stats'].dtype.kind == 'i'


def test_drop_removes_named_tables():
    assert set(yb.read_yaml(document(), drop=['stats'])) == {'rows'}
    merged = yb.read_yaml(document(), merge=True, drop=['stats'])
    assert 'pid_stats' not in merged.columns


# -------------------------------------------------------------- editing by hand
def test_hand_edited_table_still_parses():
    """Padding is cosmetic: ragged pipes and stray spaces must not break reading,
    since the document is meant to be edited in a text editor."""
    text = '''- parameters:
    groupby: [g]

# ----
- key: [A]
  rows: |
    | pid | plen | product  |
    |---|---|---|
    | p1 |  10 | prod one |
    |p2|200|   sloppy but valid   |
'''
    rows = yb.read_yaml(text)['rows']
    assert rows['product'].tolist() == ['prod one', 'sloppy but valid']
    assert rows['plen'].tolist() == [10, 200]
    assert rows['g'].tolist() == ['A', 'A']


def test_groupby_override():
    frames = yb.read_yaml(document(), groupby=['c80e3'])
    assert 'c80e3' in frames['rows'].columns


# --------------------------------------------------------------------- sources
def test_reads_from_path_and_file_object(tmp_path):
    text = document()
    path = tmp_path / 'blocks.yml'
    yb.to_yaml(frame(), groupby='c80e3', buf=str(path), **DOC_KWARGS)
    assert path.read_text() == text
    assert yb.read_yaml(str(path))['rows'].equals(yb.read_yaml(text)['rows'])
    assert yb.read_yaml(io.StringIO(text))['rows'].equals(yb.read_yaml(text)['rows'])


def test_rejects_a_non_document():
    with pytest.raises(ValueError):
        yb.read_yaml('just: a mapping\n')


if __name__ == '__main__':
    import sys
    text = document()
    print(text)
    checks = [
        ('valid yaml', lambda: isinstance(yaml.safe_load(text), list)),
        ('groupby stored once', lambda: text.count('groupby') == 1),
        ('round trip', lambda: yb.read_yaml(text, merge=True, drop=['stats'])
            [['c80e3', 'pid', 'plen', 'product']].equals(frame()[['c80e3', 'pid', 'plen', 'product']])),
        ('escaped pipe', lambda: 'has | a pipe' in yb.read_yaml(text)['rows']['product'].tolist()),
        ('drop=[stats]', lambda: set(yb.read_yaml(text, drop=['stats'])) == {'rows'}),
        ('merge suffixes', lambda: 'pid_stats' in yb.read_yaml(text, merge=True).columns),
    ]
    failures = 0
    for name, check in checks:
        try:
            ok = check()
        except Exception as exc:
            ok, name = False, f'{name} ({type(exc).__name__}: {exc})'
        failures += not ok
        print('  %-20s %s' % (name, 'ok' if ok else 'FAILED'))
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
