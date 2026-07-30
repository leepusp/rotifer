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


# --------------------------------------------------------- caller-defined tables
def tables_frame():
    return pd.DataFrame({'g': ['A', 'A', 'B'], 'pid': ['p1', 'p2', 'p3'],
                         'plen': [10, 200, 3000],
                         'product': ['kinase', 'kinase', 'permease']})


def test_mapping_spec_becomes_an_aggregated_table():
    text = yb.to_yaml(tables_frame(), groupby='g', columns=['pid'],
                      stats={'pid': 'nunique'})
    block = yaml.safe_load(text)[1]
    assert 'stats' in block
    assert yb.read_md_table(block['stats'])['pid'].tolist() == [2]


def test_callable_spec_becomes_a_per_block_table():
    text = yb.to_yaml(tables_frame(), groupby='g', columns=['pid'],
                      products=lambda b: b['product'].value_counts().reset_index())
    block = yaml.safe_load(text)[1]
    assert yb.read_md_table(block['products'])['product'].tolist() == ['kinase']


def test_tables_keep_the_order_they_were_named_in():
    """Keyword order is insertion order, and the main table comes last."""
    text = yb.to_yaml(tables_frame(), groupby='g', columns=['pid'],
                      stats={'pid': 'nunique'}, span={'plen': 'max'},
                      products=lambda b: b[['product']].drop_duplicates())
    block = yaml.safe_load(text)[1]
    assert [k for k in block if k != 'key'] == ['stats', 'span', 'products', 'rows']


def test_any_aggregated_table_can_order_the_blocks():
    """sortblocks is not limited to a table called `stats`."""
    df = pd.DataFrame({'g': ['A', 'B', 'C'], 'pid': ['p1', 'p2', 'p3'],
                       'plen': [200, 3000, 5]})
    text = yb.to_yaml(df, groupby='g', columns=['pid'], span={'plen': 'max'},
                      sortblocks={'plen': True})
    assert block_keys(text) == ['C', 'A', 'B']


def test_unrelated_tables_do_not_interfere_when_sorting():
    """Tables are never joined to one another on write, so two of them defining
    the same column name is not an error; the one named first decides."""
    df = pd.DataFrame({'g': list('AABBCC'), 'pid': [f'p{i}' for i in range(6)],
                       'plen': [10, 200, 1, 3000, 5, 5]})
    # maxima are A=200, B=3000, C=5; minima are A=10, B=1, C=5
    by_max = yb.to_yaml(df, groupby='g', columns=['pid'], hi={'plen': 'max'},
                        lo={'plen': 'min'}, sortblocks={'plen': True})
    by_min = yb.to_yaml(df, groupby='g', columns=['pid'], lo={'plen': 'min'},
                        hi={'plen': 'max'}, sortblocks={'plen': True})
    assert block_keys(by_max) == ['C', 'A', 'B']
    assert block_keys(by_min) == ['B', 'C', 'A']


def test_an_aggregation_and_a_callback_table_coexist():
    df = pd.DataFrame({'g': ['A', 'B', 'C'], 'pid': ['p1', 'p2', 'p3'],
                       'plen': [200, 3000, 5]})
    text = yb.to_yaml(df, groupby='g', columns=['pid'], a={'plen': 'max'},
                      b=lambda block: block[['pid']], sortblocks={'plen': False})
    assert block_keys(text) == ['B', 'A', 'C']
    assert [k for k in yaml.safe_load(text)[1] if k != 'key'] == ['a', 'b', 'rows']


def test_the_main_table_can_be_renamed():
    text = yb.to_yaml(tables_frame(), groupby='g', columns=['pid'], name='neighbours')
    assert 'neighbours' in yaml.safe_load(text)[1]


@pytest.mark.parametrize('kwargs,error', [
    ({'key': {'plen': 'max'}}, ValueError),          # reserved for block metadata
    ({'parameters': {'plen': 'max'}}, ValueError),
    ({'oops': 42}, TypeError),                       # neither mapping nor callable
    ({'stats': {'nosuch': 'nunique'}}, KeyError),    # column not in the frame
])
def test_rejected_table_specs(kwargs, error):
    with pytest.raises(error):
        yb.to_yaml(tables_frame(), groupby='g', columns=['pid'], **kwargs)


# ------------------------------------------------------------------- ordering
def ordering_frame():
    return pd.DataFrame({'g': ['b', 'a', 'c', 'a'], 'n': [1, 2, 2, 3],
                         'pid': ['p1', 'p2', 'p3', 'p4']})


def block_keys(text):
    return [b['key'][0] for b in yaml.safe_load(text)[1:]]


def test_blocks_default_to_first_appearance():
    text = yb.to_yaml(ordering_frame(), groupby='g', columns=['pid'])
    assert block_keys(text) == ['b', 'a', 'c']


@pytest.mark.parametrize('ascending,expected', [(True, ['a', 'b', 'c']), (False, ['c', 'b', 'a'])])
def test_sortblocks_on_a_grouping_column(ascending, expected):
    """reset_index puts the grouping columns in reach, so they can order blocks
    even though they are the index of the block frame."""
    text = yb.to_yaml(ordering_frame(), groupby='g', columns=['pid'],
                      sortblocks={'g': ascending})
    assert block_keys(text) == expected


@pytest.mark.parametrize('ascending,expected', [(True, ['b', 'c', 'a']), (False, ['a', 'c', 'b'])])
def test_sortblocks_on_a_stats_column(ascending, expected):
    """Blocks order by a column produced by stats: b sums 1, c sums 2, a sums 5."""
    text = yb.to_yaml(ordering_frame(), groupby='g', columns=['pid'],
                      stats={'n': 'sum'}, sortblocks={'n': ascending})
    assert block_keys(text) == expected


def test_sortblocks_falls_through_to_a_later_key():
    """All blocks tie on the statistic, so the grouping column decides."""
    text = yb.to_yaml(ordering_frame(), groupby='g', columns=['pid'],
                      stats={'n': 'nunique'}, sortblocks={'n': False, 'g': True})
    assert block_keys(text) == ['a', 'b', 'c']


def test_sortblocks_ignores_unknown_columns():
    """One ordering can be reused across frames of differing shape."""
    text = yb.to_yaml(ordering_frame(), groupby='g', columns=['pid'],
                      sortblocks={'nosuch': True})
    assert block_keys(text) == ['b', 'a', 'c']


def test_sortrows_orders_within_a_block_only():
    """b appears first in the frame, but sorting globally by lineage would put an
    a row first -- block order must stay first-appearance."""
    df = pd.DataFrame({'g': ['b', 'a', 'a'], 'lineage': ['L2', 'L9', 'L1'],
                       'pid': ['p1', 'p2', 'p3']})
    text = yb.to_yaml(df, groupby='g', columns=['lineage', 'pid'],
                      sortrows={'lineage': True})
    assert block_keys(text) == ['b', 'a']
    rows = yb.read_yaml(text)['rows']
    assert rows[rows.g == 'a']['lineage'].tolist() == ['L1', 'L9']


def test_sortrows_descending():
    df = pd.DataFrame({'g': ['a', 'a', 'a'], 'n': [2, 3, 1], 'pid': ['p1', 'p2', 'p3']})
    text = yb.to_yaml(df, groupby='g', columns=['n'], sortrows={'n': False})
    assert yb.read_yaml(text)['rows']['n'].tolist() == [3, 2, 1]


# ------------------------------------------------------------- combining tables
def test_merge_false_returns_one_frame_per_table():
    frames = yb.read_yaml(document())
    assert set(frames) == {'stats', 'rows'}
    assert len(frames['stats']) == 2         # one row per block, concatenated
    assert len(frames['rows']) == 3


def test_tables_are_recombined_by_key_not_by_columns():
    """The key is what makes two tables equivalent. A `stats` table and a `span`
    table are different things even when their columns match exactly, so they are
    never concatenated onto one another."""
    df = pd.DataFrame({'g': list('AABC'), 'pid': ['p1', 'p2', 'p3', 'p4'],
                       'plen': [10, 200, 3000, 5]})
    text = yb.to_yaml(df, groupby='g', columns=['pid', 'plen'],
                      stats={'pid': 'nunique'}, span={'pid': 'max'})
    frames = yb.read_yaml(text)
    assert set(frames) == {'stats', 'span', 'rows'}
    assert list(frames['stats'].columns) == list(frames['span'].columns)  # same shape
    assert len(frames['stats']) == len(frames['span']) == 3               # still apart
    merged = yb.read_yaml(text, merge=True)
    assert {'pid', 'pid_stats', 'pid_span'} <= set(merged.columns)


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


def test_grouping_columns_come_from_the_document():
    """The reader is told nothing; the parameters element describes the document,
    and the block keys supply the values."""
    df = pd.DataFrame({'left': ['a', 'a', 'b'], 'right': ['x', 'y', 'x'],
                       'pid': ['p1', 'p2', 'p3']})
    frames = yb.read_yaml(yb.to_yaml(df, groupby=['left', 'right'], columns=['pid']))
    assert list(frames['rows'].columns[:2]) == ['left', 'right']
    assert frames['rows'].sort_values('pid')['left'].tolist() == ['a', 'a', 'b']
    assert frames['rows'].sort_values('pid')['right'].tolist() == ['x', 'y', 'x']


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
