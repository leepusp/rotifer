"""Unit tests for the read_yaml pipeline steps and the shared table reader."""

import numpy as np
import pandas as pd
import pytest

from rotifer.pandas import yaml_blocks as yb


def doc(*blocks, groupby='[g]'):
    out = [f'- parameters:\n    groupby: {groupby}']
    for key, rows in blocks:
        out.append(f'- key: [{key}]\n  rows: |\n' + '\n'.join('    ' + r for r in rows))
    return '\n'.join(out) + '\n'


def read(*rows_by_block, **kw):
    blocks = [(k, ['| a | b |', '| - | - |'] + r) for k, r in rows_by_block]
    return yb.read_yaml(doc(*blocks), **kw)


# ----------------------------------------------------------------- _merge_plan
FRAMES = {'x': pd.DataFrame(), 'y': pd.DataFrame(), 'z': pd.DataFrame()}


def test_plan_true_takes_every_table_in_reverse_on_the_group_columns():
    assert yb._merge_plan(True, FRAMES, ['g']) == {'z': ['g'], 'y': ['g'], 'x': ['g']}
    assert list(yb._merge_plan(True, FRAMES, ['g'])) == ['z', 'y', 'x']


@pytest.mark.parametrize('merge', [None, False])
def test_plan_none_and_false_mean_no_merge(merge):
    assert yb._merge_plan(merge, FRAMES, ['g']) == {}


def test_plan_string_and_list_use_group_columns_in_given_order():
    assert yb._merge_plan('y', FRAMES, ['g']) == {'y': ['g']}
    plan = yb._merge_plan(['y', 'x'], FRAMES, ['g'])
    assert plan == {'y': ['g'], 'x': ['g']} and list(plan) == ['y', 'x']


def test_plan_mapping_normalises_each_join_to_a_list():
    plan = yb._merge_plan({'x': 'id', 'y': ['id', 'k']}, FRAMES, ['g'])
    assert plan == {'x': ['id'], 'y': ['id', 'k']}


@pytest.mark.parametrize('merge', ['nope', ['x', 'nope'], {'nope': 'g'}])
def test_plan_unknown_table_is_a_key_error(merge):
    with pytest.raises(KeyError, match='nope'):
        yb._merge_plan(merge, FRAMES, ['g'])


@pytest.mark.parametrize('merge', [[], {}])
def test_plan_empty_is_a_value_error(merge):
    with pytest.raises(ValueError, match='merge names no tables'):
        yb._merge_plan(merge, FRAMES, ['g'])


# ------------------------------------------------------------ typing of columns
def test_column_numeric_except_one_value_stays_text():
    df = read(('A', ['| 1 | 1 |', '| 2 | x |']), ('B', ['| 3 | 4 |']))
    assert df['a'].dtype == np.int64
    assert not pd.api.types.is_numeric_dtype(df['b'])
    assert df['b'].tolist() == ['1', 'x', '4']


def test_int_and_float_typing():
    df = read(('A', ['| 1 | 1.5 |', '| 2 | 2 |']))
    assert df['a'].dtype == np.int64
    assert df['b'].dtype == np.float64 and df['b'].tolist() == [1.5, 2.0]


def test_missing_cells_become_missing_and_int_columns_turn_float():
    df = read(('A', ['| 1 |   |', '|   | y |']))
    assert df['a'].dtype == np.float64 and df['a'].isna().tolist() == [False, True]
    assert df['b'].isna().tolist() == [True, False] and df['b'].iloc[1] == 'y'


def test_all_missing_column_stays_object_of_none():
    df = read(('A', ['| 1 |   |', '| 2 |   |']))
    assert df['b'].dtype == object and df['b'].isna().all()


def test_escaped_pipe_only_in_some_columns():
    df = read(('A', [r'| p \| q | 1 |', '| r | 2 |']), ('B', ['| s | 3 |']))
    assert df['a'].tolist() == ['p | q', 'r', 's']
    assert df['b'].tolist() == [1, 2, 3]


def test_group_columns_come_first_and_keep_their_key_type():
    df = read((7, ['| u | 1 |']), (8, ['| v | 2 |']))
    assert list(df.columns)[0] == 'g' and df['g'].tolist() == [7, 8]


def test_same_key_across_blocks_concatenates_rows():
    df = read(('A', ['| 1 | 2 |']), ('B', ['| 3 | 4 |']))
    assert df['a'].tolist() == [1, 3]


def test_many_rows_cross_the_chunk_boundary(monkeypatch):
    monkeypatch.setattr(yb, '_CHUNK', 3)
    rows = [f'| {i} | t{i} |' for i in range(10)]
    df = read(('A', rows[:4]), ('B', rows[4:]))
    assert df['a'].tolist() == list(range(10))
    assert df['b'].tolist() == [f't{i}' for i in range(10)]


# ---------------------------------------------------------------- read_md_table
def test_md_table_ragged_padding_and_typing():
    text = '\n'.join([
        'noise before',
        '|id   |  name |score|',
        '|:--  |  ---- |---: |',
        '|   1 | alpha   |  1.5 |',
        '| 22|b \\| c|  |',
        '|3 | d |4|',
        'noise after',
    ])
    df = yb.read_md_table(text)
    assert list(df.columns) == ['id', 'name', 'score']
    assert df['id'].tolist() == [1, 22, 3]
    assert df['name'].tolist() == ['alpha', 'b | c', 'd']
    assert df['score'].isna().tolist() == [False, True, False]
    assert df['score'].dtype == np.float64


def test_md_table_matches_read_yaml_on_the_same_table():
    table = ['| a | b |', '| - | - |', '| 1 | x \\| y |', '| 2 | z |']
    via_yaml = yb.read_yaml(doc(('A', table)))
    via_md = yb.read_md_table('\n'.join(table))
    pd.testing.assert_frame_equal(via_yaml.drop(columns='g'), via_md)


def test_md_table_short_row_is_padded_and_long_row_rejected():
    df = yb.read_md_table('| a | b |\n| - | - |\n| 1 |\n')
    assert df['a'].tolist() == [1] and df['b'].isna().all()
    with pytest.raises(ValueError):
        yb.read_md_table('| a |\n| - |\n| 1 | 2 |\n')


# ---------------------------------------------------------- escaping and hex
def test_pipe_in_an_aggregate_round_trips():
    df = pd.DataFrame({'g': ['A', 'A'], 's': ['a|b'] * 2, 'n': [1, 2]})
    out = yb.read_yaml(yb.to_yaml(df, 'g', stats={'s': 'max', 'n': 'sum'}), merge=None)
    assert out['stats'].loc[0, 's'] == 'a|b'
    assert out['stats'].loc[0, 'n'] == 3


def test_pipe_and_newline_in_a_column_name_round_trip():
    df = pd.DataFrame({'g': ['A', 'B'], 'x|y\nz': [1, 2], 'k': ['p', 'q']})
    out = yb.read_yaml(yb.to_yaml(df, 'g'), merge=None)['rows']
    assert list(out.columns) == ['g', 'x|y z', 'k']
    assert out['x|y z'].tolist() == [1, 2]
    assert out['k'].tolist() == ['p', 'q']


def test_pipe_and_newline_in_an_aggregate_name_round_trip():
    df = pd.DataFrame({'g': ['A', 'A'], 'v': [1, 2]})
    out = yb.read_yaml(yb.to_yaml(df, 'g', stats={'a|b\nc': ('v', 'sum')}), merge=None)
    assert out['stats'].loc[0, 'a|b c'] == 3


@pytest.mark.parametrize('ch', ['\r', '\x85', ' ', ' ', '\r\n'])
def test_line_break_characters_in_cells_become_a_space(ch):
    df = pd.DataFrame({'g': ['A', 'A'], 't': [f'one{ch}two', 'three']})
    out = yb.read_yaml(yb.to_yaml(df, 'g'), merge=None)['rows']
    assert out['t'].tolist() == ['one two', 'three']


def test_line_break_in_an_aggregate_becomes_a_space():
    df = pd.DataFrame({'g': ['A', 'A'], 's': ['a\rb', 'a\rb']})
    out = yb.read_yaml(yb.to_yaml(df, 'g', stats={'s': 'max'}), merge=None)
    assert out['stats'].loc[0, 's'] == 'a b'


def test_hexadecimal_column_stays_text():
    df = pd.DataFrame({'g': ['A', 'A'], 'h': ['0x10', '1']})
    out = yb.read_yaml(yb.to_yaml(df, 'g'), merge=None)['rows']
    assert out['h'].tolist() == ['0x10', '1']


def test_signed_integers_are_still_numbers():
    df = pd.DataFrame({'g': ['A', 'A'], 'i': ['+3', '-4']})
    out = yb.read_yaml(yb.to_yaml(df, 'g'), merge=None)['rows']
    assert out['i'].tolist() == [3, -4]
