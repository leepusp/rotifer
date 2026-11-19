#!/usr/bin/env python3
"""Tests for the shared block engine, plan_blocks, and the to_blocks built on it."""

import numpy as np
import pandas as pd
import pytest

import rotifer.pandas as rp
from rotifer.pandas.blocks import plan_blocks, BlockPlan


def frame():
    return pd.DataFrame({'g': ['B', 'A', 'B', 'C', 'A'], 'h': ['x', 'y', 'x', 'z', 'y'],
                         'n': [3, 10, 200, 4, 5], 'v': [1.5, 2.5, 3.5, 4.5, 5.5],
                         's': ['q', 'p', 'r', 's', 't']})


def test_blocks_come_out_in_first_appearance_order():
    plan = plan_blocks(frame(), 'g')
    assert isinstance(plan, BlockPlan)
    assert plan.keys == ['B', 'A', 'C']
    assert plan.groupby == ['g']
    assert {k: list(v) for k, v in plan.positions.items()} == {'B': [0, 2], 'A': [1, 4], 'C': [3]}


def test_list_groupby_gives_tuple_keys_of_matching_arity():
    assert plan_blocks(frame(), ['g']).keys == [('B',), ('A',), ('C',)]
    assert plan_blocks(frame(), ['g', 'h']).keys == [('B', 'x'), ('A', 'y'), ('C', 'z')]


def test_sortblocks_by_a_grouping_column_and_by_an_aggregate():
    assert plan_blocks(frame(), 'g', sortblocks={'g': True}).keys == ['A', 'B', 'C']
    assert plan_blocks(frame(), 'g', sortblocks={'g': False}).keys == ['C', 'B', 'A']
    plan = plan_blocks(frame(), 'g', sortblocks={'n': False}, tables={'stats': {'n': 'sum'}})
    assert plan.keys == ['B', 'A', 'C']          # 203, 15, 4
    plan = plan_blocks(frame(), 'g', sortblocks={'n': True}, tables={'stats': {'n': 'sum'}})
    assert plan.keys == ['C', 'A', 'B']
    # the tables follow the blocks
    assert list(plan.tables['stats']['n']) == [4, 15, 203]


def test_sortblocks_ignores_unknown_columns():
    assert plan_blocks(frame(), 'g', sortblocks={'nope': True}).keys == ['B', 'A', 'C']


def test_sortrows_orders_rows_inside_blocks_only():
    plan = plan_blocks(frame(), 'g', sortrows={'n': False})
    assert plan.keys == ['B', 'A', 'C']
    assert list(plan.positions['B']) == [2, 0]
    assert list(plan.positions['A']) == [1, 4]


def test_sampling_is_reproducible_with_a_seeded_rng():
    df = pd.DataFrame({'g': ['a'] * 20 + ['b'] * 3, 'n': range(23)})
    one = plan_blocks(df, 'g', sample=5, rng=np.random.default_rng(7))
    two = plan_blocks(df, 'g', sample=5, rng=np.random.default_rng(7))
    assert len(one.positions['a']) == 5 and len(one.positions['b']) == 3
    assert list(one.positions['a']) == list(two.positions['a'])
    assert set(one.positions['a']) <= set(range(20))


def test_sampled_rows_are_still_sorted_by_sortrows():
    df = pd.DataFrame({'g': ['a'] * 20, 'n': range(20)})
    plan = plan_blocks(df, 'g', sample=6, sortrows={'n': False}, rng=np.random.default_rng(0))
    got = list(plan.positions['a'])
    assert got == sorted(got, reverse=True)


def test_aggregations_by_column_by_named_tuple_and_by_value():
    plan = plan_blocks(frame(), 'g', tables={
        'stats': {'n': 'sum', 'size': ('s', 'size'), 'top': ('v', lambda x: x.max())}})
    seven = plan_blocks(frame(), 'g', tables={'stats': {'n': 7}}).tables['stats']
    assert list(seven['n']) == [7, 7, 7]
    t = plan.tables['stats']
    assert list(t.columns) == ['n', 'size', 'top']
    assert list(t['n']) == [203, 15, 4] and t['n'].dtype.kind == 'i'
    assert list(t['size']) == [2, 2, 1]
    assert list(t['top']) == [3.5, 5.5, 4.5]


def test_callable_tables_become_callbacks():
    plan = plan_blocks(frame(), 'g', tables={'names': lambda b: b[['s']]})
    assert plan.tables == {} and list(plan.callbacks) == ['names']


def test_tables_keep_the_order_the_caller_named_them():
    plan = plan_blocks(frame(), 'g', tables={'b': {'n': 'sum'}, 'a': {'v': 'max'}})
    assert list(plan.tables) == ['b', 'a']


@pytest.mark.parametrize('name', ['key', 'parameters', 'merged'])
def test_reserved_table_names(name):
    with pytest.raises(ValueError, match='reserved'):
        plan_blocks(frame(), 'g', tables={name: {'n': 'sum'}})


def test_missing_column_raises_keyerror():
    with pytest.raises(KeyError, match='nope'):
        plan_blocks(frame(), 'g', tables={'stats': {'nope': 'sum'}})
    with pytest.raises(KeyError, match='nope'):
        plan_blocks(frame(), 'g', tables={'stats': {'x': ('nope', 'sum')}})


def test_bad_table_spec_raises_typeerror():
    with pytest.raises(TypeError):
        plan_blocks(frame(), 'g', tables={'stats': 3})


def test_to_blocks_rejects_unknown_keywords():
    with pytest.raises(TypeError, match='bogus'):
        rp.to_blocks(frame(), groupby='g', columns=['n'], stats={}, sortrows={}, bogus=1)


def test_to_blocks_header_hook_receives_key_groupby_and_stats():
    seen = []

    def header(key, groupby, stats):
        seen.append((key, groupby, stats))
        return f'## {key}'

    text = rp.to_blocks(frame(), groupby=['g'], columns=['n'], stats={'n': 'sum'},
                        sortrows={}, header=header)
    assert seen == [(('B',), ('g',), {'n': 203}), (('A',), ('g',), {'n': 15}),
                    (('C',), ('g',), {'n': 4})]
    assert text.startswith("## ('B',)\n  : n  ")


def test_to_blocks_sortblocks_may_use_a_stat():
    text = rp.to_blocks(frame(), groupby='g', columns=['n'], stats={'n': 'sum'},
                        sortrows={}, sortblocks={'n': True})
    assert [ln for ln in text.splitlines() if ln.startswith('# ')] == \
        ['# g=C; n=4', '# g=A; n=15', '# g=B; n=203']


def test_to_blocks_insert_option_goes_through_to_string():
    text = rp.to_blocks(frame(), groupby='g', columns=['n'], stats={'n': 'sum'}, sortrows={},
                        insert=[(0, 'tag', 'T')])
    assert text.splitlines()[2].split(' : ')[:2] == [' ', 'tag']


# --------------------------------------------------------- deprecated paths
def test_default_header_called_with_a_dataframe_warns_and_still_works():
    with pytest.warns(FutureWarning, match='2026-11-30'):
        out = rp.default_header(frame(), groupby='g', stats={'n': 'sum'})
    assert out == '#\n# g=B; n=222'


def test_default_header_new_signature():
    assert rp.default_header('A', 'g', {'n': 1, 'x': 2}) == '#\n# g=A; n=1; x=2'
    assert rp.default_header(('A', 'x'), ('g', 'h'), None) == "#\n# ('g', 'h')=('A', 'x')"


def test_default_block_rows_warns():
    with pytest.warns(FutureWarning, match='2026-11-30'):
        out = rp.default_block_rows(frame().drop(columns='g'), groupby='g', columns=['n'],
                                    stats={'n': 'sum'}, sortrows={})
    assert list(out.columns) == ['blocks', 'n']


def test_to_blocks_apply_warns_and_runs_the_old_algorithm():
    df = frame()
    with pytest.warns(FutureWarning, match='apply'):
        old = rp.to_blocks(df, groupby='g', columns=['n'], stats={'n': 'sum'},
                           sortrows={}, apply=_old_apply)
    assert old == 'X B\nX A\nX C'


def _old_apply(block, groupby, header, colsep, **kwargs):
    return pd.DataFrame({'blocks': [f'X {block.name}']})


# ------------------------------------------------------- genome-neutral defaults
def test_groupby_is_required():
    with pytest.raises(TypeError):
        rp.to_blocks(frame())
    with pytest.raises(TypeError):
        rp.to_yaml(frame())


def test_defaults_show_every_column_but_the_grouping_ones_and_no_stats():
    text = rp.to_blocks(frame()[['g', 'n', 's']], groupby='g')
    assert text.splitlines()[:4] == ['#', '# g=B', '  : n   : s', '0 : 3   : q']
    assert rp.to_blocks(frame(), groupby=['g', 'h'], columns=['n']).splitlines()[1] \
        == "# ('g', 'h')=('B', 'x')"


def test_no_sorting_by_default():
    text = rp.to_blocks(frame(), groupby='g', columns=['n'])
    assert [ln for ln in text.splitlines() if ln.startswith('# ')] == ['# g=B', '# g=A', '# g=C']
    assert text.splitlines()[3:5] == ['0 : 3  ', '2 : 200']


def test_the_package_names_no_genome_columns():
    import pathlib
    import rotifer.pandas
    root = pathlib.Path(rotifer.pandas.__file__).parent
    for path in root.glob('*.py'):
        if path.name in ('network.py',):         # to_network is a separate concern
            continue
        for word in ('c80e3', 'lineage', 'c80i70'):
            assert word not in path.read_text(), (path.name, word)
