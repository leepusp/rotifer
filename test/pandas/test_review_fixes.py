"""Regression tests for the verified code-review findings in rotifer.pandas."""

import numpy as np
import pandas as pd
import pytest

from rotifer.pandas import to_blocks, to_yaml, to_string, pad, to_color, default_header
from rotifer.pandas import blocks as rblocks
from rotifer.pandas.blocks import plan_blocks


def frame():
    return pd.DataFrame({'grp': ['A', 'A', 'B', 'B', 'B'], 'v': [1, 2, 3, 4, 5],
                         'w': list('abcde')})


# 1 --------------------------------------------------------------------------
def test_datetime_and_timedelta_stats_print_as_timestamps():
    df = pd.DataFrame({'grp': ['A', 'A'],
                       't': pd.to_datetime(['2021-01-01', '2021-01-02']),
                       'd': pd.to_timedelta([1, 3], unit='s'), 'n': [1, 2]})
    stats = {'t': 'min', 'd': 'max', 'n': 'max'}
    out = to_blocks(df, 'grp', stats=stats)
    assert '; t=2021-01-01 00:00:00; d=0 days 00:00:03; n=2' in out
    out = to_yaml(df, 'grp', stats=stats)
    assert '2021-01-01 00:00:00' in out and '0 days 00:00:03' in out
    assert '2.0' not in out.split('rows:')[0]


# 2 --------------------------------------------------------------------------
def test_sample_is_seedable_in_both_writers():
    df = pd.DataFrame({'grp': ['A'] * 50, 'v': range(50)})
    for write in (to_blocks, to_yaml):
        a = write(df, 'grp', sample=5, random_state=7)
        assert a == write(df, 'grp', sample=5, random_state=7)
        assert a == write(df, 'grp', sample=5, random_state=np.random.default_rng(7))
        assert any(a != write(df, 'grp', sample=5, random_state=s) for s in range(8, 14))


# 3 --------------------------------------------------------------------------
def test_apply_path_forwards_unknown_kwargs():
    seen = {}

    def apply(block, groupby, header, colsep, **kw):
        seen.update(kw)
        return pd.DataFrame({'blocks': ['x']})
    with pytest.warns(FutureWarning):
        assert to_blocks(frame(), 'grp', apply=apply, extra=3) == 'x\nx'
    assert seen == {'extra': 3}
    with pytest.raises(TypeError):
        to_blocks(frame(), 'grp', extra=3)


# 4 --------------------------------------------------------------------------
def test_nan_group_keys_are_dropped():
    df = pd.DataFrame({'grp': ['A', None, 'A'], 'v': [1, 2, 3]})
    plan = plan_blocks(df, 'grp')
    assert plan.keys == ['A'] and len(plan.positions['A']) == 2
    out = to_blocks(df, 'grp')
    assert out.count('# grp=') == 1 and len(out.splitlines()) == 5


# 5 --------------------------------------------------------------------------
def test_duplicate_column_names():
    df = pd.DataFrame([['A', 1, 2], ['A', 3, 4]], columns=['grp', 'x', 'x'])
    text = to_string(df.iloc[:, 1:])
    assert text.split('\n')[0] == '  : x : x'
    out = to_blocks(df, 'grp')
    assert out.endswith(text)


# 6 --------------------------------------------------------------------------
@pytest.mark.parametrize('align', ['left', 'right', 'center'])
def test_header_aligned_like_values(align):
    just = {'left': str.ljust, 'right': str.rjust, 'center': str.center}[align]
    df = pd.DataFrame({'a': ['xyzzy', 'zz']})
    head, row1, row2 = to_string(df, align=align, index=False).split('\n')
    assert head == just('a', 5) and row1 == 'xyzzy' and row2 == just('zz', 5)
    out = to_blocks(df.assign(grp='g'), 'grp', align=align, index=False)
    assert out.split('\n')[2] == head


@pytest.mark.parametrize('call', [lambda a: to_string(frame(), align=a),
                                  lambda a: to_blocks(frame(), 'grp', align=a),
                                  lambda a: to_blocks(frame(), 'grp', align=a, insert=[])])
def test_bad_align_names_accepted_values(call):
    with pytest.raises(ValueError, match="'left', 'right' or 'center'"):
        call('middle')


# 9 --------------------------------------------------------------------------
def test_old_style_header_hook_is_explained():
    def old(df, groupby, stats, **kw):
        return str(df.name)
    with pytest.raises(TypeError, match=r'header\(key, groupby, stats\)') as info:
        to_blocks(frame(), 'grp', header=old)
    assert isinstance(info.value.__cause__, AttributeError)


def test_default_header_errors_are_not_masked(monkeypatch):
    def boom(key, groupby, stats=None, **kw):
        raise AttributeError('inner')
    monkeypatch.setattr(rblocks, 'default_header', boom)
    with pytest.raises(AttributeError, match='inner'):
        to_blocks(frame(), 'grp', header=boom)


# 10 -------------------------------------------------------------------------
def test_apply_path_warns_once_at_the_caller(recwarn):
    from rotifer.pandas.blocks import default_block_rows
    to_blocks(frame(), 'grp', apply=default_block_rows, header=default_header,
              columns=['v'], stats={'v': 'max'})
    warns = [w for w in recwarn if issubclass(w.category, FutureWarning)]
    assert len(warns) == 1
    assert warns[0].filename == __file__


# 11 -------------------------------------------------------------------------
def test_row_rank_sorts_only_the_sort_columns(monkeypatch):
    seen = []
    orig = pd.DataFrame.sort_values

    def spy(self, *a, **k):
        seen.append(list(self.columns))
        return orig(self, *a, **k)
    monkeypatch.setattr(pd.DataFrame, 'sort_values', spy)
    rblocks._row_rank(frame(), {'v': False})
    assert seen == [['v']]


# 12 -------------------------------------------------------------------------
def test_to_color_does_not_join_all_values():
    from rotifer.pandas import text
    out = to_color(pd.Series(['AC', 'CA']), colors={'A': '1', 'C': '2'})
    assert out.iloc[0].count('\033[38;5;1m') == 1 and out.iloc[1].count('\033[38;5;2m') == 1
    assert "''.join(values)" not in open(text.__file__).read()


# 16 -------------------------------------------------------------------------
def test_to_blocks_and_to_string_share_the_helpers():
    from rotifer.pandas import text
    assert rblocks._text_columns is text._text_columns
    assert rblocks._join_rows is text._join_rows
    assert not hasattr(rblocks, '_string_columns') and not hasattr(rblocks, '_render_rows')


# 18 -------------------------------------------------------------------------
def test_empty_columns_list_agrees_between_writers():
    out = to_yaml(frame(), 'grp', columns=[])
    assert '| v' not in out and '| w' not in out
    assert 'v' not in to_blocks(frame(), 'grp', columns=[], index=False).replace('grp', '')


# 20 -------------------------------------------------------------------------
def test_pad_rejects_non_pandas_input():
    with pytest.raises(TypeError, match='DataFrame or Series'):
        pad(['a', 'b'])
    df = pd.DataFrame({'a': [1, 22]})
    pad(df, header=True)
    assert df['a'].tolist() == [1, 22]          # input untouched


# 23 -------------------------------------------------------------------------
@pytest.mark.parametrize('write', [to_blocks, to_yaml])
def test_bad_aggregation_names_table_and_column(write):
    kw = {'stats': {'v': 'bogus'}}
    with pytest.raises(ValueError, match=r"table 'stats': cannot aggregate column 'v' with 'bogus'"):
        write(frame(), 'grp', **kw)


# to_string signature --------------------------------------------------------
def test_to_string_rejects_unknown_keywords():
    with pytest.raises(TypeError):
        to_string(frame(), bogus=1)
