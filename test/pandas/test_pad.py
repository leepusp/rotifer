"""Tests for rotifer.pandas.text.pad with na_rep and header=True.

The golden values are the output of the legacy devel function padding_df on the
frame below, for each of its variants. The legacy function padded a column only
to its value width, so a column name longer than its values misaligned the header;
for that column the tests assert the corrected alignment instead."""
import numpy as np
import pandas as pd
import pytest

from rotifer.pandas import text

LONG = 'a_very_long_column_name'


def frame():
    return pd.DataFrame({
        'a': ['x', 'yyyy', None, 'zz'],
        LONG: [1, 22, 333, 4],
        'f': [1.5, np.nan, 22.25, 3.0],
        'mid': ['abc', 'de', 'f', np.nan],
    })


# Golden columns of the legacy output (name, values), for the columns whose name
# is no longer than their values.
GOLDEN = {
    'left': {
        '   a': ['   x', 'yyyy', '    ', '  zz'],
        '    f': ['  1.5', '     ', '22.25', '  3.0'],
        'mid': ['abc', ' de', '  f', '   '],
    },
    'right': {
        'a   ': ['x   ', 'yyyy', '    ', 'zz  '],
        'f    ': ['1.5  ', '     ', '22.25', '3.0  '],
        'mid': ['abc', 'de ', 'f  ', '   '],
    },
}

# variant of the legacy function -> the pad side that reproduces it
VARIANTS = {
    'blast': 'left',
    'rodolfo': 'left',
    'gian_func_left': 'left',
    'gian_func_right': 'right',
    'gian_func_default': 'right',
}


@pytest.mark.parametrize('variant', sorted(VARIANTS))
def test_header_matches_legacy_golden_for_short_names(variant):
    side = VARIANTS[variant]
    out = text.pad(frame(), side=side, na_rep='', header=True)
    for name, values in GOLDEN[side].items():
        assert name in out.columns
        assert out[name].tolist() == values


@pytest.mark.parametrize('side', ['left', 'right'])
def test_long_name_column_is_aligned_with_its_header(side):
    # Legacy bug: the values were padded to 3 while the 23-character name was not,
    # so header and values misaligned. Corrected: values are padded to the name width.
    out = text.pad(frame(), side=side, header=True)
    width = len(LONG)
    name = out.columns[1]
    assert name == (LONG.rjust(width) if side == 'left' else LONG.ljust(width))
    cells = out[name].tolist()
    assert all(len(c) == width for c in cells)
    expected = ['1', '22', '333', '4']
    if side == 'left':
        assert cells == [v.rjust(width) for v in expected]
    else:
        assert cells == [v.ljust(width) for v in expected]


def test_header_false_keeps_column_names_and_value_width():
    out = text.pad(frame(), side='left', header=False)
    assert list(out.columns) == ['a', LONG, 'f', 'mid']
    assert out[LONG].tolist() == ['  1', ' 22', '333', '  4']


def test_header_both_centres_name_like_values():
    df = pd.DataFrame({'ab': ['xxxxx', 'y'], 'q': [1, 22]})
    out = text.pad(df, side='both', header=True)
    assert list(out.columns) == ['  ab ', 'q ']
    assert out.to_dict('list') == {'  ab ': ['xxxxx', '  y  '], 'q ': ['1 ', '22']}


def test_na_rep_replaces_missing_before_padding():
    assert text.pad(pd.Series(['a', None, 'ccc']), na_rep='.').tolist() == ['  a', '  .', 'ccc']
    out = text.pad(pd.DataFrame({'a': ['x', None]}), na_rep='-', side='right')
    assert out['a'].tolist() == ['x', '-']


def test_na_rep_defaults_to_empty():
    assert text.pad(pd.Series(['a', None])).tolist() == ['a', ' ']
