"""Tests for rotifer.pandas.display.optimize_memory_usage."""
import warnings

import numpy as np
import pandas as pd
import pytest

from rotifer.pandas import display


def test_empty_frame_is_returned_unchanged_copy():
    # Bug: the text ratio divided by the row count, so an empty frame raised ZeroDivisionError.
    df = pd.DataFrame({'i': pd.Series([], dtype='int64'), 's': pd.Series([], dtype=object)})
    out = display.optimize_memory_usage(df)
    assert out is not df
    pd.testing.assert_frame_equal(out, df)
    assert out.dtypes.astype(str).to_dict() == {'i': 'int64', 's': 'object'}


def test_nullable_int_stays_nullable_and_keeps_na():
    df = pd.DataFrame({'n': pd.array([1, None, 300], dtype='Int64')})
    out = display.optimize_memory_usage(df)
    assert str(out['n'].dtype) == 'Int16'
    assert out['n'].isna().tolist() == [False, True, False]
    assert out['n'].tolist()[0] == 1 and out['n'].tolist()[2] == 300


def test_nullable_float_is_downcast():
    df = pd.DataFrame({'f': pd.array([1.5, None, 2.0], dtype='Float64')})
    out = display.optimize_memory_usage(df)
    assert str(out['f'].dtype) == 'Float32'
    assert out['f'].isna().tolist() == [False, True, False]


def test_text_becomes_category_only_when_few_distinct_values():
    df = pd.DataFrame({'lo': ['a', 'a', 'b', 'b'], 'hi': ['w', 'x', 'y', 'z'],
                       'o_lo': pd.Series(['a', 'a', 'b', 'b'], dtype=object)})
    out = display.optimize_memory_usage(df)
    assert out.dtypes.astype(str).to_dict() == {'lo': 'category', 'hi': 'str', 'o_lo': 'category'}


def test_other_dtypes_are_left_alone():
    df = pd.DataFrame({'b': [True, False], 'd': pd.to_datetime(['2020-01-01', '2021-01-01'])})
    out = display.optimize_memory_usage(df)
    assert out.dtypes.equals(df.dtypes)


def test_no_warning_is_emitted_for_any_dtype():
    df = pd.DataFrame({
        'i': [1, 2, 3, 4], 'f': [1.5, np.nan, 2.0, 4.0], 'lo': ['a', 'a', 'b', 'b'],
        'hi': ['w', 'x', 'y', 'z'], 'o_hi': pd.Series(['w', 'x', 'y', 'z'], dtype=object),
        'n': pd.array([1, None, 2, 3], dtype='Int64'), 'b': [True, False, True, True]})
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        display.optimize_memory_usage(df)


def test_input_is_not_mutated():
    df = pd.DataFrame({'i': [1, 2, 3, 70000], 'f': [1.5, 2.5, np.nan, 4.0],
                       'lo': ['a', 'a', 'b', 'b'], 'n': pd.array([1, None, 2, 3], dtype='Int64')})
    before = df.copy()
    display.optimize_memory_usage(df)
    pd.testing.assert_frame_equal(df, before)
    assert df.dtypes.astype(str).to_dict() == {'i': 'int64', 'f': 'float64', 'lo': 'str', 'n': 'Int64'}


def test_duplicate_column_labels_are_converted_by_position():
    df = pd.DataFrame([[1, 'a'], [2, 'b']], columns=[0, 0])
    out = display.optimize_memory_usage(df)
    assert out.shape == (2, 2)
    assert str(out.dtypes.iloc[0]) == 'int8'
