#!/usr/bin/env python3
"""Tests for the deprecated devel copies of padding_df (issue #TBD).

Three devel modules each held a copy of padding_df. Each is now a thin shim over
rotifer.pandas.pad(df, side=..., na_rep='', header=True), which reproduces the
legacy output. What has to hold, and why:

  * calling a shim emits a FutureWarning that names rotifer.pandas.pad and the
    removal date 2026-11-30, and points at the CALLER's line (stacklevel=2), not
    at the shim: a user needs to find the call site that must change
  * the output of a shim equals rotifer.pandas.pad with the matching side, so
    the shim is a pure alias and cannot drift from the new function
  * gian_func.padding_df takes how='right' (default) or how='left' and maps it
    to side=how

The frame has a NaN, an int column and column names shorter than their values.
That is the case where the legacy copies and pad agree. A name longer than all
its values is a known, intended difference (pad pads the values to the header
width), so it is not covered here.

Modules that cannot be imported in this environment are skipped, not failed.

Run under pytest: PYTHONPATH=lib pytest test/devel/test_pandas_shims.py
"""

import importlib
import warnings

import numpy as np
import pandas as pd
import pytest

from rotifer.pandas import pad

REMOVAL = '2026-11-30'

# (module, function, kwargs passed to the shim, equivalent pad call kwargs)
SHIMS = [
    ('rotifer.devel.beta.blast', 'padding_df', {}, {'side': 'left'}),
    ('rotifer.devel.alpha.rodolfo', 'padding_df', {}, {'side': 'left'}),
    ('rotifer.devel.alpha.gian_func', 'padding_df', {'how': 'right'}, {'side': 'right'}),
    ('rotifer.devel.alpha.gian_func', 'padding_df', {'how': 'left'}, {'side': 'left'}),
    # Placeholder: add the epsoares.to_network shim entry here when it lands
]


def frame():
    return pd.DataFrame({
        'a': ['x', 'yyyy', None, 'zz'],
        'n': [1, 22, 333, 4],
        'f': [1.5, np.nan, 22.25, 3.0],
        'mid': ['abc', 'de', 'f', np.nan],
    })


def load(module, function):
    try:
        mod = importlib.import_module(module)
    except ImportError as exc:  # missing optional dependency in this venv
        pytest.skip(f'{module} not importable: {exc}')
    return getattr(mod, function)


def shim_id(param):
    module, function, kwargs, _ = param
    suffix = ','.join(f'{k}={v}' for k, v in kwargs.items())
    return f'{module.rsplit(".", 1)[-1]}.{function}({suffix})'


@pytest.mark.parametrize('param', SHIMS, ids=shim_id)
def test_call_warns_at_the_caller(param):
    module, function, kwargs, _ = param
    shim = load(module, function)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        shim(frame(), **kwargs)
    deprecations = [x for x in w if issubclass(x.category, FutureWarning)]
    assert len(deprecations) == 1, [str(x.message) for x in w]
    warning = deprecations[0]
    assert warning.filename == __file__
    message = str(warning.message)
    assert 'rotifer.pandas.pad' in message
    assert REMOVAL in message
    assert f'{module}.{function}' in message


@pytest.mark.parametrize('param', SHIMS, ids=shim_id)
def test_output_equals_pad(param):
    module, function, kwargs, pad_kwargs = param
    shim = load(module, function)
    df = frame()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        got = shim(df, **kwargs)
    expected = pad(df, na_rep='', header=True, **pad_kwargs)
    pd.testing.assert_frame_equal(got, expected)
    # the input frame is not modified
    pd.testing.assert_frame_equal(df, frame())
