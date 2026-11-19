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
import importlib.util
import sys
import warnings
from importlib.machinery import ModuleSpec
from unittest import mock

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


# ---- to_network: rotifer.pandas.to_network and epsoares.to_network ----------

def network_input():
    pytest.importorskip('networkx')
    from rotifer.genome.data import NeighborhoodDF
    rows = [
        (1, 'p1', 1, 1, 0, 'X'), (1, 'p2', 2, 1, 1, 'A+B'), (1, 'p3', 3, 1, 0, 'C'),
        (2, 'p4', 1, -1, 0, 'C'), (2, 'p5', 2, -1, 1, 'A'), (2, 'p6', 3, -1, 0, 'X'),
    ]
    df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'strand', 'query', 'pfam'])
    df['nucleotide'] = ['n1'] * 3 + ['n2'] * 3
    df['end'] = df.start + 1
    df['type'] = 'CDS'
    df['assembly'] = 'a'
    df['internal_id'] = range(6)
    df['feature_order'] = [1, 2, 3, 1, 2, 3]
    df['topology'] = 'linear'
    return NeighborhoodDF(df)


def test_pandas_to_network_warns_at_the_caller():
    import rotifer.pandas
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        shim = rotifer.pandas.to_network
    assert len(w) == 1 and issubclass(w[0].category, FutureWarning)
    assert w[0].filename == __file__
    message = str(w[0].message)
    assert 'rotifer.pandas.to_network' in message and REMOVAL in message
    assert 'NeighborhoodDF.to_network' in message
    df = network_input()
    pd.testing.assert_frame_equal(shim(df), df.to_network())
    pd.testing.assert_frame_equal(shim(df, strand=False), df.to_network(strand=False))


def test_importing_rotifer_pandas_does_not_warn():
    import subprocess
    import sys
    code = 'import rotifer.pandas; from rotifer.pandas import pad'
    result = subprocess.run([sys.executable, '-W', 'error::FutureWarning', '-c', code],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_pandas_unknown_attribute_still_raises():
    import rotifer.pandas
    with pytest.raises(AttributeError):
        rotifer.pandas.no_such_function


# epsoares imports heavy optional third-party packages at module level (and so do
# the rotifer.devel modules it imports). The to_network shim uses none of them, so
# the fixture below stubs only the names that importlib cannot find in this venv.
# Installed packages are never stubbed, and rotifer.* and pandas/numpy are never
# stubbed. The stubs are mocks because the import itself reads attributes from them
# (for example `from joblib import Parallel, delayed`).
EPSOARES_THIRD_PARTY = (
    'pyhmmer', 'pygraphviz', 'joblib', 'ete3', 'seaborn', 'matplotlib', 'tqdm',
    'psutil', 'Bio', 'IPython',
)


def make_stub(name):
    stub = mock.MagicMock(name=name)
    stub.__spec__ = ModuleSpec(name, None)
    return stub


@pytest.fixture
def fresh_epsoares(monkeypatch):
    """Import rotifer.devel.alpha.epsoares fresh, with missing optional deps stubbed.

    Cached rotifer.devel modules are dropped first so the import re-runs, and
    everything the import adds is removed afterwards; monkeypatch restores the
    rest of sys.modules.
    """
    stubbed = [name for name in EPSOARES_THIRD_PARTY
               if importlib.util.find_spec(name) is None]
    before = set(sys.modules)
    for name in list(sys.modules):
        if name == 'rotifer.devel' or name.startswith('rotifer.devel.'):
            monkeypatch.delitem(sys.modules, name)
    for name in stubbed:
        monkeypatch.setitem(sys.modules, name, make_stub(name))
    try:
        yield importlib.import_module('rotifer.devel.alpha.epsoares')
    finally:
        for name in set(sys.modules) - before:
            if name == 'rotifer.devel' or name.startswith('rotifer.devel.'):
                sys.modules.pop(name, None)


def test_epsoares_to_network_warns_and_matches_the_method(fresh_epsoares):
    to_network = fresh_epsoares.to_network
    df = network_input()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        got = to_network(df)
    deprecations = [x for x in w if issubclass(x.category, FutureWarning)]
    assert len(deprecations) == 1, [str(x.message) for x in w]
    assert deprecations[0].filename == __file__
    message = str(deprecations[0].message)
    assert 'NeighborhoodDF.to_network' in message and REMOVAL in message
    pd.testing.assert_frame_equal(got, df.to_network())
