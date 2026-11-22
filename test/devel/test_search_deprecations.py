#!/usr/bin/env python3
"""Tests for the deprecation of the search helpers in rotifer.devel.alpha (issue #141).

The helpers below were superseded by rotifer.devel.beta.hmmer and
rotifer.devel.beta.blast and will be removed on November 30, 2026. Until then
each one has to tell its callers so, and only its callers, because the
warning is the only signal they get.

What has to hold, and why:

  * each public deprecated function emits a FutureWarning whose message names
    the removal date and the replacement, so the warning says what to switch to
  * the epsoares wrappers pass their arguments and return value through to the
    implementation unchanged, so existing callers keep their results
  * calling the implementation from inside epsoares (hmmscan, add_arch_to_df,
    igem_pipeline) does not warn, so one deprecated call warns once, and the
    internal callers are not reported as users of the old names
  * the epsoares implementations live under a private name (_hmmscan and
    friends), so nothing in the module reaches the public deprecated wrapper
"""
import subprocess
import warnings

import pytest

DATE = 'November 30, 2026'


def _warnings_from(fn, *args, **kwargs):
    """Call fn, swallowing any exception it raises, and return the FutureWarnings it emitted."""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        try:
            fn(*args, **kwargs)
        except Exception:
            pass
    return [w for w in rec if issubclass(w.category, FutureWarning)]


@pytest.fixture
def epsoares():
    return pytest.importorskip('rotifer.devel.alpha.epsoares')


@pytest.fixture
def rodolfo():
    return pytest.importorskip('rotifer.devel.alpha.rodolfo')


@pytest.fixture
def gian_func():
    return pytest.importorskip('rotifer.devel.alpha.gian_func')


# --- epsoares: wrappers around private implementations ------------------------

EPSOARES_FUNCS = ['pyhmmer_to_df', 'hmmscan_linear', 'hmmsearch', 'hmmscan']


@pytest.mark.parametrize('name', EPSOARES_FUNCS)
def test_epsoares_public_function_warns_with_date_and_replacement(epsoares, monkeypatch, name):
    monkeypatch.setattr(epsoares, f'_{name}', lambda *a, **k: 'ok')
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        getattr(epsoares, name)('arg')
    msgs = [str(w.message) for w in rec if issubclass(w.category, FutureWarning)]
    assert len(msgs) == 1
    assert DATE in msgs[0]
    assert f'rotifer.devel.beta.hmmer.{name}' in msgs[0]


@pytest.mark.parametrize('name', EPSOARES_FUNCS)
def test_epsoares_wrapper_passes_arguments_and_result_through(epsoares, monkeypatch, name):
    seen = {}

    def impl(*args, **kwargs):
        seen['args'], seen['kwargs'] = args, kwargs
        return 'result'

    monkeypatch.setattr(epsoares, f'_{name}', impl)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        assert getattr(epsoares, name)(1, 2, columns=['a'], rename=False) == 'result'
    assert seen == {'args': (1, 2), 'kwargs': {'columns': ['a'], 'rename': False}}


def test_epsoares_implementation_has_private_name(epsoares):
    for name in EPSOARES_FUNCS:
        assert callable(getattr(epsoares, f'_{name}'))


def test_epsoares_internal_call_does_not_warn(epsoares):
    # With no models, _hmmsearch takes its fallback path through the private
    # _pyhmmer_to_df and returns the empty table; that internal call must not warn.
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        out = epsoares._hmmsearch(models_path=[], query_db='unused.fasta')
    assert len(out) == 0
    assert not [w for w in rec if issubclass(w.category, FutureWarning)]


def test_epsoares_add_arch_to_df_calls_private_hmmscan(epsoares, monkeypatch):
    # add_arch_to_df runs hmmscan internally: the deprecated public name must
    # not be the one it reaches, or every caller of add_arch_to_df would warn.
    import pandas as pd
    calls = []

    def fake_hmmscan(*a, **k):
        calls.append(True)
        raise RuntimeError('stop after the call')

    monkeypatch.setattr(epsoares, '_hmmscan', fake_hmmscan)
    monkeypatch.setattr(epsoares, 'hmmscan', lambda *a, **k: pytest.fail('public hmmscan was called'))
    df = pd.DataFrame({'pid': ['p1']})
    with pytest.raises(RuntimeError, match='stop after the call'):
        epsoares.add_arch_to_df(df, run_hmmscan=True)
    assert calls == [True]


# --- rodolfo and gian_func: warning at the top of the body --------------------

PLAIN = [
    ('rodolfo', 'psiblast', 'blast'),
    ('rodolfo', 'search2aln', 'blast'),
    ('rodolfo', 'split_by_model', 'blast'),
    ('gian_func', 'psiblast', 'blast'),
    ('gian_func', 'search2aln', 'blast'),
]


@pytest.mark.parametrize('module,name,newmodule', PLAIN)
def test_plain_public_function_warns_with_date_and_replacement(request, monkeypatch, module, name, newmodule):
    mod = request.getfixturevalue(module)
    # The warning is the first statement of the body: the call only has to get
    # that far. Stop every external program before it can run.
    def no_run(*a, **k):
        raise RuntimeError('external program not run in tests')
    monkeypatch.setattr(subprocess, 'Popen', no_run)
    monkeypatch.setattr(subprocess, 'run', no_run)
    monkeypatch.setattr(subprocess, 'call', no_run)
    recorded = _warnings_from(getattr(mod, name), None, None)
    assert len(recorded) == 1, [str(w.message) for w in recorded]
    msg = str(recorded[0].message)
    assert DATE in msg
    assert f'rotifer.devel.alpha.{module}.{name} is deprecated' in msg
    assert f'rotifer.devel.beta.{newmodule}.{name}' in msg
