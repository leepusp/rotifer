#!/usr/bin/env python3
"""Golden-output tests for rotifer.pandas.to_string.

The expected strings were captured from the pre-pyarrow implementation, so they
pin the exact rendering across the pyarrow and pandas backends alike. Both are
exercised: _to_string_arrow is the default, _to_string_pandas the fallback used
when pyarrow is unavailable.

Run under pytest, or standalone: python test/pandas/test_to_string.py
"""

import numpy as np
import pandas as pd
import pytest

import rotifer.pandas as rp

# name -> (DataFrame, to_string kwargs)
def cases():
    df = pd.DataFrame({'a': ['x', 'yy', 'zzz'], 'b': [1, 22, 333], 'f': [1.5, np.nan, 3.25]})
    strindex = df.copy()
    strindex.index = ['i1', 'i22', 'i333']
    strindex.index.name = 'idx'
    return {
        'basic':      (df, dict()),
        'right':      (df, dict(align='right')),
        'noindex':    (df, dict(index=False)),
        'ff':         (df, dict(float_format='%.3f')),
        'ff_narep':   (df, dict(float_format='%.2f', na_rep='NULL')),
        'sep':        (df, dict(sep='|')),
        'insert':     (df, dict(insert=[(0, 'new', ['p', 'q', 'r'])])),
        'longheader': (pd.DataFrame({'averyverylongcolumnname': ['a', 'b']}), dict()),
        'strindex':   (strindex, dict()),
        'nastr':      (pd.DataFrame({'a': ['x', None, 'zzz']}), dict()),
        'nastr_rep':  (pd.DataFrame({'a': ['x', None, 'zzz']}), dict(na_rep='-')),
        'empty':      (pd.DataFrame({'a': [], 'b': []}), dict()),
        'unicode':    (pd.DataFrame({'a': ['héllo', 'ok'], 'b': ['ß', 'xx']}), dict()),
        'onecol':     (pd.DataFrame({'a': ['x', 'yy']}), dict(index=False)),
    }

GOLDEN = {
    'basic':      '  : a   : b   : f   \n0 : x   : 1   : 1.5 \n1 : yy  : 22  :     \n2 : zzz : 333 : 3.25',
    'right':      '  :   a :   b :    f\n0 :   x :   1 :  1.5\n1 :  yy :  22 :     \n2 : zzz : 333 : 3.25',
    'noindex':    'a   : b   : f   \nx   : 1   : 1.5 \nyy  : 22  :     \nzzz : 333 : 3.25',
    'ff':         '  : a   : b   : f    \n0 : x   : 1   : 1.500\n1 : yy  : 22  :      \n2 : zzz : 333 : 3.250',
    'ff_narep':   '  : a   : b   : f   \n0 : x   : 1   : 1.50\n1 : yy  : 22  : NULL\n2 : zzz : 333 : 3.25',
    'sep':        ' |a  |b  |f   \n0|x  |1  |1.5 \n1|yy |22 |    \n2|zzz|333|3.25',
    'insert':     '  : new : a   : b   : f   \n0 : p   : x   : 1   : 1.5 \n1 : q   : yy  : 22  :     \n2 : r   : zzz : 333 : 3.25',
    'longheader': '  : averyverylongcolumnname\n0 : a                      \n1 : b                      ',
    'strindex':   'idx  : a   : b   : f   \ni1   : x   : 1   : 1.5 \ni22  : yy  : 22  :     \ni333 : zzz : 333 : 3.25',
    'nastr':      '  : a  \n0 : x  \n1 :    \n2 : zzz',
    'nastr_rep':  '  : a  \n0 : x  \n1 : -  \n2 : zzz',
    'empty':      ' : a : b\n',
    'unicode':    '  : a     : b \n0 : héllo : ß \n1 : ok    : xx',
    'onecol':     'a \nx \nyy',
}


@pytest.fixture(params=['arrow', 'pandas'])
def backend(request, monkeypatch):
    """to_string dispatches to _to_string_arrow, falling back to _to_string_pandas
    on ImportError. Point the former at the latter to exercise the fallback."""
    if request.param == 'pandas':
        monkeypatch.setattr(rp, '_to_string_arrow', rp._to_string_pandas)
    return request.param


@pytest.mark.parametrize('name', sorted(GOLDEN))
def test_to_string_matches_golden(name, backend):
    df, kwargs = cases()[name]
    assert rp.to_string(df, **kwargs) == GOLDEN[name]


def test_both_backends_agree():
    """The fallback must be byte-identical to the pyarrow path, not merely close."""
    for name, (df, kwargs) in cases().items():
        just = str.rjust if kwargs.get('align') == 'right' else str.ljust
        arrow = rp.to_string(df, **kwargs)
        original = rp._to_string_arrow
        try:
            rp._to_string_arrow = rp._to_string_pandas
            fallback = rp.to_string(df, **kwargs)
        finally:
            rp._to_string_arrow = original
        assert arrow == fallback, name


if __name__ == '__main__':
    import sys
    failures = 0
    for impl in ('arrow', 'pandas'):
        original = rp._to_string_arrow
        if impl == 'pandas':
            rp._to_string_arrow = rp._to_string_pandas
        try:
            bad = [n for n, (d, k) in cases().items() if rp.to_string(d, **k) != GOLDEN[n]]
        finally:
            rp._to_string_arrow = original
        failures += len(bad)
        print('%-7s %s' % (impl, 'ALL %d MATCH' % len(GOLDEN) if not bad else 'MISMATCH %s' % bad))
    sys.exit(1 if failures else 0)
