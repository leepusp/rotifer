#!/usr/bin/env python3
"""Golden-output tests for rotifer.pandas.to_string.

The expected strings were captured from the pre-pyarrow implementation, so they
pin the exact rendering.

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


@pytest.mark.parametrize('name', sorted(GOLDEN))
def test_to_string_matches_golden(name):
    df, kwargs = cases()[name]
    assert rp.to_string(df, **kwargs) == GOLDEN[name]


if __name__ == '__main__':
    import sys
    bad = [n for n, (d, k) in cases().items() if rp.to_string(d, **k) != GOLDEN[n]]
    print('ALL %d MATCH' % len(GOLDEN) if not bad else 'MISMATCH %s' % bad)
    sys.exit(1 if bad else 0)
