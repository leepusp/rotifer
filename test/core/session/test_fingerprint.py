"""Tests for rotifer.core.session._fingerprint. Run with pytest or as a script."""
import sys

import numpy as np
import pandas as pd

from rotifer.core.session._fingerprint import fingerprint
from rotifer.core.session.codecs import Codec


def test_equal_frames_same_fingerprint():
    a = pd.DataFrame({'x': [1, 2, 3], 'y': list('abc')})
    b = a.copy()
    fa = fingerprint(a)
    assert fa == fingerprint(b)
    assert fa.startswith('pandas:')


def test_changed_value_changes_fingerprint():
    a = pd.DataFrame({'x': [1, 2, 3]})
    b = a.copy()
    b.loc[1, 'x'] = 99
    assert fingerprint(a) != fingerprint(b)


def test_column_rename_and_dtype_change_detected():
    a = pd.DataFrame({'x': [1, 2, 3]})
    assert fingerprint(a) != fingerprint(a.rename(columns={'x': 'z'}))
    assert fingerprint(a) != fingerprint(a.astype('float64'))


def test_index_matters():
    a = pd.DataFrame({'x': [1, 2]})
    assert fingerprint(a) != fingerprint(a.set_axis([5, 6]))


def test_series():
    s = pd.Series([1, 2, 3], name='n')
    assert fingerprint(s) == fingerprint(s.copy())
    assert fingerprint(s) != fingerprint(s.rename('m'))


def test_numpy():
    a = np.arange(6).reshape(2, 3)
    assert fingerprint(a).startswith('numpy:')
    assert fingerprint(a) == fingerprint(a.copy())
    assert fingerprint(a) == fingerprint(np.asfortranarray(a))
    assert fingerprint(a) != fingerprint(a.reshape(3, 2))
    assert fingerprint(a) != fingerprint(a.astype('int32'))


def test_other_objects_use_pickle():
    fp = fingerprint({'a': [1, 2]})
    assert fp.startswith('pickle:')
    assert fp == fingerprint({'a': [1, 2]})
    assert fp != fingerprint({'a': [1, 3]})


def test_unhashable_returns_none():
    assert fingerprint(lambda: 1) is None


def test_codec_fingerprint_wins():
    codec = Codec('x', 1, (), None, None, fingerprint=lambda o: 'custom:1')
    assert fingerprint([1], codec) == 'custom:1'
    plain = Codec('y', 1, (), None, None)
    assert fingerprint([1], plain).startswith('pickle:')


if __name__ == '__main__':
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn()
                print('  ok  ', name)
            except Exception as exc:
                failures += 1
                print('  FAIL', name, repr(exc))
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
