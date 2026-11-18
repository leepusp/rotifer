"""Tests for rotifer.core.session.codecs. Run with pytest or as a script."""
import os
import subprocess
import sys
import textwrap

import pandas as pd
import pytest

from rotifer.core.session import codecs
from rotifer.core.session.codecs import (Codec, decode, encode, get_codec, lookup,
                                    register_alias, register_codec)


class Point:
    def __init__(self, x, y):
        self.x, self.y = x, y

    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y)


class Base:
    pass


class Child(Base):
    pass


class SubFrame(pd.DataFrame):
    @property
    def _constructor(self):
        return SubFrame


def roundtrip(obj, tmp_path, **kw):
    codec, meta = encode(obj, tmp_path / 'var', **kw)
    return codec, meta, decode(codec.name, tmp_path / 'var', meta)


def created(tmp_path):
    return sorted(p.name for p in tmp_path.iterdir() if not p.name.startswith('.'))


# ---------------------------------------------------------------- pickle

@pytest.mark.parametrize('obj', [{'a': [1, 2, (3, 4)]}, [1, 'x', None], Point(1, 2), 'text'])
def test_pickle_roundtrip(tmp_path, obj):
    codec, meta, back = roundtrip(obj, tmp_path)
    assert codec.name == 'pickle'
    assert back == obj
    assert meta['version_bound'] is False
    assert meta['files'] == ['var.pkl'] == created(tmp_path)
    assert 'fallback_from' not in meta


def test_gzip(tmp_path):
    obj = list(range(1000))
    codec, meta, back = roundtrip(obj, tmp_path, compress='gzip')
    assert back == obj
    assert meta['compress'] == 'gzip'
    assert meta['files'] == ['var.pkl.gz'] == created(tmp_path)
    with open(tmp_path / 'var.pkl.gz', 'rb') as f:
        assert f.read(2) == b'\x1f\x8b'


def test_bad_compress_value(tmp_path):
    with pytest.raises(ValueError):
        encode([1], tmp_path / 'var', compress='zip')


# --------------------------------------------------------------- parquet

def test_parquet_plain_frame(tmp_path):
    df = pd.DataFrame({
        'i': [1, 2, 3],
        'f': [0.5, None, 2.0],
        's': ['a', 'b', 'c'],
        'cat': pd.Categorical(['u', 'v', 'u']),
        'n': pd.array([1, None, 3], dtype='Int64'),
    }, index=pd.Index([10, 20, 30], name='id'))
    codec, meta, back = roundtrip(df, tmp_path)
    assert codec.name == 'parquet-dataframe'
    assert meta['files'] == ['var.parquet'] == created(tmp_path)
    assert df.equals(back)
    assert list(df.dtypes) == list(back.dtypes)
    assert back.index.name == 'id'
    assert str(back['n'].dtype) == 'Int64'
    assert isinstance(back['cat'].dtype, pd.CategoricalDtype)


def test_parquet_fingerprint(tmp_path):
    codec = get_codec('parquet-dataframe')
    df = pd.DataFrame({'a': [1, 2], 'b': ['x', 'y']})
    assert codec.fingerprint(df).startswith('pandas:')


def test_mixed_object_column_falls_back(tmp_path):
    df = pd.DataFrame({'m': [1, 'a', [2, 3]]})
    codec, meta, back = roundtrip(df, tmp_path)
    assert codec.name == 'pickle'
    assert meta['fallback_from'] == 'parquet-dataframe'
    assert meta['fallback_reason']
    assert back['m'].tolist() == [1, 'a', [2, 3]]
    assert meta['files'] == ['var.pkl'] == created(tmp_path)   # no parquet leftover


def test_roundtrip_mismatch_falls_back(tmp_path):
    df = pd.DataFrame({'a': [1, 2]})
    orig = codecs._read_parquet

    # make the stored file differ from the input by lying in pandas.read_parquet
    real = pd.read_parquet
    pd.read_parquet = lambda *a, **k: real(*a, **k).assign(a=[9, 9])
    try:
        codec, meta = encode(df, tmp_path / 'var')
    finally:
        pd.read_parquet = real
    assert codecs._read_parquet is orig
    assert codec.name == 'pickle'
    assert meta['fallback_from'] == 'parquet-dataframe'
    assert meta['fallback_reason'] == 'round-trip mismatch'
    assert created(tmp_path) == ['var.pkl']


def test_dataframe_subclass_uses_pickle(tmp_path):
    df = SubFrame({'a': [1, 2]})
    assert lookup(df).name == 'pickle'
    codec, meta, back = roundtrip(df, tmp_path)
    assert codec.name == 'pickle' and 'fallback_from' not in meta
    assert type(back) is SubFrame and back.equals(df)


# ---------------------------------------------------------- registration

def _dummy(name='dummy', types=(complex,), **kw):
    return Codec(name, 1, types, lambda o, p: {'files': []}, lambda p, m: None, **kw)


@pytest.fixture
def clean_registry():
    names, by_type = dict(codecs._codecs), dict(codecs._by_type)
    aliases = dict(codecs._aliases)
    yield
    codecs._codecs.clear(); codecs._codecs.update(names)
    codecs._by_type.clear(); codecs._by_type.update(by_type)
    codecs._aliases.clear(); codecs._aliases.update(aliases)


def test_register_and_lookup(clean_registry):
    c = _dummy()
    register_codec(c)
    assert get_codec('dummy') is c
    assert lookup(1 + 2j) is c
    with pytest.raises(ValueError):
        register_codec(c)
    with pytest.raises(ValueError):
        register_codec(_dummy('other'))            # same type, other name
    c2 = _dummy()
    register_codec(c2, replace=True)
    assert get_codec('dummy') is c2


def test_get_codec_unknown():
    with pytest.raises(KeyError):
        get_codec('nope')


def test_subclass_lookup_rules(clean_registry):
    register_codec(_dummy('exact', (Base,)))
    assert lookup(Base()).name == 'exact'
    assert lookup(Child()).name == 'pickle'        # no opt-in
    register_codec(_dummy('subs', (Base,), subclasses=True), replace=True)
    assert lookup(Child()).name == 'subs'


def test_lookup_never_fails():
    class Weird:
        def __reduce__(self):
            raise RuntimeError('no')
    assert lookup(Weird()).name == 'cloudpickle'
    assert lookup(3).name == 'pickle'


def test_typed_codec_failure_falls_back(tmp_path, clean_registry):
    def boom(o, p):
        raise OSError('disk')
    register_codec(Codec('boomer', 1, (complex,), boom, lambda p, m: None))
    codec, meta = encode(1j, tmp_path / 'var')
    assert codec.name == 'pickle'
    assert meta['fallback_from'] == 'boomer'
    assert 'OSError' in meta['fallback_reason']
    assert decode(codec.name, tmp_path / 'var', meta) == 1j


# --------------------------------------------------------------- aliases

def test_alias_loads_old_class_path(tmp_path, clean_registry):
    import types
    mod = types.ModuleType('oldmod_for_alias')
    sys.modules['oldmod_for_alias'] = mod
    try:
        class Moved:
            def __init__(self): self.v = 7
        Moved.__module__, Moved.__qualname__ = 'oldmod_for_alias', 'Moved'
        mod.Moved = Moved
        codec, meta = encode(Moved(), tmp_path / 'var')
    finally:
        del sys.modules['oldmod_for_alias']
    with pytest.raises(ModuleNotFoundError):
        decode('pickle', tmp_path / 'var', meta)
    register_alias('oldmod_for_alias.Moved', __name__ + '.Point')
    obj = decode('pickle', tmp_path / 'var', meta)
    assert type(obj) is Point and obj.v == 7


# -------------------------------------------------------------- __main__

def test_uses_main_scanner():
    import pickle
    assert not codecs._uses_main(pickle.dumps({'a': Point(1, 2)}, 5))
    assert not codecs._uses_main(pickle.dumps(pickle.dumps, 5))
    # strings that merely say "__main__" are not references
    assert not codecs._uses_main(pickle.dumps(['__main__', 'x'], 5))
    fake = type('F', (), {})
    fake.__module__ = '__main__'
    sys.modules['__main__'].__dict__['_F_for_scan'] = fake
    fake.__qualname__ = '_F_for_scan'
    try:
        for proto in (2, 4, 5):
            assert codecs._uses_main(pickle.dumps([fake(), 1, fake], proto))
    finally:
        del sys.modules['__main__'].__dict__['_F_for_scan']


SAVE_SCRIPT = textwrap.dedent('''
    import sys
    from pathlib import Path
    from rotifer.core.session.codecs import encode
    import json

    def double(x):
        return 2 * x

    class Thing:
        def __init__(self, v):
            self.v = v
        def triple(self):
            return 3 * self.v

    out = Path(sys.argv[1])
    result = {}
    for label, obj in [("func", double), ("cls", Thing), ("inst", Thing(5)),
                       ("lam", lambda x: x + 1), ("nested", {"k": [Thing(2)]})]:
        codec, meta = encode(obj, out / label)
        result[label] = [codec.name, meta]
    json.dump(result, open(out / "meta.json", "w"))
''')

LOAD_SCRIPT = textwrap.dedent('''
    import sys, json
    from pathlib import Path
    from rotifer.core.session.codecs import decode
    out = Path(sys.argv[1])
    meta = json.load(open(out / "meta.json"))
    load = lambda k: decode(meta[k][0], out / k, meta[k][1])
    assert load("func")(4) == 8
    assert load("cls")(3).triple() == 9
    assert load("inst").triple() == 15
    assert load("lam")(1) == 2
    assert load("nested")["k"][0].triple() == 6
    print("LOADED")
''')


def _run(script, tmp_path, name):
    path = tmp_path / name
    path.write_text(script)
    env = dict(os.environ)
    import rotifer
    lib = os.path.dirname(os.path.dirname(os.path.abspath(rotifer.__file__)))
    env['PYTHONPATH'] = os.pathsep.join([lib] + env.get('PYTHONPATH', '').split(os.pathsep))
    return subprocess.run([sys.executable, str(path), str(tmp_path / 'out')],
                          env=env, capture_output=True, text=True, cwd=str(tmp_path))


def test_main_defined_values_use_cloudpickle(tmp_path):
    pytest.importorskip('cloudpickle')
    (tmp_path / 'out').mkdir()
    r = _run(SAVE_SCRIPT, tmp_path, 'save_script.py')
    assert r.returncode == 0, r.stderr
    import json
    meta = json.load(open(tmp_path / 'out' / 'meta.json'))
    for label, (codec, m) in meta.items():
        assert codec == 'cloudpickle', label
        assert m['version_bound'] is True
        assert m['files'] == [label + '.cloudpickle.pkl']
    # a fresh process, where nothing named double/Thing exists
    r = _run(LOAD_SCRIPT, tmp_path, 'load_script.py')
    assert r.returncode == 0, r.stderr
    assert 'LOADED' in r.stdout


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q']))


def test_lazy_codec_provider_and_registration(monkeypatch):
    import types
    mod = types.ModuleType("fake_provider_mod")
    monkeypatch.setitem(sys.modules, "fake_provider_mod", mod)
    codecs.register_codec_provider("fake.codec", "fake_provider_mod")

    def on_import():
        codecs.register_codec(Codec("fake.codec", 1, (), lambda o, p: {"files": []},
                                    lambda p, m: 1))
    # the provider module registers the codec when imported
    import importlib
    real = importlib.import_module
    monkeypatch.setattr(importlib, "import_module",
                        lambda n, *a, **k: (on_import(), mod)[1] if n == "fake_provider_mod"
                        else real(n, *a, **k))
    try:
        assert get_codec("fake.codec").name == "fake.codec"
    finally:
        codecs._codecs.pop("fake.codec", None)
    with pytest.raises(KeyError):
        get_codec("never.heard.of")


def test_parquet_datetimeindex_freq_roundtrips_exactly(tmp_path):
    df = pd.DataFrame({'a': [1, 2, 3]}, index=pd.date_range('2020', periods=3, freq='D'))
    codec, meta, back = roundtrip(df, tmp_path)
    assert back.index.freq == df.index.freq
    assert codec.name == 'parquet-dataframe' or 'fallback_from' in meta


_CALLS = []


class _Counting:
    def __reduce__(self):
        _CALLS.append(1)
        return (_Counting, ())


def test_pickle_encode_returns_fingerprint_without_second_dump(tmp_path):
    _CALLS.clear()
    codec, meta = encode([_Counting()], tmp_path / "v")
    assert codec.name == "pickle"
    assert meta["fingerprint"].startswith("pickle:")
    assert len(_CALLS) == 1
