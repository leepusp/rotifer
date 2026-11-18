#!/usr/bin/env python3
"""Tests for rotifer.core.session._legacy (read-only reader of old shelve sessions).

Run under pytest, or standalone: python test/core/session/test_legacy.py
"""

import dbm
import dbm.dumb
import pickle
import shelve
import sys

import pandas as pd
import pytest

from rotifer.core.session._legacy import is_legacy_session, legacy_target, read_legacy

try:
    import dbm.gnu
except ImportError:
    dbm.gnu = None

BACKENDS = ["dumb", "gnu"]


def _open_raw(backend, p, flag="c"):
    if backend == "dumb":
        return dbm.dumb.open(str(p), flag)
    if dbm.gnu is None:
        pytest.skip("dbm.gnu not available")
    return dbm.gnu.open(str(p), flag)


def _values():
    return {"n": 3, "l": [1, "a", 2.5], "d": {"k": [1, 2]},
            "df": pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})}


@pytest.fixture(params=BACKENDS)
def session_file(request, tmp_path):
    p = tmp_path / "old"
    raw = _open_raw(request.param, p)
    shelf = shelve.Shelf(raw)
    for k, v in _values().items():
        shelf[k] = v
    shelf.close()
    return p


def _snapshot(p):
    return {f.name: (f.read_bytes(), f.stat().st_mtime_ns)
            for f in p.parent.iterdir()}


def test_detection(session_file, tmp_path):
    assert is_legacy_session(session_file)
    assert not is_legacy_session(tmp_path)
    assert not is_legacy_session(tmp_path / "missing")
    txt = tmp_path / "plain.txt"
    txt.write_text("hello")
    assert not is_legacy_session(txt)


def test_legacy_target(tmp_path):
    assert legacy_target(tmp_path / "x") == tmp_path / "x.session"


def test_read_values_warning_and_untouched(session_file):
    before = _snapshot(session_file)
    with pytest.warns(FutureWarning, match="2026-11-30"):
        values, errors = read_legacy(session_file)
    assert errors == {}
    exp = _values()
    assert set(values) == set(exp)
    for k in ("n", "l", "d"):
        assert values[k] == exp[k]
    pd.testing.assert_frame_equal(values["df"], exp["df"])
    assert _snapshot(session_file) == before


def test_read_while_another_reader_open(session_file):
    backend = "gnu" if dbm.whichdb(str(session_file)) == "dbm.gnu" else "dumb"
    other = _open_raw(backend, session_file, "r")
    try:
        with pytest.warns(FutureWarning):
            values, errors = read_legacy(session_file)
    finally:
        other.close()
    assert errors == {} and values["n"] == 3


def test_corrupt_value_lands_in_errors(session_file):
    backend = "gnu" if dbm.whichdb(str(session_file)) == "dbm.gnu" else "dumb"
    raw = _open_raw(backend, session_file, "w")
    raw["bad"] = b"\x00not a pickle"
    raw["gone"] = pickle.dumps(pickle.Pickler)  # fine
    raw["missing"] = b"\x80\x04\x95\x1c\x00\x00\x00\x00\x00\x00\x00\x8c\x0bno_such_mod\x94\x8c\x03Foo\x94\x93\x94."
    raw.close()
    with pytest.warns(FutureWarning):
        values, errors = read_legacy(session_file)
    assert set(errors) == {"bad", "missing"}
    assert all(": " in m for m in errors.values())
    assert values["n"] == 3 and "df" in values and "gone" in values


def test_errors_for_bad_paths(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_legacy(tmp_path / "nope")
    with pytest.raises(ValueError):
        read_legacy(tmp_path)
    txt = tmp_path / "plain.txt"
    txt.write_text("hello")
    with pytest.raises(ValueError):
        read_legacy(txt)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
