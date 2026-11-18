#!/usr/bin/env python3
"""Tests for rotifer.core.functions.save_session and restore_session.

These wrappers now write a rotifer.core.session directory. The main regression
covered here: an old shelve file (gdbm) used to make save_session fail with
EAGAIN; restore_session now converts it to `<path>.session/` and leaves the
old file alone.

Run under pytest, or standalone: python test/core/session/test_wrappers.py
"""

import os
import shelve
import sys
import warnings
import threading
from pathlib import Path

import pandas as pd
import pytest

from rotifer.core.functions import restore_session, save_session
from rotifer.core.session import Session

try:
    import dbm.gnu
except ImportError:
    dbm.gnu = None


def _names(path):
    return sorted(Session(path, create=False).ls()["name"])


def test_round_trip_plain_values_and_dataframe(tmp_path):
    p = tmp_path / "session"
    ns = {"n": 3, "text": "abc", "df": pd.DataFrame({"x": [1, 2], "y": ["a", "b"]})}
    report = save_session(str(p), namespace=ns)
    assert sorted(report.saved) == ["df", "n", "text"]

    new_ns = {}
    report = restore_session(str(p), namespace=new_ns)
    assert sorted(report.loaded) == ["df", "n", "text"]
    assert new_ns["n"] == 3
    assert new_ns["text"] == "abc"
    pd.testing.assert_frame_equal(new_ns["df"], ns["df"])


def test_modules_and_imported_functions_are_not_saved(tmp_path):
    from os.path import join
    p = tmp_path / "session"
    ns = {"value": 1, "os": os, "join": join, "sys": sys}
    save_session(str(p), namespace=ns)
    assert _names(str(p)) == ["value"]


def test_exclude_applies_on_save_and_restore(tmp_path):
    p = tmp_path / "session"
    ns = {"a": 1, "b": 2, "c": 3}
    save_session(str(p), exclude=["b"], namespace=ns)
    assert _names(str(p)) == ["a", "c"]

    new_ns = {}
    restore_session(str(p), exclude=["c"], namespace=new_ns)
    assert new_ns == {"a": 1}


def test_prune_removes_variables_deleted_from_namespace(tmp_path):
    p = tmp_path / "session"
    ns = {"a": 1, "b": 2}
    save_session(str(p), namespace=ns)
    assert _names(str(p)) == ["a", "b"]

    del ns["b"]
    save_session(str(p), namespace=ns)
    assert _names(str(p)) == ["a"]


def test_failed_save_is_reported_and_others_are_saved(tmp_path):
    p = tmp_path / "session"
    ns = {"ok": 1, "bad": threading.Lock()}  # a lock cannot be pickled
    report = save_session(str(p), namespace=ns)
    assert "ok" in report.saved
    assert "bad" in report.failed
    assert _names(str(p)) == ["ok"]


@pytest.mark.skipif(dbm.gnu is None, reason="dbm.gnu not available")
def test_old_gdbm_shelve_is_converted_and_left_unchanged(tmp_path):
    import shelve
    p = tmp_path / "old_session"
    sh = shelve.Shelf(dbm.gnu.open(str(p), "c"))
    sh["a"] = 1
    sh.close()
    before = p.read_bytes()

    new_ns = {}
    with pytest.warns(FutureWarning):
        report = restore_session(str(p), namespace=new_ns)
    assert new_ns["a"] == 1
    assert report.loaded == ["a"]
    assert Path(str(p) + ".session").is_dir()
    assert p.read_bytes() == before


@pytest.mark.skipif(dbm.gnu is None, reason="dbm.gnu not available")
def test_save_session_works_after_old_gdbm_shelve_exists(tmp_path):
    import shelve
    old = tmp_path / "old_session"
    sh = shelve.Shelf(dbm.gnu.open(str(old), "c"))
    sh["a"] = 1
    sh.close()

    fresh = tmp_path / "fresh"
    report = save_session(str(fresh), namespace={"x": 2})
    assert report.saved == ["x"]
    new_ns = {}
    restore_session(str(fresh), namespace=new_ns)
    assert new_ns == {"x": 2}


def test_restore_session_missing_path_raises_and_creates_nothing(tmp_path):
    with pytest.raises(FileNotFoundError):
        restore_session(str(tmp_path / "nope"), namespace={})
    assert not (tmp_path / "nope").exists()


def test_restore_session_still_converts_legacy(tmp_path):
    old = str(tmp_path / "old")
    with shelve.open(old) as sh:
        sh["a"] = 7
    ns = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        restore_session(old, namespace=ns)
    assert ns == {"a": 7}


def test_save_session_exclude_type_prunes_stale(tmp_path):
    p = str(tmp_path / "w")
    ns = {"a": 1, "c": 5}
    save_session(p, namespace=ns)
    ns["c"] = os
    save_session(p, namespace=ns)
    assert Session(p).ls().name.tolist() == ["a"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
