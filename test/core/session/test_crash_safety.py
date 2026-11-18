"""Crash safety, conversion atomicity, prune scope and diff limits of Session."""
import os
import shelve
import warnings

import pandas as pd
import pytest

import rotifer.core.session._session as ss
from rotifer.core.session import Session


def objects(s):
    return sorted(p.name for p in (s.path / "objects").iterdir())


class Interrupt:
    def __reduce__(self):
        raise KeyboardInterrupt()


class Bad:
    def __reduce__(self):
        raise TypeError("nope")


def test_crash_before_manifest_keeps_old_session_loadable(tmp_path, monkeypatch):
    s = Session(tmp_path / "s")
    s.save({"x": pd.DataFrame({"a": [1, 2]})})
    before = objects(s)

    def boom(*a, **k):
        raise RuntimeError("crash")
    monkeypatch.setattr(ss, "write_manifest", boom)
    with pytest.raises(RuntimeError):
        s.save({"x": [1, 2, 3], "y": 5})      # codec change parquet -> pickle
    monkeypatch.undo()
    assert s.get("x").equals(pd.DataFrame({"a": [1, 2]}))
    assert objects(s) == before               # no leftovers either


def test_old_files_deleted_only_after_manifest(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"x": 1})
    old = s.info("x")["files"]
    s.save({"x": 2})
    new = s.info("x")["files"]
    assert old != new and not any((s.path / f).exists() for f in old)
    assert s.get("x") == 2


def test_failed_variable_old_files_survive_until_manifest(tmp_path, monkeypatch):
    s = Session(tmp_path / "s")
    s.save({"x": [1]})
    old = s.info("x")["files"]
    seen = {}
    real = ss.write_manifest

    def spy(path, manifest):
        seen["exists"] = all((path / f).exists() for f in old)
        real(path, manifest)
    monkeypatch.setattr(ss, "write_manifest", spy)
    rep = s.save({"x": Bad()})
    assert "x" in rep.failed and seen["exists"]
    assert not any((s.path / f).exists() for f in old)


def test_keyboard_interrupt_cleans_up_and_reraises(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"x": [1, 2, 3]})
    before = objects(s)
    with pytest.raises(KeyboardInterrupt):
        s.save({"a": 1, "x": Interrupt()})
    assert objects(s) == before
    assert s.get("x") == [1, 2, 3]
    assert "a" not in s.ls().name.tolist()


def test_remove_deletes_files_after_manifest(tmp_path, monkeypatch):
    s = Session(tmp_path / "s")
    s.save({"x": 1})
    before = objects(s)

    def boom(*a, **k):
        raise RuntimeError("crash")
    monkeypatch.setattr(ss, "write_manifest", boom)
    with pytest.raises(RuntimeError):
        s.remove("x")
    monkeypatch.undo()
    assert objects(s) == before and s.get("x") == 1
    s.remove("x")
    assert objects(s) == []


def test_prune_deletes_files_after_manifest(tmp_path, monkeypatch):
    s = Session(tmp_path / "s")
    s.save({"x": 1, "y": 2})
    before = objects(s)
    monkeypatch.setattr(ss, "write_manifest",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("crash")))
    with pytest.raises(RuntimeError):
        s.save({"y": 2}, prune=True)
    monkeypatch.undo()
    assert objects(s) == before and s.get("x") == 1


def test_sweep_removes_orphans_but_not_referenced(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"x": 1})
    (s.path / "objects" / "orphan.pkl").write_bytes(b"junk")
    (s.path / "objects" / ".tmp-left.tmp").write_bytes(b"junk")
    s.save({"y": 2})
    names = objects(s)
    assert "orphan.pkl" not in names and ".tmp-left.tmp" not in names
    assert s.get("x") == 1 and s.get("y") == 2


def test_corrupt_file_fails_that_variable_only(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"x": [1, 2, 3], "y": [4]})
    f = s.path / s.info("x")["files"][0]
    data = bytearray(f.read_bytes())
    data[-3] ^= 0xFF
    f.write_bytes(bytes(data))
    with pytest.raises(ValueError, match="corrupt"):
        s.get("x")
    out = {}
    rep = s.load(ns=out)
    assert "x" in rep.failed and "corrupt" in rep.failed["x"]
    assert out == {"y": [4]}


def test_failed_conversion_leaves_no_target(tmp_path):
    old = str(tmp_path / "old")
    with shelve.open(old) as sh:
        sh["a"] = [1, 2]
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        with pytest.raises(FutureWarning):
            Session(old)
    assert not os.path.exists(old + ".session")
    assert not [p for p in os.listdir(tmp_path) if ".tmp-" in p]
    s = Session(old)                     # now converts properly
    assert s.get("a") == [1, 2] and s.conversion_report.saved == ["a"]
    assert s.path == tmp_path / "old.session"


def test_incomplete_conversion_from_older_code_is_redone(tmp_path):
    old = str(tmp_path / "old")
    with shelve.open(old) as sh:
        sh["a"] = [1, 2]
    Session(old + ".session")            # empty rev-0 session, as an aborted conversion left
    s = Session(old)
    assert s.get("a") == [1, 2]


def test_convert_classmethod_is_atomic(tmp_path):
    old = str(tmp_path / "old")
    with shelve.open(old) as sh:
        sh["a"] = 1
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        with pytest.raises(FutureWarning):
            Session.convert(old, tmp_path / "new")
    assert not (tmp_path / "new").exists()
    assert [p for p in os.listdir(tmp_path) if ".tmp-" in p] == []


def test_prune_removes_entries_now_filtered_or_excluded(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"a": 1, "b": 2, "c": 3})
    ns = {"a": 1, "b": os, "c": 3}               # b is a module now
    rep = s.save(ns, prune=True)
    assert rep.removed == ["b"]
    rep = s.save(ns, exclude=["c"], prune=True)
    assert rep.removed == ["c"]
    assert s.ls().name.tolist() == ["a"]


def test_prune_stays_inside_names(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"a1": 1, "b1": 2})
    s.save({"a1": 1}, names="a*", prune=True)
    assert sorted(s.ls().name) == ["a1", "b1"]


def test_diff_content_set_values_are_unknown_not_changed(tmp_path):
    s = Session(tmp_path / "s")
    ns = {"st": {"abc", "def", "ghi"}, "lst": [frozenset("xy")], "d": {"k": {1, 2}},
          "plain": [1, 2]}
    s.save(ns)
    r = s.diff(ns, level="content")
    assert sorted(r.unknown) == ["d", "lst", "st"]
    assert r.unchanged == ["plain"] and r.changed == []


def test_fingerprint_not_recomputed_with_second_pickle(tmp_path):
    from rotifer.core.session.codecs import encode

    class Once:
        n = 0

        def __reduce__(self):
            Once.n += 1
            return (Once, ())
    globals()["Once"] = Once
    Once.__module__ = __name__
    Once.__qualname__ = "Once"
    s = Session(tmp_path / "s")
    s.save({"o": [Once()]})
    assert Once.n == 1
    assert s.info("o")["fingerprint"].startswith("pickle:")
