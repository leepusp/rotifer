#!/usr/bin/env python3
"""Tests for rotifer.core.session.Session. Run under pytest, or standalone:
python test/core/session/test_session.py
"""
import json
import shutil
import sys

import filelock
import pandas as pd
import pytest

from rotifer.core.session import Session, codecs, describe, open_session
from rotifer.core.session._files import SessionLocked
from rotifer.core.session._manifest import manifest_schema


def make_ns():
    return {
        "n": 42,
        "xs": [1, 2, 3],
        "d": {"a": 1, "b": [2, 3]},
        "df": pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}),
    }


def manifest(path):
    return json.loads((path / "manifest.json").read_text())


def test_roundtrip(tmp_path):
    ns = make_ns()
    s = Session(tmp_path / "s")
    rep = s.save(ns)
    assert sorted(rep.saved) == ["d", "df", "n", "xs"]
    assert not rep.failed
    assert manifest(s.path)["variables"]["df"]["codec"] == "parquet-dataframe"
    out = {}
    rep = s.load(ns=out)
    assert sorted(rep.loaded) == ["d", "df", "n", "xs"]
    assert out["n"] == 42 and out["xs"] == [1, 2, 3] and out["d"] == ns["d"]
    assert out["df"].equals(ns["df"])
    assert rep.to_dict()["action"] == "load"
    json.dumps(rep.to_dict())


def test_names_and_exclude(tmp_path):
    ns = {"a1": 1, "a2": 2, "b1": 3, "c": 4}
    s = Session(tmp_path / "s")
    rep = s.save(ns, names=["a*", "b1"], exclude="a2")
    assert sorted(rep.saved) == ["a1", "b1"]
    assert s.ls()["name"].tolist() == ["a1", "b1"]
    out = {}
    assert s.load(names="a*", ns=out).loaded == ["a1"]
    assert out == {"a1": 1}
    rep = s.save(ns, names=["zzz"])
    assert rep.skipped["zzz"] == "not in the namespace"


def test_default_filters(tmp_path):
    def local_function():
        return 1
    local_function.__module__ = "__main__"
    ns = {"_hidden": 1, "In": [], "Out": {}, "exit": 1, "quit": 1,
          "get_ipython": 1, "open": open, "sys": sys, "json_dumps": json.dumps,
          "Path": __import__("pathlib").Path, "mine": local_function, "ok": 5}
    s = Session(tmp_path / "s")
    rep = s.save(ns)
    assert sorted(rep.saved) == ["mine", "ok"]
    out = {}
    s.load(ns=out)
    assert out["mine"]() == 1


def test_add_update_and_prune(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"a": 1, "b": 2})
    rep = s.save({"a": 10})
    assert rep.removed == []
    assert sorted(s.ls()["name"]) == ["a", "b"]
    assert s.get("a") == 10
    rep = s.save({"a": 10}, prune=True)
    assert rep.removed == ["b"]
    assert s.ls()["name"].tolist() == ["a"]
    assert len(list((s.path / "objects").iterdir())) == 1
    # prune only inside the selection
    s.save({"a": 1, "b": 2})
    s.save({}, names="b*", prune=True)
    assert s.ls()["name"].tolist() == ["a"]


def test_update_removes_old_files_when_codec_changes(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"x": pd.DataFrame({"a": [1]})})
    assert s.info("x")["codec"] == "parquet-dataframe"
    s.save({"x": [1]})
    assert s.info("x")["codec"] == "pickle"
    files = [p.name for p in (s.path / "objects").iterdir()]
    assert len(files) == 1 and files[0].endswith(".pkl")


class Unsaveable:
    def __reduce__(self):
        raise RuntimeError("no way")


def test_failure_leaves_nothing_then_recovers(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"v": 1, "w": 2})
    rep = s.save({"v": Unsaveable()})
    assert "v" in rep.failed and "RuntimeError" in rep.failed["v"]
    assert s.ls()["name"].tolist() == ["w"]
    assert manifest(s.path)["failed"]["v"].startswith("RuntimeError")
    assert all(not p.name.startswith("v-") for p in (s.path / "objects").iterdir())
    with pytest.raises(KeyError):
        s.get("v")
    out = {}
    s.load(ns=out)
    assert "v" not in out
    rep = s.save({"v": 3})
    assert rep.saved == ["v"]
    assert manifest(s.path)["failed"] == {}
    assert s.get("v") == 3


def test_names_differing_by_case(tmp_path):
    s = Session(tmp_path / "s")
    a, b = pd.DataFrame({"a": [1]}), pd.DataFrame({"b": [2]})
    s.save({"df": a, "DF": b})
    assert s.get("df").equals(a) and s.get("DF").equals(b)
    files = [p.name.lower() for p in (s.path / "objects").iterdir()]
    assert len(set(files)) == len(files) == 2


def test_ls_reads_no_object_files(tmp_path):
    s = Session(tmp_path / "s")
    s.save(make_ns())
    shutil.rmtree(s.path / "objects")
    df = s.ls()
    assert list(df.columns) == ["name", "type", "codec", "bytes", "saved",
                                "version_bound", "kind", "length", "shape"]
    assert len(df) == 4
    row = df.set_index("name").loc["df"]
    assert row["kind"] == "dataframe" and list(row["shape"]) == [3, 2]
    assert Session(tmp_path / "empty").ls().empty


def test_get_info_remove(tmp_path):
    s = Session(tmp_path / "s")
    s.save(make_ns())
    ns = {}
    assert s.get("xs") == [1, 2, 3] and ns == {}
    info = s.info("df")
    assert info["type"].endswith("DataFrame") and info["bytes"] > 0
    assert len(info["sha256"]) == 64 and info["fingerprint"]
    rep = s.remove("xs")
    assert rep.removed == ["xs"]
    assert "xs" not in s.ls()["name"].tolist()
    assert not list((s.path / "objects").glob("xs-*"))
    with pytest.raises(KeyError):
        s.remove("xs")
    with pytest.raises(KeyError):
        s.info("xs")


def test_open_session(tmp_path):
    with pytest.raises(FileNotFoundError):
        open_session(tmp_path / "nope")
    Session(tmp_path / "s").save({"a": 1})
    assert open_session(tmp_path / "s").get("a") == 1


def test_lock_held_raises(tmp_path):
    s = Session(tmp_path / "s", lock_timeout=0.2)
    s.save({"a": 1})
    with filelock.FileLock(str(s.path / ".lock")):
        with pytest.raises(SessionLocked):
            s.save({"a": 2})
        with pytest.raises(SessionLocked):
            s.remove("a")
        assert s.get("a") == 1      # readers do not lock
    s.save({"a": 2})


def test_diff_levels(tmp_path):
    s = Session(tmp_path / "s")
    df = pd.DataFrame({"a": [1, 2, 3]})
    s.save({"df": df, "n": 1, "lst": [1, 2], "gone": 0})
    live = {"df": df.copy(), "n": 1, "lst": [1, 3], "new": 5}
    r = s.diff(live)
    assert r.added == ["new"] and r.removed == ["gone"]
    assert sorted(r.unchanged) == ["df", "lst", "n"] and not r.changed
    r = s.diff(live, level="summary")
    assert sorted(r.unchanged) == ["df", "lst", "n"] and not r.changed
    live["df"] = pd.DataFrame({"a": [1, 2, 3, 4]})
    r = s.diff(live, level="summary")
    assert r.changed == ["df"]
    live["df"] = df.copy()
    r = s.diff(live, level="content")
    assert r.changed == ["lst"] and sorted(r.unchanged) == ["df", "n"]
    live["df"] = pd.DataFrame({"a": [1, 2, 4]})
    assert s.diff(live, level="content").changed == ["df", "lst"]
    r = s.diff(live, level="content", max_bytes=10)
    assert "df" in r.unknown and "lst" in r.unknown and r.changed == []
    assert r.to_dict()["level"] == "content"
    with pytest.raises(ValueError):
        s.diff(live, level="everything")


def test_environment_warning_on_load(tmp_path):
    s = Session(tmp_path / "s")
    s.save({"a": 1})
    m = manifest(s.path)
    m["environment"]["packages"]["pandas"] = "0.0"
    (s.path / "manifest.json").write_text(json.dumps(m))
    rep = s.load(ns={})
    assert any("pandas" in w for w in rep.warnings)


def test_manifest_valid_and_revision(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    s = Session(tmp_path / "s")
    assert manifest(s.path)["revision"] == 0
    s.save(make_ns())
    s.save({"a": 1}, prune=True)
    m = manifest(s.path)
    jsonschema.validate(m, manifest_schema())
    assert m["revision"] == 2


def test_describe_does_not_pickle():
    class NoPickle:
        def __reduce__(self):
            raise RuntimeError("describe must not pickle")
    assert describe(NoPickle())["kind"] == "object"
    assert codecs.find_typed_codec(NoPickle()) is None
    assert codecs.find_typed_codec(pd.DataFrame()).name == "parquet-dataframe"


# ------------------------------------------------------------------ legacy

def _make_shelve(path, backend):
    import importlib
    import shelve
    mod = importlib.import_module(backend)
    with shelve.Shelf(mod.open(str(path), "n")) as db:
        db["n"] = 5
        db["df"] = pd.DataFrame({"a": [1, 2]})
        db["l"] = [1, 2, 3]


def _snapshot(directory, stem):
    snap = {}
    for p in sorted(directory.iterdir()):
        if p.name.startswith(stem) and p.is_file():
            snap[p.name] = (p.read_bytes(), p.stat().st_mtime_ns)
    return snap


@pytest.mark.parametrize("backend", ["dbm.dumb", "dbm.gnu"])
def test_legacy_conversion(tmp_path, backend):
    pytest.importorskip("rotifer.core.session._legacy")
    pytest.importorskip(backend)
    old = tmp_path / "old"
    _make_shelve(old, backend)
    before = _snapshot(tmp_path, "old")
    assert before
    with pytest.warns(FutureWarning):
        s = Session(old)
    assert s.path == tmp_path / "old.session"
    assert s.converted_from == str(old)
    assert s.conversion_report.converted_from == str(old)
    assert s.get("n") == 5 and s.get("l") == [1, 2, 3]
    assert s.get("df").equals(pd.DataFrame({"a": [1, 2]}))
    assert _snapshot(tmp_path, "old") == {k: v for k, v in before.items()}
    # second open: already converted, opens the directory
    s2 = Session(old)
    assert s2.path == s.path and s2.get("n") == 5


def test_legacy_explicit_convert(tmp_path):
    pytest.importorskip("rotifer.core.session._legacy")
    old = tmp_path / "old"
    _make_shelve(old, "dbm.dumb")
    with pytest.warns(FutureWarning):
        s = Session.convert(old, tmp_path / "new")
    assert s.path == tmp_path / "new" and s.get("n") == 5


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
