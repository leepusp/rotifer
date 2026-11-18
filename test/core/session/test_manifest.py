#!/usr/bin/env python3
"""Tests for rotifer.core.session._manifest and the schema CLI.

Run under pytest, or standalone: python test/core/session/test_manifest.py
"""

import json
import subprocess
import sys

import pytest

from rotifer.core.session._manifest import (
    Manifest, ManifestError, VariableEntry, manifest_schema, read_manifest,
    write_manifest)
import rotifer.core.session._manifest as mod


def sample():
    entry = VariableEntry(
        type="pandas.core.frame.DataFrame", codec="parquet", codec_version=1,
        files=["objects/df-abc.parquet"], bytes=10, sha256="ff",
        fingerprint="1234", saved="2026-01-01T00:00:00+00:00",
        summary={"kind": "dataframe", "shape": [2, 1], "columns": ["a"],
                 "dtypes": {"a": "int64"}, "preview": "x"},
        meta={"compression": None})
    return Manifest(
        revision=2, created="2026-01-01T00:00:00+00:00",
        updated="2026-01-02T00:00:00+00:00",
        environment={"python": "3.11.0", "rotifer": None,
                     "packages": {"pandas": "2.0", "pyarrow": None},
                     "tools": {}},
        context={"history_session": 3},
        variables={"df": entry}, failed={"bad": "cannot pickle"})


def test_schema_validates_sample():
    jsonschema = pytest.importorskip("jsonschema")
    jsonschema.Draft202012Validator.check_schema(manifest_schema())
    jsonschema.validate(sample().to_dict(), manifest_schema())
    jsonschema.validate(Manifest().to_dict(), manifest_schema())


def test_schema_rejects_bad_manifest():
    jsonschema = pytest.importorskip("jsonschema")
    bad = sample().to_dict()
    del bad["variables"]["df"]["codec"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, manifest_schema())


def test_dict_round_trip():
    m = sample()
    assert Manifest.from_dict(json.loads(json.dumps(m.to_dict()))) == m


def test_disk_round_trip(tmp_path):
    m = sample()
    write_manifest(tmp_path, m)
    assert read_manifest(tmp_path) == m
    assert [p.name for p in tmp_path.iterdir()] == ["manifest.json"]


def test_newer_version_is_an_error():
    data = sample().to_dict()
    data["version"] = 99
    with pytest.raises(ManifestError, match="newer"):
        Manifest.from_dict(data)


def test_wrong_format_is_an_error():
    data = sample().to_dict()
    data["format"] = "other"
    with pytest.raises(ManifestError):
        Manifest.from_dict(data)


def test_upgrade_chain_applied(monkeypatch):
    monkeypatch.setattr(mod, "VERSION", 2)
    monkeypatch.setitem(mod._UPGRADES, 1, lambda d: {**d, "revision": 41})
    data = sample().to_dict()
    data["version"] = 1
    m = Manifest.from_dict(data)
    assert (m.version, m.revision) == (2, 41)


def test_unknown_variable_keys_are_ignored():
    data = sample().to_dict()
    data["variables"]["df"]["future_field"] = 1
    assert Manifest.from_dict(data).variables["df"].codec == "parquet"


def test_cli_schema():
    jsonschema = pytest.importorskip("jsonschema")
    out = subprocess.run([sys.executable, "-m", "rotifer.core.session", "schema"],
                         capture_output=True, text=True, check=True)
    jsonschema.Draft202012Validator.check_schema(json.loads(out.stdout))


def test_cli_unknown_command():
    out = subprocess.run([sys.executable, "-m", "rotifer.core.session", "nope"],
                         capture_output=True, text=True)
    assert out.returncode == 2 and out.stdout == "" and "usage" in out.stderr


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
