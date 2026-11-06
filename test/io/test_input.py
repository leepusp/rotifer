"""
Tests for the input helpers shared by the parsers in rotifer.io.
"""

import bz2
import gzip
import io

import pytest

from rotifer.io import each_source, open_input, split_fields, strip_suffix, written_as_fraction

@pytest.mark.parametrize("suffix,opener", [("", open), (".gz", gzip.open), (".bz2", bz2.open)])
def test_open_input_decompresses(tmp_path, suffix, opener):
    path = str(tmp_path / f"table.tsv{suffix}")
    with opener(path, "wt") as out:
        out.write("a\tb\n")
    with open_input(path) as handle:
        assert handle.read() == "a\tb\n"

def test_open_input_reads_standard_input(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO("x\n"))
    with open_input("-") as handle:
        assert handle.read() == "x\n"

def test_each_source(tmp_path):
    for name in ("a", "b"):
        (tmp_path / name).write_text(name)
    paths = [str(tmp_path / "a"), str(tmp_path / "b")]
    assert [(name, handle.read()) for name, handle in each_source(paths)] == list(zip(paths, "ab"))
    handle = io.StringIO("c")
    assert [h.read() for _, h in each_source(handle)] == ["c"]

def test_split_fields_drops_trailing_empty_fields():
    assert split_fields("a\t\tb\t\t", "\t") == ["a", "", "b"]
    assert split_fields("", "\t") == []

def test_strip_suffix():
    assert strip_suffix("/x/PF00001.hmmsearch", [r"\.hmm\w+"]) == "PF00001"
    assert strip_suffix("/x/name.out", [r"\.txt"]) == "name.out"

@pytest.mark.parametrize("value,expected", [("0.4", True), (".5", True), ("1e-2", True), ("10", False), (10, False)])
def test_written_as_fraction(value, expected):
    assert written_as_fraction(value) is expected
