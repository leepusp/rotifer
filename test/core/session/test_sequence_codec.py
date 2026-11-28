"""Tests for the rotifer.sequence session codec. Run with pytest or as a script."""
import os
import subprocess
import sys
import textwrap

import pandas as pd
import pytest

pytest.importorskip("pyarrow")
from rotifer.core.session import codecs, describe
from rotifer.devel.beta.sequence import sequence

FASTA = ">Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n>Seq3\nACFH--GHS\n"


def make():
    aln = sequence(FASTA)
    aln.add_consensus()
    aln.df["score"] = [float(i) for i in range(len(aln.df))]
    return aln


def test_lookup():
    assert codecs.lookup(sequence(FASTA)).name == "rotifer.sequence"


def test_round_trip(tmp_path):
    aln = make()
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    assert codec.name == "rotifer.sequence"
    assert meta["tables"] == {"df": "parquet"}
    assert sorted(meta["files"]) == ["aln.df.parquet", "aln.json"]
    back = codecs.decode(codec.name, tmp_path / "aln", meta)
    assert back == aln
    assert back.df.equals(aln.df)
    assert list(back.df.columns) == list(aln.df.columns)
    assert back.name == aln.name and back.file_path == aln.file_path
    assert back.input_format == aln.input_format
    assert back.cursors is not aln.cursors
    for key in aln.cursors:
        assert back.cursors[key] is not aln.cursors[key]


def test_ndf_round_trip(tmp_path):
    aln = make()
    aln.ndf = pd.DataFrame({"pid": ["a", "b"], "start": [1, 5]})
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    assert meta["tables"] == {"df": "parquet", "ndf": "parquet"}
    back = codecs.decode(codec.name, tmp_path / "aln", meta)
    assert back.ndf.equals(aln.ndf)


def test_no_module_path_in_files(tmp_path):
    aln = make()
    aln.ndf = pd.DataFrame({"pid": ["a"]})
    _, meta = codecs.encode(aln, tmp_path / "aln")
    for name in meta["files"]:
        assert b"rotifer.devel" not in (tmp_path / name).read_bytes()
    assert "rotifer.devel" not in str(meta)


def test_pickle_fallback(tmp_path):
    aln = make()
    aln.df["mixed"] = [[1, "a", None, (1, 2)][i % 4] for i in range(len(aln.df))]
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    assert codec.name == "rotifer.sequence"
    assert meta["tables"]["df"] == "pickle"
    assert "aln.df.pkl" in meta["files"]
    assert not (tmp_path / "aln.df.parquet").exists()
    back = codecs.decode(codec.name, tmp_path / "aln", meta)
    assert back == aln
    assert back.df.equals(aln.df)
    assert list(back.df["mixed"]) == list(aln.df["mixed"])


def test_path_attribute_keeps_type(tmp_path):
    from pathlib import Path
    for value in (Path("/x/y.fa"), "/x/y.fa", None):
        aln = make()
        aln.file_path = value
        codec, meta = codecs.encode(aln, tmp_path / "aln")
        back = codecs.decode(codec.name, tmp_path / "aln", meta)
        assert back.file_path == value and type(back.file_path) is type(value)


def test_unserializable_attribute_raises_at_write(tmp_path):
    from rotifer.devel.beta.sequence import _sequence_session_write
    aln = make()
    aln.name = object()
    with pytest.raises(TypeError):
        _sequence_session_write(aln, tmp_path / "aln")


def test_tampered_annotation_detected(tmp_path):
    aln = make()
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    other = make()
    other.df["score"] = 99.0
    other.df.to_parquet(tmp_path / meta["table_files"]["df"], engine="pyarrow")
    with pytest.raises(ValueError):
        codecs.decode(codec.name, tmp_path / "aln", meta)


def test_tampered_ndf_detected(tmp_path):
    aln = make()
    aln.ndf = pd.DataFrame({"pid": ["a", "b"], "start": [1, 5]})
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    pd.DataFrame({"pid": ["a", "b"], "start": [1, 6]}).to_parquet(
        tmp_path / meta["table_files"]["ndf"], engine="pyarrow")
    with pytest.raises(ValueError):
        codecs.decode(codec.name, tmp_path / "aln", meta)


def test_table_lookup_uses_recorded_names(tmp_path):
    aln = make()
    codec, meta = codecs.encode(aln, tmp_path / "my.df.v2")
    assert meta["table_files"] == {"df": "my.df.v2.df.parquet"}
    back = codecs.decode(codec.name, tmp_path / "my.df.v2", meta)
    assert back.df.equals(aln.df)


def test_index_freq_falls_back_to_pickle(tmp_path):
    aln = make()
    aln.ndf = pd.DataFrame({"v": range(3)},
                           index=pd.date_range("2020-01-01", periods=3, freq="D"))
    codec, meta = codecs.encode(aln, tmp_path / "aln")
    back = codecs.decode(codec.name, tmp_path / "aln", meta)
    assert back.ndf.equals(aln.ndf)
    assert back.ndf.index.freq == aln.ndf.index.freq


def test_describe():
    info = codecs.lookup(make()).describe(make())
    assert info["kind"] == "sequence"
    assert info["length"] == 3
    assert info["alignment_length"] == 9
    assert info["ids"][0] == "Seq1"


def test_session_round_trip(tmp_path):
    session = pytest.importorskip("rotifer.core.session")
    if not hasattr(session, "Session"):
        pytest.skip("rotifer.core.session.Session not available yet")
    aln = make()
    s = session.Session(tmp_path / "s.session")
    report = s.save({"aln": aln})
    assert "aln" in report.saved
    assert s.info("aln")["codec"] == "rotifer.sequence"
    ns = {}
    session.Session(tmp_path / "s.session").load(ns=ns)
    assert ns["aln"] == aln
    assert ns["aln"].df.equals(aln.df)


def test_sequence_codec_loads_in_fresh_process(tmp_path):
    save = textwrap.dedent("""
        import sys
        from rotifer.devel.beta.sequence import sequence
        from rotifer.core.session import Session
        s = sequence(">a\\nAC-G\\n>b\\nACTG\\n")
        r = Session(sys.argv[1]).save({"s": s})
        assert r.saved == ["s"], r.to_dict()
    """)
    load = textwrap.dedent("""
        import sys
        from rotifer.core.session import Session
        assert "rotifer.devel.beta.sequence" not in sys.modules
        v = Session(sys.argv[1]).get("s")
        print(type(v).__name__, len(v.df))
    """)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
    path = str(tmp_path / "sess")
    subprocess.run([sys.executable, "-c", save, path], check=True, env=env)
    out = subprocess.run([sys.executable, "-c", load, path], check=True, env=env,
                         capture_output=True, text=True)
    assert out.stdout.split() == ["sequence", "2"]


def test_describe_uses_typed_codec_describe():
    pytest.importorskip("pyarrow")
    from rotifer.devel.beta.sequence import sequence
    d = describe(sequence(">A\nACD\n>B\nACE\n"))
    assert d["kind"] == "sequence" and d["length"] == 2


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
