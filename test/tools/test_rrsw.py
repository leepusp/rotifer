"""
Tests for rotifer.tools.rrsw and the rrsw command (formerly Perl), which
the cluster uses to mirror UniProt (``rrsw uniprot``).

These tests run offline. While porting, the Python and Perl versions
were also compared on the real servers: ``rrsw -n uniprot`` selected the
same 130 files, and a small FTP mirror plus its unpacked tree were
identical.
"""

import gzip
import os
import stat
import subprocess
import tarfile

import pytest

import golden
from rotifer.tools import rrsw

def test_default_table_targets():
    table = rrsw.load_table()
    assert table[0] == rrsw.TABLE_COLUMNS
    rows = {row[0]: row for row in table[1:]}
    assert rows["uniprot"][:5] == ["uniprot", "ftp", "ftp.uniprot.org",
                                   "pub/databases/uniprot/current_release", "unpack"]
    assert rows["blast"][5]["callback"] == ["process_blast_dbs.sh"]
    assert set(row[4] for row in table[1:] if row[4]) <= set(rrsw.PROCESSORS) | {"none"}
    assert set(row[1] for row in table[1:]) <= set(rrsw.METHODS)

def test_all_target_configurations_load():
    for name in os.listdir(rrsw.config_dir()):
        if name.endswith(".yml"):
            rrsw.load_yaml(os.path.join(rrsw.config_dir(), name))
    assert rrsw.load_yaml(os.path.join(rrsw.config_dir(), "eggNOG.yml"))["exclude"][0] == r".*/all\..*|^all\..*"

def test_select_targets():
    table = [rrsw.TABLE_COLUMNS, ["gene"] + [None] * 5, ["uniprot"] + [None] * 5, ["pdb"] + [None] * 5]
    assert [row[0] for row in rrsw.select_targets(table, ["UniProt"])] == ["NAME", "uniprot"]
    assert [row[0] for row in rrsw.select_targets(table, ["all"], exclude=["pdb"])] == ["NAME", "gene", "uniprot"]

def test_parse_listing():
    lines = [
        "drwxr-xr-x    4 ftp      ftp          4096 Oct 15 10:00 knowledgebase",
        "-rw-r--r--    1 ftp      ftp          3412 Oct 15  2025 README",
        "lrwxrwxrwx    1 ftp      ftp            16 Oct 15 10:00 current -> release-2026_04",
        "total 12",
    ]
    assert rrsw.parse_listing(lines) == [("knowledgebase", "d"), ("README", "f"),
                                         ("current", "l release-2026_04")]

def test_ftp_filters(tmp_path):
    config = tmp_path / "filters.yml"
    config.write_text("---\ninclude:\n - README\n - uniref\nexclude:\n - \\.xml\\.gz$\n")
    ftp = rrsw.Ftp("nonexistent-target", options={"config": str(config)})
    assert not ftp._ignore("/pub/README")
    assert not ftp._ignore("/pub/uniref/uniref50.fasta.gz")
    assert ftp._ignore("/pub/uniref/uniref50.xml.gz")
    assert ftp._ignore("/pub/news.html")

def _downloaded_tree(root):
    source = root / "src" / "local" / "fake"
    (source / "sub" / "deep").mkdir(parents=True)
    (source / "README").write_text("hello\n")
    (source / "table.tsv.gz").write_bytes(gzip.compress(b"a\nb\n"))
    (source / "sub" / "deep" / "notes.txt").write_text("deep\n")
    inner = root / "tarball" / "inner"
    inner.mkdir(parents=True)
    (inner / "t.txt").write_text("t\n")
    with tarfile.open(source / "sub" / "archive.tar.gz", "w:gz") as tar:
        tar.add(root / "tarball", arcname="tarball")
    return source

def test_unpack(tmp_path):
    source = _downloaded_tree(tmp_path)
    destination = tmp_path / "fake"
    destination.mkdir()
    (destination / "old.txt").write_text("previous copy\n")
    rrsw.Unpack("fake").execute(str(source) + "/", str(destination))
    assert (destination / "table.tsv").read_text() == "a\nb\n"      # decompressed
    assert os.readlink(destination / "README") == str(source / "README")  # linked
    assert os.readlink(destination / "sub" / "deep" / "notes.txt") == str(source / "sub" / "deep" / "notes.txt")
    assert (destination / "tarball" / "inner" / "t.txt").read_text() == "t\n"  # extracted at the top
    assert stat.S_IMODE(os.stat(destination).st_mode) == 0o775
    previous = [p for p in tmp_path.iterdir() if p.name.startswith("fake.")]
    assert len(previous) == 1 and (previous[0] / "old.txt").exists()

def test_unpack_test_mode_changes_nothing(tmp_path, capsys):
    source = _downloaded_tree(tmp_path)
    rrsw.Unpack("fake", test=True).execute(str(source), str(tmp_path / "fake"))
    assert not (tmp_path / "fake").exists()
    assert "table.tsv.gz -dc >table.tsv" in capsys.readouterr().out

def test_unpack_runs_shell_callbacks(tmp_path, monkeypatch):
    import rotifer
    scripts = tmp_path / "user" / "lib" / "rrsw"
    scripts.mkdir(parents=True)
    callback = scripts / "mark.sh"
    callback.write_text('#!/bin/bash\necho "$1 $(basename $2) $(basename $3)" > callback.txt\n')
    callback.chmod(0o755)
    monkeypatch.setitem(rotifer.config, "user", str(tmp_path / "user"))
    source = _downloaded_tree(tmp_path)
    destination = tmp_path / "fake"
    rrsw.Unpack("fake", options={"callback": ["mark.sh"]}).execute(str(source), str(destination))
    assert (destination / "callback.txt").read_text() == "0 fake fake\n"

def test_perl_callbacks_are_rejected(tmp_path):
    with pytest.raises(ValueError):
        rrsw.Unpack("fake", options={"callback": ["sub { 1 }"]}).execute(str(tmp_path), str(tmp_path / "out"))

def test_blast_callback_replaces_gaps(tmp_path):
    """process_blast_dbs.sh (formerly a perl one-liner, broken by eval)."""
    script = os.path.join(golden.REPO, "lib", "rotifer", "rrsw", "process_blast_dbs.sh")
    work = tmp_path / "blast"
    (work / "db" / "FASTA").mkdir(parents=True)
    (work / "db" / "FASTA" / "pdbaa").write_text(">a-b\nAC-GT\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "esl-sfetch").write_text("#!/bin/sh\n")
    (bin_dir / "esl-sfetch").chmod(0o755)
    env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
    subprocess.run(["bash", script, "0", str(tmp_path), str(tmp_path / "blast")], cwd=work, env=env, check=True)
    assert (work / "db" / "FASTA" / "pdbaa").read_text() == ">a-b\nACXGT\n"

def test_command_lists_targets(tmp_path):
    out = golden.run("rrsw -l uniprot", tmp_path)
    assert "ftp.uniprot.org" in out and "gene" not in out

def test_command_without_targets_shows_help(tmp_path):
    result = subprocess.run(["bash", "-c", "rrsw"], cwd=golden.GOLDEN, env=golden.environ(tmp_path),
                            capture_output=True, text=True)
    assert result.returncode == 1 and "usage: rrsw" in result.stdout

def test_command_reports_progress(monkeypatch):
    import logging
    from rotifer.cli import rrsw as command
    monkeypatch.setattr(rrsw, "rrsw", lambda *args, **kwargs: None)
    logger = logging.getLogger("rotifer.tools.rrsw")
    monkeypatch.setattr(logger, "level", logging.NOTSET)
    assert command.main(["-n", "uniprot"]) == 0
    assert logger.isEnabledFor(logging.INFO)
