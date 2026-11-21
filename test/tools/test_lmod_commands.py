"""
Tests for the commands that replace the shell functions of the ROTIFER
Lmod module (etc/lmod/modules/rotifer/0.3.lua): acc2pfam, acc2profiledb,
dom_count, dom_per_prot, distinct_dom_per_prot, tdesc and
aln_order_by_tree.

Golden outputs (test/data/golden/lmod) were recorded by running the
former Perl pipelines on recorded hmmscan, RPS-BLAST and Phobius reports.
"""

import collections
import io
import os
import shutil

import pytest

import golden
from rotifer.cli import write_padded
from rotifer.tools import alignments, domains, tables

def _padded(rows):
    out = io.StringIO()
    write_padded(rows, out=out)
    return out.getvalue()

def test_acc2pfam_rows_match_perl_pipeline():
    rows = domains.pfam_rows(os.path.join(golden.INPUTS, "hmmscan.out"), os.path.join(golden.INPUTS, "phobius.out"))
    assert _padded(rows) == golden.expected("lmod", "acc2pfam")

def test_acc2profiledb_rows_match_perl_pipeline():
    rows = domains.profiledb_rows(os.path.join(golden.INPUTS, "rpsblast.out"), os.path.join(golden.INPUTS, "phobius.out"))
    assert _padded(rows) == golden.expected("lmod", "acc2profiledb")

def _recorded_run(command, stdout):
    """Stand-in for hmmscan, rpsblast and phobius: copy recorded reports."""
    report = {"hmmscan": "hmmscan.out", "rpsblast": "rpsblast.out", "phobius": "phobius.out"}[command[0]]
    shutil.copy(os.path.join(golden.INPUTS, report), stdout)

def test_acc2pfam_runs_searches(monkeypatch):
    monkeypatch.setattr(domains, "_run", _recorded_run)
    df = domains.acc2pfam(fasta=golden.read("inputs/set.faa"), database="Pfam-A.hmm")
    expected = [line.split() for line in golden.expected("lmod", "acc2pfam").splitlines()[1:]]
    assert list(df.columns) == ["ID", "domain", "start", "end", "evalue"]
    assert df[["ID", "domain"]].values.tolist() == [row[:2] for row in expected]
    assert list(df.start) == [int(row[2]) for row in expected]

def test_acc2profiledb_runs_searches(monkeypatch):
    monkeypatch.setattr(domains, "_run", _recorded_run)
    df = domains.acc2profiledb(fasta=golden.read("inputs/set.faa"), database="allprofiles")
    assert len(df) == len(golden.expected("lmod", "acc2profiledb").splitlines()) - 1

def test_fetch_fasta_uses_rotifer_db(monkeypatch):
    import rotifer.db
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    calls = []
    def proteins(accessions, methods, progress):
        calls.append((accessions, methods, progress))
        return [SeqRecord(Seq("MKV"), id="NP_1.1", description="NP_1.1 a protein")]
    monkeypatch.setattr(rotifer.db, "proteins", proteins)
    assert domains.fetch_fasta("NP_1.1") == ">NP_1.1 a protein\nMKV\n"
    assert calls == [(["NP_1.1"], ["esl_sfetch", "entrez"], False)]

def _data_rows(text, skip=1):
    """Data rows as a multiset (tgroup printed them in hash order)."""
    return collections.Counter(tuple(line.split("\t")) for line in text.splitlines()[skip:])

@pytest.mark.parametrize("function", ["dom_count", "dom_per_prot", "distinct_dom_per_prot"])
def test_counts_match_perl_pipeline(function):
    table = getattr(domains, function)(os.path.join(golden.INPUTS, "pfam.arch"))
    rows = collections.Counter(tuple(row) for row in table.astype(str).values.tolist())
    # The Perl output's first row is tgroup's count of the header line
    assert rows == _data_rows(golden.expected("lmod", function))
    assert list(table.columns)[-1] == "count"

def test_tdesc():
    assert tables.tdesc(os.path.join(golden.INPUTS, "table.tsv")) == \
        [(int(i), name) for i, name in (line.split("\t") for line in golden.expected("lmod", "tdesc").splitlines())]

def test_aln_order_by_tree():
    ordered = alignments.aln_order_by_tree(os.path.join(golden.INPUTS, "aln.fa"), os.path.join(golden.INPUTS, "tree.nwk"))
    assert "".join(f">{i}\n{s}\n" for i, s in ordered) == golden.expected("lmod", "aln_order_by_tree")

@pytest.mark.parametrize("cmdline,case", [
    ("tdesc inputs/table.tsv", "tdesc"),
    ("aln_order_by_tree inputs/aln.fa inputs/tree.nwk", "aln_order_by_tree"),
])
def test_commands_match_perl_output(cmdline, case, tmp_path):
    assert golden.run(cmdline, tmp_path) == golden.expected("lmod", case)

@pytest.mark.parametrize("function", ["dom_count", "dom_per_prot", "distinct_dom_per_prot"])
def test_count_commands(function, tmp_path):
    out = golden.run(f"{function} inputs/pfam.arch", tmp_path)
    assert out.splitlines()[0].split("\t")[-1] == "count"
    assert _data_rows(out) == _data_rows(golden.expected("lmod", function))

@pytest.mark.parametrize("command", ["acc2pfam", "acc2profiledb"])
def test_search_commands_have_help(command, tmp_path):
    assert command in golden.run(f"{command} --help", tmp_path)
