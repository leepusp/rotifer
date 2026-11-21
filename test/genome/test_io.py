"""
Tests for rotifer.genome.io: reading genome files into SeqRecords.

What has to hold, and why
-------------------------
* ``parse`` accepts any callable as ``assembly`` (a lambda, a
  ``functools.partial``, a callable object, a builtin). Before, only plain
  functions were recognised and anything else gave ``Assembly:Assembly:...``.
  A truthy value that is neither callable nor a string raises TypeError.
* ``_get_basename`` never fails on a handle without a file name
  (``io.StringIO``); it returns ``'__stream__'``. Paths, ``sys.stdin`` and
  real files keep their old results.
* ``gff`` reads a file with no ``##FASTA`` section (features, no sequence)
  and still returns ``[]`` for an empty file.
* ``gff(organism=taxid)`` takes the organism name and taxonomy from
  ``rotifer.db.ncbi.TaxonomyCursor``. Organism names raise ValueError.
* Files opened by ``parse`` and ``gff`` are closed; handles passed in by
  the caller are left open.
* ``parse`` accepts a list of inputs, renames number-only ids to
  ``<basename>:<id>`` and fills an empty description.

Fixture description
-------------------
Everything is built in ``tmp_path``: a two-record GenBank file named like an
NCBI assembly file (``GCF_000005845.2_toy_genomic.gbff``), and small GFF
files written from strings. The taxonomy cursor is replaced by a stub, so
no test reaches the network or the local ETE database.

Run under pytest: PYTHONPATH=lib pytest test/genome/test_io.py
"""
import functools
import io
import sys
import warnings

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord

from rotifer.genome import io as rgio

GFF_BODY = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
    "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1\n"
)
GFF_FASTA = GFF_BODY + "##FASTA\n>chr1\nATGAAATAA\n>chr2\nGGGG\n"


def make_record(rid, name=None, description="toy", dbxref=None):
    """Return a small DNA record, optionally with an Assembly dbxref."""
    rec = SeqRecord(Seq("ATGAAATAA"), id=rid, name=name or rid, description=description)
    rec.annotations["molecule_type"] = "DNA"
    if dbxref:
        rec.dbxrefs = [dbxref]
    return rec


@pytest.fixture
def gbff(tmp_path):
    """Write a two-record GenBank file named like an NCBI assembly file."""
    path = tmp_path / "GCF_000005845.2_toy_genomic.gbff"
    SeqIO.write([make_record("NC_1"), make_record("NC_2")], path, "genbank")
    return path


@pytest.fixture
def gff_file(tmp_path):
    """Write a GFF file with a ##FASTA section."""
    path = tmp_path / "with_fasta.gff"
    path.write_text(GFF_FASTA)
    return path


@pytest.fixture
def gff_nofasta(tmp_path):
    """Write a GFF file without a ##FASTA section."""
    path = tmp_path / "no_fasta.gff"
    path.write_text(GFF_BODY)
    return path


class Rule:
    """A callable object, to test that any callable works as ``assembly``."""

    def __call__(self, path, s):
        """Return a fixed identifier."""
        return "OBJ"


def dbxrefs(records):
    """Return the dbxrefs lists of ``records``."""
    return [r.dbxrefs for r in records]


# --- parse: assembly rules ------------------------------------------------------

def test_default_rule_reads_assembly_from_filename(gbff):
    """The default rule takes GCF_000005845.2 from the file name."""
    assert dbxrefs(rgio.parse(str(gbff), "genbank")) == [["Assembly:GCF_000005845.2"]] * 2


@pytest.mark.parametrize("rule", [
    lambda p, s: "X",
    functools.partial(lambda p, s, x: x, x="X"),
    Rule(),
], ids=["lambda", "partial", "object"])
def test_any_callable_is_an_assembly_rule(gbff, rule):
    """Lambdas, partials and callable objects are all called, no doubled prefix."""
    expected = "OBJ" if isinstance(rule, Rule) else "X"
    assert dbxrefs(rgio.parse(str(gbff), "genbank", assembly=rule)) == [[f"Assembly:{expected}"]] * 2


def test_string_assembly_is_used_for_every_record(gbff):
    """A string rule is used as is."""
    assert dbxrefs(rgio.parse(str(gbff), "genbank", assembly="GCA_1.1")) == [["Assembly:GCA_1.1"]] * 2


def test_false_assembly_keeps_record_value(tmp_path):
    """None keeps the file's Assembly dbxref, or the bare sequence name if it has none."""
    path = tmp_path / "x.gbk"
    SeqIO.write([make_record("A1", dbxref="Assembly:GCF_9.1"), make_record("B1")], path, "genbank")
    assert dbxrefs(rgio.parse(str(path), "genbank", assembly=None)) == [["Assembly:GCF_9.1"], ["B1"]]


@pytest.mark.parametrize("bad", [5, 1.5, ["x"], object()])
def test_non_callable_non_string_assembly_raises(gbff, bad):
    """A truthy value that is not callable or a string raises TypeError."""
    with pytest.raises(TypeError, match="assembly must be"):
        list(rgio.parse(str(gbff), "genbank", assembly=bad))


# --- parse: inputs and record fixes --------------------------------------------

def test_parse_list_of_paths(gbff, tmp_path):
    """A list of paths is read in order."""
    other = tmp_path / "GCA_000000001.1_other.gbk"
    SeqIO.write([make_record("OT_1")], other, "genbank")
    out = list(rgio.parse([str(gbff), str(other)], "genbank"))
    assert [r.id for r in out] == ["NC_1", "NC_2", "OT_1"]
    assert [r.dbxrefs for r in out] == [["Assembly:GCF_000005845.2"]] * 2 + [["Assembly:GCA_000000001.1"]]


def test_parse_open_handle_stays_open(gbff):
    """A handle given by the caller is read but not closed."""
    with open(gbff) as h:
        out = list(rgio.parse(h, "genbank"))
        assert not h.closed
    assert [r.id for r in out] == ["NC_1", "NC_2"]
    assert out[0].dbxrefs == ["Assembly:GCF_000005845.2"]


def test_parse_stringio_uses_stream_placeholder(gbff):
    """A StringIO has no name, so the default rule returns '__stream__'."""
    h = io.StringIO(gbff.read_text())
    assert dbxrefs(rgio.parse(h, "genbank")) == [["Assembly:__stream__"]] * 2


def test_number_only_id_is_renamed(tmp_path):
    """An id of digits becomes <basename>:<id>, for id and name."""
    path = tmp_path / "mygenome.gbk"
    SeqIO.write([make_record("12345", name="12345")], path, "genbank")
    (rec,) = rgio.parse(str(path), "genbank")
    assert (rec.id, rec.name) == ("mygenome:12345", "mygenome:12345")
    assert rec.dbxrefs == ["Assembly:mygenome"]


def test_empty_description_is_filled(monkeypatch):
    """An empty description becomes '<name>, whole genome shotgun sequence'."""
    monkeypatch.setattr(rgio.SeqIO, "parse", lambda h, f: iter([make_record("NC_7", description="")]))
    (out,) = rgio.parse(io.StringIO(""), "genbank")
    assert out.description == "NC_7, whole genome shotgun sequence"


def test_parse_closes_files_it_opens(gbff):
    """Reading a path to the end leaves no open file (no ResourceWarning)."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", ResourceWarning)
        list(rgio.parse(str(gbff), "genbank"))
        import gc
        gc.collect()


# --- _get_basename -------------------------------------------------------------

def test_basename_of_paths(tmp_path):
    """Paths (existing or not) give the name without directory and last extension."""
    assert rgio._get_basename("/no/such/dir/a_b.gbff") == "a_b"
    assert rgio._get_basename("a.gbff.gz") == "a.gbff"
    p = tmp_path / "real.gbk"
    p.write_text("x")
    assert rgio._get_basename(str(p)) == "real"
    assert rgio._get_basename(p) == "real"


def test_basename_of_handles(tmp_path):
    """Handles with a name use it; stdin and nameless handles use placeholders."""
    p = tmp_path / "real.gbk"
    p.write_text("x")
    with open(p) as h:
        assert rgio._get_basename(h) == "real"
    with open(p, "rb") as h:
        assert rgio._get_basename(h) == "real"
    assert rgio._get_basename(sys.stdin) == "__stdin__"
    assert rgio._get_basename(io.StringIO("x")) == "__stream__"
    assert rgio._get_basename(io.BytesIO(b"x")) == "__stream__"


def test_basename_of_unknown_object_raises():
    """Something that is neither a path nor a handle raises TypeError."""
    with pytest.raises(TypeError, match="Unknown file-like object"):
        rgio._get_basename(3.5)


# --- gff -----------------------------------------------------------------------

def test_gff_without_fasta_section(gff_nofasta):
    """Features come back with no sequence data."""
    (rec,) = rgio.gff(str(gff_nofasta))
    assert rec.id == "chr1"
    assert [f.type for f in rec.features] == ["gene"]
    assert [f.type for f in rec.features[0].sub_features] == ["CDS"]
    assert rec.annotations["topology"] == "linear"
    with pytest.raises(UndefinedSequenceError):
        str(rec.seq)


def test_gff_with_fasta_section(gff_file):
    """Sequences are attached, and a sequence without features is returned too."""
    out = rgio.gff(str(gff_file))
    assert [r.id for r in out] == ["chr1", "chr2"]
    assert str(out[0].seq) == "ATGAAATAA"
    assert str(out[1].seq) == "GGGG"
    assert len(out[0].features) == 1 and len(out[1].features) == 0


def test_gff_empty_file(tmp_path):
    """An empty file gives an empty list."""
    p = tmp_path / "empty.gff"
    p.write_text("")
    assert rgio.gff(str(p)) == []
    assert rgio.gff(io.StringIO("")) == []


def test_gff_annotations(gff_nofasta):
    """Topology, keywords, description and seqid are applied."""
    kw = ["a", "b"]
    (rec,) = rgio.gff(str(gff_nofasta), topology="circular", keywords=kw,
                      description=lambda p, s: "d:" + s.id, seqid=lambda s, p: "new")
    assert rec.annotations["topology"] == "circular"
    assert rec.annotations["keywords"] == ["a", "b"]
    assert rec.description == "d:new" and rec.id == "new"
    assert kw == ["a", "b"]
    assert re_date(rec.annotations["date"])


def re_date(text):
    """Return True when ``text`` looks like 05-MAR-2025."""
    import re
    return bool(re.match(r"^\d{2}-[A-Z]{3}-\d{4}$", text))


def test_gff_from_stringio(gff_file):
    """A text handle without a name works, and stays open."""
    h = io.StringIO(gff_file.read_text())
    out = rgio.gff(h)
    assert not h.closed
    assert [r.id for r in out] == ["chr1", "chr2"]


def test_gff_through_parse(gff_nofasta):
    """parse(..., 'gff') reads the file and applies the assembly rule."""
    (rec,) = rgio.parse(str(gff_nofasta), "gff")
    assert rec.id == "chr1"
    assert rec.dbxrefs == ["Assembly:no_fasta"]


class StubTaxonomyCursor:
    """Stand-in for rotifer.db.ncbi.TaxonomyCursor; records the requests."""

    calls = []

    def __init__(self, *args, **kwargs):
        """Accept the arguments of the real cursor."""

    def fetchall(self, accessions, *args, **kwargs):
        """Return one row for taxid 562, nothing for other taxids."""
        StubTaxonomyCursor.calls.append(list(accessions))
        if list(accessions) != [562]:
            return pd.DataFrame(columns=["taxid", "organism", "classification"])
        return pd.DataFrame({
            "taxid": [562],
            "organism": ["Escherichia coli"],
            "classification": ["Bacteria; Pseudomonadota; Gammaproteobacteria; Escherichia coli"],
        })


@pytest.fixture
def stub_taxonomy(monkeypatch):
    """Replace TaxonomyCursor so that no database or network is used."""
    import rotifer.db.ncbi as ncbi
    StubTaxonomyCursor.calls = []
    monkeypatch.setattr(ncbi, "TaxonomyCursor", StubTaxonomyCursor)
    return StubTaxonomyCursor


@pytest.mark.parametrize("taxid", [562, "562"])
def test_gff_organism_taxid(gff_nofasta, stub_taxonomy, taxid):
    """The organism, source and taxonomy annotations come from the cursor."""
    (rec,) = rgio.gff(str(gff_nofasta), organism=taxid, strain=" K-12")
    assert rec.annotations["organism"] == "Escherichia coli K-12"
    assert rec.annotations["source"] == "Escherichia coli K-12"
    assert rec.annotations["taxonomy"] == [
        "Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Escherichia coli"]
    assert stub_taxonomy.calls == [[562]]


def test_gff_organism_unknown_taxid(gff_nofasta, stub_taxonomy):
    """A taxid the cursor does not find raises ValueError."""
    with pytest.raises(ValueError, match="not found"):
        rgio.gff(str(gff_nofasta), organism=999)


def test_gff_organism_name_raises(gff_nofasta, stub_taxonomy):
    """Organism names are not supported and the cursor is not called."""
    with pytest.raises(ValueError, match="taxid"):
        rgio.gff(str(gff_nofasta), organism="Escherichia coli")
    assert stub_taxonomy.calls == []


def test_gff_closes_files_it_opens(gff_file):
    """A path is opened and closed inside gff (no ResourceWarning)."""
    import gc
    with warnings.catch_warnings():
        warnings.simplefilter("error", ResourceWarning)
        rgio.gff(str(gff_file))
        gc.collect()


def test_gff_leaves_caller_handle_open(gff_file):
    """A handle given by the caller is not closed."""
    with open(gff_file) as h:
        rgio.gff(h)
        assert not h.closed
