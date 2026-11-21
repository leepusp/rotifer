#!/usr/bin/env python3
"""Byte-for-byte equivalence of rotifer.genome with the code it replaced.

What has to hold, and why:

  * seqrecords_to_dataframe and io.parse give the same table and the same
    SeqRecord identifiers as before. Every NCBI, SQLite and parallel cursor
    builds its tables with them, so any difference reaches users
  * NeighborhoodDF.neighbors and vicinity give the same rows, columns,
    dtypes, index and values for every option the cursors pass (before,
    after, min_block_distance, strand None or 'same', fttype 'same'), on
    linear and circular replicons, with queries near the origin and next to
    each other
  * boundaries, group_strand and same_strand are unchanged
  * select_neighbors (now a shim over trim_blocks) gives the old output
    wherever the old code was right: windows that do not overlap, a default
    RangeIndex. Its overlap and index bugs are pinned in test_trim_blocks.py

Two checks compare two results, strictest last:

  1. pandas.testing.assert_frame_equal with every check on, except the frame
     class, which is a different (frozen) NeighborhoodDF class on each side
  2. _canonical_bytes, a readable encoding of everything a caller can see
     (class name, column names and dtypes, index, the type and repr of every
     cell), so a failure shows which part differs

The reference is _legacy_data.py, _legacy_utils.py and _legacy_io.py, frozen
copies of the modules before the refactor, and the files in data/genome, the
old modules' output on fixed generated genomes. The golden files are compared
as raw bytes, so they keep guarding the output after the legacy copies are
deleted. Regenerate them only from the legacy code:
ROTIFER_WRITE_GENOME_GOLDEN=1 pytest test/genome/test_genome_golden.py

Do not edit these tests to make a change pass. A change that breaks them
changes user-visible output and needs its own test and an entry in the PR's
"Output that changes" list.

The data are generated with a seeded generator; nothing touches the network.

Run under pytest: PYTHONPATH=lib pytest test/genome/test_genome_golden.py
"""

import importlib.util
import itertools
import os
import pathlib
import random
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from rotifer.genome import io as current_io
from rotifer.genome import utils as current_utils
from rotifer.genome.data import NeighborhoodDF

HERE = pathlib.Path(__file__).resolve().parent
GOLDEN = HERE.parent / "data" / "genome"
WRITE_GOLDEN = os.environ.get("ROTIFER_WRITE_GENOME_GOLDEN") == "1"

pytestmark = pytest.mark.filterwarnings("ignore")


def _load(name, filename):
    """Import a frozen module by path under `name`; test/ is not a package."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    legacy_utils = _load("_legacy_genome_utils", "_legacy_utils.py")
    legacy_data = _load("_legacy_genome_data", "_legacy_data.py")
    legacy_io = _load("_legacy_genome_io", "_legacy_io.py")


# --- data --------------------------------------------------------------------

GENE = "ATG" + "GCT" * 98 + "TAA"  # 300 nt, translates to 99 residues + stop


def make_record(rng, name, topology, ngenes):
    """One replicon of `ngenes` genes with the feature mix NCBI files have.

    Each gene has a CDS and a gene feature; about one CDS in 15 is a
    pseudogene, one in 10 a tRNA instead, and strands are random. Circular
    replicons may get a CDS that runs through the origin (a two-part
    location). Taxonomy and an Assembly dbxref are set like RefSeq entries.
    """
    seq = Seq(GENE * ngenes)
    rec = SeqRecord(seq, id=name, name=name, description=f"{name} toy replicon")
    rec.annotations = {
        "molecule_type": "DNA", "topology": topology, "organism": "Escherichia coli",
        "taxonomy": ["Bacteria", "Pseudomonadota", "Gammaproteobacteria"],
        "data_file_division": "BCT", "date": "01-JAN-2026",
    }
    rec.dbxrefs = ["Assembly:GCF_000005845.2"]
    rec.features.append(SeqFeature(FeatureLocation(0, len(seq), 1), type="source",
                                   qualifiers={"db_xref": ["taxon:562"]}))
    for i in range(ngenes):
        strand = int(rng.choice([1, -1]))
        loc = FeatureLocation(i * 300, i * 300 + 300, strand)
        tag = f"{name}_{i:04d}"
        draw = rng.random()
        if draw < 0.1:
            rec.features.append(SeqFeature(loc, type="tRNA", qualifiers={"locus_tag": [tag], "product": ["tRNA-Ala"]}))
        else:
            q = {"locus_tag": [tag], "protein_id": [f"WP_{name}_{i}.1"], "product": [f"protein {i}"]}
            if draw > 0.93:
                q = {"locus_tag": [tag], "pseudo": [""]}
            rec.features.append(SeqFeature(loc, type="CDS", qualifiers=q))
        rec.features.append(SeqFeature(loc, type="gene", qualifiers={"locus_tag": [tag], "gene": [f"g{i}"]}))
    if topology == "circular" and rng.random() < 0.5:
        end = len(seq)
        loc = CompoundLocation([FeatureLocation(end - 150, end, 1), FeatureLocation(0, 150, 1)])
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers={"locus_tag": [f"{name}_ori"], "protein_id": [f"WP_{name}_ori.1"]}))
    return rec


def make_records(seed):
    """1 to 3 replicons of 8 to 40 genes, some circular, from one seed."""
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 4))
    return [make_record(rng, f"NC_{seed:03d}{k}", rng.choice(["circular", "linear"]), int(rng.integers(8, 41)))
            for k in range(n)]


def with_queries(df, seed, adjacent=False):
    """Mark 1 to 3 CDS rows per replicon as queries (query == 1).

    With adjacent=True every chosen CDS also gets its next CDS marked, so
    windows overlap. Picks include the first and last CDS often enough to
    exercise the circular wrap.
    """
    rng = random.Random(seed)
    df = df.copy()
    df["query"] = 0
    for _, group in df[df.type == "CDS"].groupby("nucleotide", sort=True):
        idx = list(group.index)
        picks = set(rng.sample(idx, min(len(idx), rng.randint(1, 3))))
        if rng.random() < 0.3:
            picks.add(idx[0] if rng.random() < 0.5 else idx[-1])
        if adjacent:
            picks |= {idx[min(idx.index(p) + 1, len(idx) - 1)] for p in list(picks)}
        df.loc[sorted(picks), "query"] = 1
    return df


def legacy_frame(df):
    """The same table as a frozen-module NeighborhoodDF, so legacy methods run."""
    return legacy_data.NeighborhoodDF(pd.DataFrame(df))


SEEDS = range(24)


# --- comparison --------------------------------------------------------------

def _canonical_bytes(obj):
    """Encode every caller-visible property of a result as bytes.

    A frame gives its class name, column names with dtype, string storage
    and NA value, its index class, dtype and name, then one line per row
    with the label and the type and repr of every cell, so '', nan, None
    and <NA> never collide. A Series is encoded as a one-column frame.
    """
    if isinstance(obj, pd.Series):
        obj = obj.to_frame()
    parts = [type(obj).__qualname__]
    for name, dtype in obj.dtypes.items():
        parts.append(f"{name!r}\tdtype={dtype!r} storage={getattr(dtype, 'storage', None)!r} "
                     f"na_value={getattr(dtype, 'na_value', None)!r}")
    index = obj.index
    parts.append(f"index={type(index).__qualname__} dtype={index.dtype!r} name={index.name!r}")
    for label, row in zip(index, obj.itertuples(index=False, name=None)):
        cells = "\t".join(f"{type(v).__qualname__}:{v!r}" for v in row)
        parts.append(f"{label!r}\t{cells}")
    return "\n".join(parts).encode("utf-8")


def assert_same_frame(new, old):
    """Fail unless `new` and `old` show the same data down to the bytes."""
    assert type(new).__qualname__ == type(old).__qualname__
    pd.testing.assert_frame_equal(
        pd.DataFrame(new), pd.DataFrame(old), check_dtype=True, check_index_type=True,
        check_column_type=True, check_frame_type=True, check_names=True, check_exact=True,
        check_flags=True, check_freq=True,
    )
    assert _canonical_bytes(new) == _canonical_bytes(old)


def check_golden(name, data):
    """Compare `data` with data/genome/<name> as raw bytes (or write it)."""
    path = GOLDEN / name
    if WRITE_GOLDEN:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    assert path.read_bytes() == data


# --- seqrecords_to_dataframe and io.parse --------------------------------------

S2D_OPTIONS = [
    {},
    {"exclude_type": ["gene"]},
    {"autopid": "auto"},
    {"autopid": "copy", "exclude_type": ["tRNA"]},
    {"assembly": "GCA_999.1"},
]


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("options", S2D_OPTIONS, ids=str)
def test_seqrecords_to_dataframe_matches_legacy(seed, options):
    """The feature table is identical for every option set and replicon mix."""
    new = current_utils.seqrecords_to_dataframe(make_records(seed), **options)
    old = legacy_utils.seqrecords_to_dataframe(make_records(seed), **options)
    assert_same_frame(new, old)


def test_seqrecords_to_dataframe_empty_matches_legacy():
    """No records, or records without features, give the same empty table."""
    for arg in ([], None):
        assert_same_frame(current_utils.seqrecords_to_dataframe(arg), legacy_utils.seqrecords_to_dataframe(arg))


def _record_bytes(records):
    """Identifier, name, description and dbxrefs of each record, as bytes."""
    return "\n".join(f"{r.id!r}\t{r.name!r}\t{r.description!r}\t{r.dbxrefs!r}" for r in records).encode()


ASSEMBLY_RULES = {
    "default": {},
    "string": {"assembly": "GCA_000001.1"},
    "function": {"assembly": lambda path, rec: rec.id + "_asm"},
    "none": {"assembly": None},
}


@pytest.fixture(scope="module")
def genbank_file(tmp_path_factory):
    """Three seeds' replicons written to an NCBI-style .gbff file name."""
    path = tmp_path_factory.mktemp("gb") / "GCF_000005845.2_ASM584v2_genomic.gbff"
    SeqIO.write(make_records(1) + make_records(2) + make_records(3), path, "genbank")
    return path


@pytest.mark.parametrize("rule", ASSEMBLY_RULES)
def test_parse_matches_legacy(genbank_file, rule):
    """io.parse renames records and sets Assembly dbxrefs as before."""
    kwargs = ASSEMBLY_RULES[rule]
    new = list(current_io.parse(str(genbank_file), "genbank", **kwargs))
    old = list(legacy_io.parse(str(genbank_file), "genbank", **kwargs))
    assert _record_bytes(new) == _record_bytes(old)
    assert_same_frame(current_utils.seqrecords_to_dataframe(new), legacy_utils.seqrecords_to_dataframe(old))


def test_parse_golden(genbank_file):
    """The parsed records and their table match the stored golden bytes."""
    records = list(current_io.parse(str(genbank_file), "genbank"))
    check_golden("parse_records.txt", _record_bytes(records))
    check_golden("parse_table.txt", _canonical_bytes(current_utils.seqrecords_to_dataframe(records)))


# --- neighbors and vicinity ----------------------------------------------------

NEIGHBOR_OPTIONS = [
    dict(before=b, after=a, min_block_distance=d, strand=s)
    for (b, a), d, s in itertools.product([(0, 0), (1, 1), (4, 3)], [-1, 3], [None, "same"])
]


def _table(seed, adjacent=False):
    """The generated genome table for `seed`, with queries marked."""
    return with_queries(current_utils.seqrecords_to_dataframe(make_records(seed)), seed, adjacent)


@pytest.mark.parametrize("seed", SEEDS[:8])
@pytest.mark.parametrize("adjacent", [False, True])
def test_neighbors_matches_legacy(seed, adjacent):
    """neighbors() gives identical tables for every option set the cursors use."""
    df = _table(seed, adjacent)
    old_df = legacy_frame(df)
    for options in NEIGHBOR_OPTIONS:
        new = df.neighbors(df["query"] == 1, **options)
        old = old_df.neighbors(old_df["query"] == 1, **options)
        assert_same_frame(new, old)


@pytest.mark.parametrize("seed", SEEDS[:8])
def test_neighbors_rule_targets_match_legacy(seed):
    """String rules, lists of rules and pid selections select the same targets."""
    df = _table(seed)
    old_df = legacy_frame(df)
    pids = df.pid.dropna().unique()[:3].tolist()
    for targets in (["query == 1"], "query == 1", ["query == 1", "strand == 1"], df.pid.isin(pids)):
        old_targets = old_df.pid.isin(pids) if isinstance(targets, pd.Series) else targets
        assert_same_frame(df.neighbors(targets, before=2, after=2), old_df.neighbors(old_targets, before=2, after=2))


def test_neighbors_without_targets_matches_legacy():
    """No matching target gives the same empty table."""
    df = _table(0)
    assert_same_frame(df.neighbors(["pid == 'nothing'"]), legacy_frame(df).neighbors(["pid == 'nothing'"]))


@pytest.mark.parametrize("seed", SEEDS)
def test_vicinity_matches_legacy(seed):
    """vicinity() gives the same region table, including origin handling."""
    df = _table(seed, adjacent=seed % 2 == 1)
    old_df = legacy_frame(df)
    for b, a, d in [(0, 0, 0), (2, 2, 0), (3, 1, 2), (4, 4, -1)]:
        assert_same_frame(df.vicinity(df["query"] == 1, b, a, d), old_df.vicinity(old_df["query"] == 1, b, a, d))


def test_neighbors_golden():
    """neighbors() on fixed genomes matches the stored golden bytes."""
    parts = []
    for seed in (5, 11):
        for adjacent in (False, True):
            df = _table(seed, adjacent)
            for options in NEIGHBOR_OPTIONS[::3]:
                parts.append(f"# seed={seed} adjacent={adjacent} {options}".encode())
                parts.append(_canonical_bytes(df.neighbors(df["query"] == 1, **options)))
    check_golden("neighbors.txt", b"\n".join(parts))


# --- helpers -------------------------------------------------------------------

@pytest.mark.parametrize("seed", SEEDS[:12])
def test_boundaries_and_strand_helpers_match_legacy(seed):
    """boundaries, group_strand and same_strand give the same tables."""
    df = _table(seed)
    old_df = legacy_frame(df)
    assert_same_frame(df.boundaries(), old_df.boundaries())
    assert_same_frame(df.group_strand(), old_df.group_strand())
    assert_same_frame(df.same_strand(), old_df.same_strand())


# --- select_neighbors ----------------------------------------------------------

SELECT_OPTIONS = [
    dict(before=b, after=a, strand=s, new_query=nq)
    for (b, a), s, nq in itertools.product([(0, 0), (1, 1), (2, 1), (3, 3)], [False, True], [True, False])
]


def _no_overlap(df, before, after):
    """True when no two query windows of the same block share a row."""
    rows = np.flatnonzero(df["query"].to_numpy() == 1)
    blocks = df["block_id"].to_numpy()
    for i, j in zip(rows, rows[1:]):
        if blocks[i] == blocks[j] and j - after <= i + before:
            return False
    return True


@pytest.mark.parametrize("seed", SEEDS)
def test_select_neighbors_matches_legacy(seed):
    """select_neighbors keeps its old output where the old code was right.

    Only option sets whose windows do not overlap are compared; overlapping
    windows duplicated rows before and are pinned in test_trim_blocks.py.
    """
    df = _table(seed).reset_index(drop=True)
    old_df = legacy_frame(df)
    compared = 0
    for options in SELECT_OPTIONS:
        if not _no_overlap(df, options["before"], options["after"]):
            continue
        assert_same_frame(df.select_neighbors(**options), old_df.select_neighbors(**options))
        compared += 1
    assert compared > 0


@pytest.mark.parametrize("seed", SEEDS[:6])
def test_select_neighbors_column_value_matches_legacy(seed):
    """column/value selection picks the same rows as before."""
    df = _table(seed).reset_index(drop=True)
    pid = df.loc[df.type == "CDS", "pid"].dropna().iloc[len(df) % 5]
    assert_same_frame(df.select_neighbors(column="pid", value=pid, before=2, after=2),
                      legacy_frame(df).select_neighbors(column="pid", value=pid, before=2, after=2))
