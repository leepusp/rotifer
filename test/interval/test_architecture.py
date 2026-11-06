"""
Tests for rotifer.interval.architecture and the domain2architecture and
architecture2table commands, formerly Perl scripts. The golden cases in
test/data/golden hold the Perl output on recorded HMMER, RPS-BLAST and
Phobius tables, including nested architectures.
"""

import io
import os

import pytest

import golden
from rotifer.interval import architecture as arch
from rotifer.io.hmmer import read_hmmer
from rotifer.io.phobius import read_phobius

@pytest.mark.parametrize("case,cmdline", golden.cases("domain2architecture"))
def test_domain2architecture_matches_perl_output(case, cmdline, tmp_path):
    golden.check("domain2architecture", case, cmdline, tmp_path)

@pytest.mark.parametrize("case,cmdline", golden.cases("architecture2table"))
def test_architecture2table_matches_perl_output(case, cmdline, tmp_path):
    golden.check("architecture2table", case, cmdline, tmp_path)

def _table(text):
    import pandas as pd
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)

def test_domain2architecture_from_files():
    df = arch.domain2architecture([os.path.join(golden.INPUTS, "hmmscan.version.tsv"),
                                   os.path.join(golden.INPUTS, "phobius.tsv")], evalue=0.0101)
    assert df.equals(_table(golden.expected("domain2architecture", "pfam_phobius")))

def test_domain2architecture_from_dataframes():
    """The library chains with the parsers, without files or subprocesses."""
    import pandas as pd
    regions = pd.concat([
        read_hmmer(os.path.join(golden.INPUTS, "hmmscan.out"), clean={"model": "version"}),
        read_phobius(os.path.join(golden.INPUTS, "phobius.out"), evalue="0.0101"),
    ])
    df = arch.domain2architecture(regions, evalue=0.0101)
    expected = _table(golden.expected("domain2architecture", "pfam_phobius"))
    assert list(df.ID) == list(expected.ID)
    assert list(df.architecture) == list(expected.architecture)
    assert list(df.regions) == list(expected.regions)

def test_architecture2table():
    df = arch.architecture2table(os.path.join(golden.INPUTS, "pfam.arch"))
    expected = _table(golden.expected("architecture2table", "pfam"))
    assert list(df.columns) == arch.TABLE_COLUMNS
    assert list(df.domain) == list(expected.domain)
    assert list(df.start) == [int(x) for x in expected.start]
    assert arch.architecture2table(arch.domain2architecture(
        os.path.join(golden.INPUTS, "hmmscan.version.tsv"), evalue=0.0101)).shape[1] == 5

def test_library_chain_replaces_shell_pipeline():
    """
    Devel code used to run 'hmmer2table | domain2architecture |
    architecture2table' and read the result with pandas.read_csv; the
    library chain must return the same table.
    """
    import pandas as pd
    chained = arch.architecture2table(arch.domain2architecture(
        [os.path.join(golden.INPUTS, "hmmscan.version.tsv"), os.path.join(golden.INPUTS, "phobius.tsv")],
        evalue=0.0101))
    expected = pd.read_csv(io.StringIO(golden.expected("architecture2table", "pfam")), sep="\t")
    pd.testing.assert_frame_equal(chained, expected)

@pytest.mark.parametrize("overlap,expected", [
    ("0.4", "A"),    # 10 residues overlap > 40% of B (20 residues): B discarded
    ("10", "A+B"),   # 10 residues are tolerated
])
def test_overlap_modes(overlap, expected):
    regions = [{"domain": "A", "start": "1", "end": "100", "length": 100, "evalue": "1e-10"},
               {"domain": "B", "start": "91", "end": "110", "length": 20, "evalue": "1e-5"}]
    assert arch.build_architecture(regions, overlap=overlap)[0] == expected

def test_nesting():
    regions = [{"domain": "A", "start": "1", "end": "100", "length": 100, "evalue": "1e-10"},
               {"domain": "B", "start": "20", "end": "40", "length": 21, "evalue": "1e-5"}]
    assert arch.build_architecture(regions, ignore_overlap={"all"})[0] == "A[B]"
