"""
Tests for rotifer.io.blast and the blast2table command.

blast2table used to be a Perl script built on BioPerl. The golden cases
in test/data/golden/blast2table hold its output on recorded RPS-BLAST,
PSI-BLAST and multi-query BLASTP reports; the Python command must print
the same bytes, including merged (redundant) HSPs and CDD name cleaning.
"""

import os

import pytest

import golden
from rotifer.io import blast

@pytest.mark.parametrize("case,cmdline", golden.cases("blast2table"))
def test_command_matches_perl_output(case, cmdline, tmp_path):
    golden.check("blast2table", case, cmdline, tmp_path)

def test_read_blast_swapped_cdd():
    df = blast.read_blast(os.path.join(golden.INPUTS, "rpsblast.out"), swap=True, clean=["cdd"])
    rows = [line.split("\t") for line in golden.expected("blast2table", "rps_swap_cdd").splitlines()]
    assert list(df.columns) == blast.COLUMNS
    assert len(df) == len(rows)
    assert list(df.hit_name) == [row[0] for row in rows]
    assert list(df.query_name) == [row[1] for row in rows]
    assert list(df.hit_start) == [int(row[2]) for row in rows]

def test_read_blast_replaces_command_plus_read_csv():
    """
    Devel code (psiblast functions) used to run blast2table and load its
    output with pandas.read_csv; read_blast must return the same table.
    """
    import io
    import pandas as pd
    df = blast.read_blast(os.path.join(golden.INPUTS, "psiblast.out"))
    expected = pd.read_csv(io.StringIO(golden.expected("blast2table", "psiblast")), sep="\t", names=blast.COLUMNS)
    pd.testing.assert_frame_equal(df, expected)

def test_psiblast_rounds_are_merged():
    df = blast.read_blast(os.path.join(golden.INPUTS, "psiblast.out"))
    assert set(df.iteration) == {1}
    assert not df.duplicated(["hit_name", "query_name", "hit_start", "hit_end"]).any()

def test_parse_blast_text():
    with open(os.path.join(golden.INPUTS, "psiblast.out")) as handle:
        results = list(blast.parse_blast_text(handle))
    assert len(results) == 2  # one result per PSI-BLAST round
    assert results[0]["query_name"] == "NP_414555.1"
    assert results[0]["query_length"] == 638
    first = results[0]["hits"][0]
    assert first["name"] == "NP_414555.1" and first["length"] == 638
    assert first["hsps"][0]["query_start"] == 1 and first["hsps"][0]["query_end"] == 638

def test_unknown_cleaner_is_rejected():
    with pytest.raises(ValueError):
        blast.read_blast(os.path.join(golden.INPUTS, "psiblast.out"), clean=["nope"])
