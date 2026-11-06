"""
Tests for rotifer.io.phobius and the phobius2table command.

phobius2table used to be a Perl script; the golden cases in
test/data/golden/phobius2table hold its output on a recorded Phobius
report, which the Python command must reproduce byte for byte.
"""

import os

import pytest

import golden
from rotifer.io import phobius

@pytest.mark.parametrize("case,cmdline", golden.cases("phobius2table"))
def test_command_matches_perl_output(case, cmdline, tmp_path):
    golden.check("phobius2table", case, cmdline, tmp_path)

def test_read_phobius():
    df = phobius.read_phobius(os.path.join(golden.INPUTS, "phobius.out"), evalue=0.0101)
    rows = [line.split("\t") for line in golden.expected("phobius2table", "evalue").splitlines()]
    assert len(df) == len(rows)
    assert list(df.columns) == phobius.COLUMNS
    assert df.iloc[0].tolist() == [rows[0][0], rows[0][1], int(rows[0][2]), int(rows[0][3]), 0.0101]
    assert set(df.domain) == {"SP", "TM"}

@pytest.mark.parametrize("identifier,expected", [
    ("NP_414542.1", "NP_414542.1"),
    ("NP_414542.1/10-200", "NP_414542.1"),
    ("gi|16127995|ref|NP_414542.1|", "16127995"),
    ("sp|P0AD86|THRL_ECOLI", "P0AD86"),
    ("|", ""),
])
def test_accession(identifier, expected):
    assert phobius.accession(identifier) == expected
