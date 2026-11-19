"""
Tests for rotifer.io.dalidat and the dat2fa command (formerly a Perl
script, used by sbin/update_pdb to build the DaliLite FASTA database).
"""

import os

import pytest

import golden
from rotifer.io import dalidat

@pytest.mark.parametrize("case,cmdline", golden.cases("dat2fa"))
def test_command_matches_perl_output(case, cmdline, tmp_path):
    golden.check("dat2fa", case, cmdline, tmp_path)

def test_parse_dat():
    with open(os.path.join(golden.INPUTS, "dali", "a0rrA.dat")) as handle:
        records = list(dalidat.parse_dat(handle))
    fasta = golden.expected("dat2fa", "files").splitlines()
    header = next(line for line in fasta if line.startswith(">a0rrA"))
    identifier, description, sequence = records[0]
    assert f">{identifier} {description}" == header
    assert sequence == fasta[fasta.index(header) + 1]
