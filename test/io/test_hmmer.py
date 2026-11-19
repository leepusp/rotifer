"""
Tests for rotifer.io.hmmer and the hmmer2table command.

hmmer2table used to be a Perl script. The golden cases in
test/data/golden/hmmer2table hold its output on recorded hmmscan and
hmmsearch reports; the Python command must print the same bytes, and
read_hmmer must load the same rows.
"""

import io
import os

import pytest

import golden
from rotifer.io import hmmer

@pytest.mark.parametrize("case,cmdline", golden.cases("hmmer2table"))
def test_command_matches_perl_output(case, cmdline, tmp_path):
    golden.check("hmmer2table", case, cmdline, tmp_path)

def _golden_table(case, header=False):
    import pandas as pd
    text = golden.expected("hmmer2table", case)
    return pd.read_csv(io.StringIO(text), sep="\t", header=0 if header else None,
                       names=hmmer.COLUMNS, dtype=hmmer.DTYPES)

def test_read_hmmer_hmmscan():
    df = hmmer.read_hmmer(os.path.join(golden.INPUTS, "hmmscan.out"))
    assert list(df.columns) == hmmer.COLUMNS
    assert df.equals(_golden_table("hmmscan"))

def test_read_hmmer_hmmsearch_with_cutoff():
    df = hmmer.read_hmmer(os.path.join(golden.INPUTS, "hmmsearch.out"), evalue=1e-3)
    assert df.equals(_golden_table("hmmsearch_evalue", header=True))

def test_clean_accepts_functions_and_names():
    path = os.path.join(golden.INPUTS, "hmmscan.out")
    by_name = hmmer.read_hmmer(path, clean={"model": "version"})
    by_function = hmmer.read_hmmer(path, clean={"query": lambda name: name.split(".")[0]})
    assert list(by_name.model) == list(by_function.model)

def test_unknown_cleaner_is_rejected():
    with pytest.raises(ValueError):
        hmmer.read_hmmer(os.path.join(golden.INPUTS, "hmmscan.out"), clean={"model": "nope"})

def test_force_basename(tmp_path):
    report = tmp_path / "PF00072.hmmsearch"
    report.write_text(golden.read("inputs/hmmsearch.out"))
    df = hmmer.read_hmmer(str(report), force_basename=[r"\.hmmsearch"])
    assert set(df.model) == {"PF00072"}
