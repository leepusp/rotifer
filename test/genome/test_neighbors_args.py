#!/usr/bin/env python3
"""Tests for the argument checks of NeighborhoodDF.neighbors and vicinity.

What has to hold, and why:

  * strand accepts None, 'any' and 'same'; None and 'any' both mean "ignore
    strand" and give the same rows, because every cursor passes strand=None
    by default
  * fttype accepts only 'same'; the values the docstrings used to advertise
    ('+' and '-' for strand, 'any' for fttype) never worked, so they raise
    ValueError instead of being ignored or returning an empty table. The
    message names the argument, the value and the allowed values
  * a string target rule that cannot be evaluated raises ValueError naming
    the rule, chained from the original error, instead of being logged and
    skipped
  * a target that is neither a string nor a Series raises TypeError
  * when no row matches, the result is still the empty table the cursors
    rely on (neighbors) or an empty frame (vicinity), not an error

The data are one hand-built linear replicon of 10 CDS (strands + + - - + + -
+ + +, internal_id 0 to 9) with query 1 on rows 3 and 4; nothing touches the
network.

Run under pytest: PYTHONPATH=lib pytest test/genome/test_neighbors_args.py
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from rotifer.genome.utils import seqrecords_to_dataframe


def make_table(strands=(1, 1, -1, -1, 1, 1, -1, 1, 1, 1), queries=(3, 4)):
    """One linear replicon of 300 nt CDS with the given strands.

    Rows listed in `queries` get query == 1, the others 0. Returns the
    NeighborhoodDF that seqrecords_to_dataframe builds from it.
    """
    rec = SeqRecord(Seq("A" * 300 * len(strands)), id="N1", name="N1", description="toy")
    rec.annotations = {"molecule_type": "DNA", "topology": "linear",
                       "organism": "E", "taxonomy": ["B"]}
    rec.dbxrefs = ["Assembly:GCF_1.1"]
    for i, s in enumerate(strands):
        rec.features.append(SeqFeature(
            FeatureLocation(i * 300, i * 300 + 300, s), type="CDS",
            qualifiers={"locus_tag": [f"N1_{i}"], "protein_id": [f"P_N1_{i}.1"]}))
    df = seqrecords_to_dataframe([rec])
    df["query"] = 0
    df.loc[list(queries), "query"] = 1
    return df


@pytest.fixture
def table():
    """The 10-CDS replicon with queries on rows 3 and 4."""
    return make_table()


@pytest.mark.parametrize("strand", ["+", "-", "both", 1, True])
def test_bad_strand_raises(table, strand):
    """Only None, 'any' and 'same' are strand values; the rest raises ValueError."""
    with pytest.raises(ValueError, match=r"strand=.* not supported.*None.*'any'.*'same'"):
        table.neighbors("query == 1", strand=strand)


@pytest.mark.parametrize("method", ["neighbors", "vicinity"])
@pytest.mark.parametrize("fttype", ["any", "different", None, 1])
def test_bad_fttype_raises(table, method, fttype):
    """Only fttype='same' works; 'any' used to log an error and return nothing."""
    with pytest.raises(ValueError, match=r"fttype=.* not supported.*'same'"):
        getattr(table, method)("query == 1", fttype=fttype)


def test_error_message_names_argument_and_value(table):
    """The message holds the argument name, the value given and the allowed values."""
    with pytest.raises(ValueError) as err:
        table.neighbors("query == 1", strand="+")
    assert str(err.value) == "strand='+' is not supported; use one of None, 'any', 'same'"
    with pytest.raises(ValueError) as err:
        table.neighbors("query == 1", fttype="any")
    assert str(err.value) == "fttype='any' is not supported; use one of 'same'"


def test_strand_none_equals_any(table):
    """strand=None (the cursors' default) and strand='any' return the same table."""
    none = table.neighbors("query == 1", before=2, after=2, strand=None)
    any_ = table.neighbors("query == 1", before=2, after=2, strand="any")
    assert len(none) == 6
    pd.testing.assert_frame_equal(pd.DataFrame(none), pd.DataFrame(any_), check_exact=True)


def test_strand_same_keeps_runs_with_a_target(table):
    """strand='same' splits at strand changes and keeps runs that hold a target."""
    out = table.neighbors("query == 1", before=1, after=1, strand="same")
    assert out["pid"].tolist() == ["P_N1_2.1", "P_N1_3.1", "P_N1_4.1", "P_N1_5.1"]
    assert out["block_id"].tolist() == ["N1:601-1200"] * 2 + ["N1:1201-1800"] * 2
    assert out["query"].tolist() == [0, 1, 1, 0]


@pytest.mark.parametrize("method", ["neighbors", "vicinity"])
def test_bad_rule_raises_value_error_naming_the_rule(table, method):
    """A rule that eval cannot evaluate raises ValueError, chained from the cause."""
    with pytest.raises(ValueError, match="nonexistent_col == 1") as err:
        getattr(table, method)("nonexistent_col == 1")
    assert err.value.__cause__ is not None


def test_bad_rule_in_a_list_raises(table):
    """One failing rule among good ones still raises."""
    with pytest.raises(ValueError, match="syntax error here"):
        table.neighbors(["query == 1", "syntax error here"])


@pytest.mark.parametrize("method", ["neighbors", "vicinity"])
@pytest.mark.parametrize("target", [5, 1.5, None, ["query == 1", 3]])
def test_target_of_wrong_type_raises(table, method, target):
    """A target that is neither a string nor a Series raises TypeError."""
    with pytest.raises(TypeError, match="targets must be"):
        getattr(table, method)(target)


def test_no_match_still_returns_the_empty_table(table):
    """neighbors gives the same empty table as before when no row is a target."""
    out = table.neighbors("pid == 'nothing'")
    expected = seqrecords_to_dataframe([])
    assert len(out) == 0
    pd.testing.assert_frame_equal(pd.DataFrame(out), pd.DataFrame(expected), check_exact=True)
    assert list(out.columns) == list(expected.columns)


def test_vicinity_no_match_returns_an_empty_frame(table):
    """vicinity returns an empty frame when no row is a target."""
    assert table.vicinity("pid == 'nothing'").empty
