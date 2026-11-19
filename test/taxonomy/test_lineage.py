#!/usr/bin/env python3
"""Exact values for what rotifer.taxonomy.utils.lineage() now does differently.

What has to hold, and why:

  * a missing classification (None, NaN, pd.NA) gives a missing lineage: NaN
    in a Series of any string dtype, None in a list. Under pandas 3 the old
    code raised TypeError on it; under pandas 2 it gave ''
  * duplicate index labels keep every row, in order; the old code merged
    the rows that shared a label
  * an unsorted unique index keeps the input order; the result is the old
    one reindexed to the input index, byte for byte
  * `case` applies to preferred_taxa too, and names are stripped, on the
    input side after splitting and on the preferred side before matching
  * preferred_taxa=None reads the shipped list; an empty iterable keeps
    nothing, where the old code used the shipped list for it
  * a separator of more than one character is split on literally; the old
    Series path read it as a regular expression ('..' matched any two
    characters)
  * a str is rejected, data has no default, a missing column is a KeyError
  * a generator gives a list, and the input is never modified
  * load_preferred_taxa() reads every file found, skips comments and blank
    lines, strips and lowercases; loadPreferredTaxa() warns at the caller's
    line and returns the same set

The unchanged behaviour is pinned byte for byte in test_lineage_equivalence.py.

Run under pytest: PYTHONPATH=lib pytest test/taxonomy/test_lineage.py
"""

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest

from rotifer.taxonomy import utils as current
from rotifer.taxonomy.utils import lineage, load_preferred_taxa, loadPreferredTaxa

HERE = pathlib.Path(__file__).resolve().parent


def _load(name, filename):
    """Import a sibling test module by path; test/ is not a package."""
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = _load("_legacy_taxonomy_utils", "_legacy_utils.py")
harness = _load("_lineage_equivalence_harness", "test_lineage_equivalence.py")

ECOLI = "Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Escherichia"
ECOLI_LINEAGE = "bacteria>pseudomonadota>gammaproteobacteria>enterobacterales"
ARCHAEA = "Archaea; Methanobacteriota; Methanobacteria"
ARCHAEA_LINEAGE = "archaea>methanobacteria"


def check(result, expected, **kwargs):
    """Compare two Series with every strict check on, except what `kwargs` relaxes."""
    options = dict(check_dtype=True, check_index_type=True, check_exact=True, check_names=True)
    options.update(kwargs)
    pd.testing.assert_series_equal(result, expected, **options)


def expected_series(values, index=None, name=None):
    """A `str` Series built the way lineage() should return it."""
    return pd.Series(values, index=index, name=name, dtype="str")


# --- missing values ----------------------------------------------------------

@pytest.mark.parametrize("missing", [None, np.nan, pd.NA], ids=["None", "nan", "NA"])
def test_missing_in_object_series(missing):
    """A missing value in an object Series gives NaN, and the neighbours are kept."""
    data = pd.Series([ECOLI, missing, ARCHAEA], dtype=object)
    check(lineage(data), expected_series([ECOLI_LINEAGE, np.nan, ARCHAEA_LINEAGE], index=pd.Index([0, 1, 2])))


@pytest.mark.parametrize("missing", [None, np.nan, pd.NA], ids=["None", "nan", "NA"])
def test_missing_in_str_series(missing):
    """A missing value in a `str` dtype Series gives NaN."""
    data = pd.Series([ECOLI, missing, ARCHAEA], dtype="str")
    check(lineage(data), expected_series([ECOLI_LINEAGE, np.nan, ARCHAEA_LINEAGE], index=pd.Index([0, 1, 2])))


@pytest.mark.parametrize("missing", [None, np.nan, pd.NA], ids=["None", "nan", "NA"])
def test_missing_in_string_python_series(missing):
    """A missing value in a string[python] Series gives NaN in a `str` Series."""
    data = pd.Series([ECOLI, missing, ARCHAEA], dtype=pd.StringDtype("python"))
    check(lineage(data), expected_series([ECOLI_LINEAGE, np.nan, ARCHAEA_LINEAGE], index=pd.Index([0, 1, 2])))


def test_all_missing_series():
    """A Series with nothing but missing values gives NaN everywhere, in `str` dtype."""
    result = lineage(pd.Series([None, np.nan, pd.NA], dtype=object))
    check(result, expected_series([np.nan] * 3, index=pd.Index([0, 1, 2])))


def test_missing_in_list():
    """Every kind of missing value in a list gives None, in place."""
    result = lineage([ECOLI, None, np.nan, pd.NA, ARCHAEA, ""])
    assert result == [ECOLI_LINEAGE, None, None, None, ARCHAEA_LINEAGE, ""]
    assert [type(x) for x in result] == [str, type(None), type(None), type(None), str, str]


def test_missing_is_not_the_empty_lineage():
    """An empty classification gives '' while a missing one gives NaN."""
    result = lineage(pd.Series(["", None]))
    assert result.iloc[0] == ""
    assert pd.isna(result.iloc[1])


def test_non_string_values_go_through_str():
    """Numbers are read as their str, so '1' and 1 give the same lineage."""
    assert lineage([1, "1"], preferred_taxa=["1"]) == ["1", "1"]
    check(lineage(pd.Series([1, 2], dtype=object), preferred_taxa=["2"]),
          expected_series(["", "2"], index=pd.Index([0, 1])))


# --- the index ---------------------------------------------------------------

def test_duplicate_index_keeps_every_row():
    """Rows that share an index label are all kept, in input order."""
    data = pd.Series([ECOLI, ARCHAEA, "Bacteria; Bacilli", ""], index=["a", "b", "a", "b"], name="c")
    expected = expected_series([ECOLI_LINEAGE, ARCHAEA_LINEAGE, "bacteria>bacilli", ""],
                               index=pd.Index(["a", "b", "a", "b"]), name="c")
    check(lineage(data), expected)
    assert len(lineage(data)) == 4


def test_duplicate_integer_index_with_missing():
    """Duplicate integer labels and a missing value together keep all rows."""
    data = pd.Series([ARCHAEA, None, ECOLI], index=[3, 3, 3])
    check(lineage(data), expected_series([ARCHAEA_LINEAGE, np.nan, ECOLI_LINEAGE], index=pd.Index([3, 3, 3])))


def test_duplicate_index_in_dataframe():
    """A DataFrame with duplicate labels keeps all rows, named after `column`."""
    frame = pd.DataFrame({"cls": [ARCHAEA, ECOLI]}, index=pd.Index([7, 7], name="pid"))
    check(lineage(frame, column="cls"),
          expected_series([ARCHAEA_LINEAGE, ECOLI_LINEAGE], index=pd.Index([7, 7], name="pid"), name="cls"))


@pytest.mark.parametrize("index", [[30, 10, 20], ["z", "a", "m"]], ids=["int", "str"])
def test_unsorted_unique_index_equals_legacy_reindexed(index):
    """With a unique unsorted index the result is the old one reindexed, down to the bytes."""
    data = pd.Series([ECOLI, ARCHAEA, "Bacteria; Bacilli"], index=pd.Index(index, name="k"), name="c")
    old = legacy.lineage(data)
    assert list(old.index) == sorted(index)
    harness.assert_same_bytes(lineage(data), old.reindex(data.index))
    assert list(lineage(data).index) == index


def test_unsorted_unique_index_keeps_input_order():
    """Rows stay where they were, whatever the order of the labels."""
    data = pd.Series([ARCHAEA, ECOLI], index=[9, 1])
    assert lineage(data).tolist() == [ARCHAEA_LINEAGE, ECOLI_LINEAGE]


def test_rangeindex_comes_back_as_int64_index():
    """A RangeIndex is returned as a plain int64 Index, as before."""
    result = lineage(pd.Series([ECOLI, ARCHAEA]))
    assert type(result.index) is pd.Index
    assert result.index.dtype == np.dtype("int64")


# --- preferred taxa ----------------------------------------------------------

def test_preferred_taxa_follow_case():
    """Mixed-case preferred names match, and the result is lowercase."""
    result = lineage([ECOLI, ARCHAEA], preferred_taxa=["Bacteria", "ENTEROBACTERALES", "Methanobacteria"])
    assert result == ["bacteria>enterobacterales", "methanobacteria"]


def test_preferred_taxa_in_a_set_and_generator():
    """preferred_taxa may be any iterable, a set or a generator included."""
    assert lineage([ECOLI], preferred_taxa={"Bacteria"}) == ["bacteria"]
    assert lineage([ECOLI], preferred_taxa=(x for x in ["Bacteria"])) == ["bacteria"]


def test_case_none_matches_exactly():
    """With case=None names are compared and returned as written."""
    assert lineage([ECOLI], preferred_taxa=["Bacteria", "Escherichia"], case=None) == ["Bacteria>Escherichia"]
    assert lineage([ECOLI], preferred_taxa=["bacteria"], case=None) == [""]
    check(lineage(pd.Series([ECOLI]), preferred_taxa=["Bacteria"], case=None),
          expected_series(["Bacteria"], index=pd.Index([0])))


def test_case_upper_applies_to_both_sides():
    """A custom case function is applied to the data and to preferred_taxa."""
    assert lineage([ECOLI], preferred_taxa=["bacteria"], case=str.upper) == ["BACTERIA"]


@pytest.mark.parametrize("empty", [[], set(), (), iter(())], ids=["list", "set", "tuple", "iterator"])
def test_empty_preferred_taxa_keep_nothing(empty):
    """An empty iterable keeps nothing; it no longer means the default list."""
    assert lineage([ECOLI, ARCHAEA], preferred_taxa=empty) == ["", ""]
    check(lineage(pd.Series([ECOLI]), preferred_taxa=empty if not hasattr(empty, "__next__") else []),
          expected_series([""], index=pd.Index([0])))


def test_preferred_taxa_none_reads_the_shipped_list():
    """preferred_taxa=None gives the same result as passing the shipped list."""
    default = load_preferred_taxa()
    assert lineage([ECOLI], preferred_taxa=None) == [ECOLI_LINEAGE]
    assert lineage([ECOLI, ARCHAEA]) == lineage([ECOLI, ARCHAEA], preferred_taxa=default)


def test_preferred_taxa_are_stripped():
    """Spaces around preferred names do not stop a match."""
    assert lineage([ECOLI], preferred_taxa=["  Bacteria ", "\tEscherichia\n"]) == ["bacteria>escherichia"]


# --- whitespace and separators -----------------------------------------------

def test_names_are_stripped_after_splitting():
    """Spaces around names do not stop a match, whatever the separator."""
    assert lineage(["Bacteria ;  Bacilli "], preferred_taxa=["bacteria", "bacilli"], insep=";") == ["bacteria>bacilli"]
    check(lineage(pd.Series([" Bacteria ;\tBacilli"]), preferred_taxa=["bacteria", "bacilli"], insep=";"),
          expected_series(["bacteria>bacilli"], index=pd.Index([0])))


def test_outsep_and_insep():
    """Custom separators are used on the way in and on the way out."""
    assert lineage(["Bacteria|Bacilli"], preferred_taxa=["bacteria", "bacilli"], insep="|", outsep=" / ") == [
        "bacteria / bacilli"]


def test_multicharacter_insep_is_literal():
    """insep='..' splits on two dots; the old Series path read it as a regex and kept nothing."""
    data = pd.Series(["Bacteria..Pseudomonadota"])
    check(lineage(data, insep=".."), expected_series(["bacteria>pseudomonadota"], index=pd.Index([0])))
    assert list(legacy.lineage(data, insep="..")) == [""]


def test_order_is_the_original_order():
    """Kept names come out in the order of the classification, not of preferred_taxa."""
    assert lineage([ECOLI], preferred_taxa=["Enterobacterales", "Bacteria"]) == ["bacteria>enterobacterales"]


# --- errors ------------------------------------------------------------------

def test_str_is_rejected():
    """A single str raises TypeError instead of being read character by character."""
    with pytest.raises(TypeError):
        lineage(ECOLI)


def test_data_is_required():
    """Calling lineage() with no data raises TypeError."""
    with pytest.raises(TypeError):
        lineage()


def test_missing_column_is_a_keyerror():
    """A DataFrame without `column` raises KeyError."""
    with pytest.raises(KeyError):
        lineage(pd.DataFrame({"other": [ECOLI]}))
    with pytest.raises(KeyError):
        lineage(pd.DataFrame({"classification": [ECOLI]}), column="cls")


# --- containers --------------------------------------------------------------

def test_generator_gives_a_list():
    """A generator is read item by item and gives a list."""
    result = lineage(x for x in [ECOLI, None, ARCHAEA])
    assert type(result) is list
    assert result == [ECOLI_LINEAGE, None, ARCHAEA_LINEAGE]


def test_tuple_and_list_give_lists():
    """A tuple and a list give a list."""
    assert lineage((ECOLI,)) == [ECOLI_LINEAGE]
    assert type(lineage((ECOLI,))) is list


def test_empty_inputs():
    """Empty inputs give empty results of the old types."""
    assert lineage([]) == []
    result = lineage(pd.Series([], dtype=str))
    assert len(result) == 0 and result.dtype == object


def test_dataframe_gives_a_series_named_after_column():
    """A DataFrame gives a Series named after `column`, with the frame's index."""
    frame = pd.DataFrame({"classification": [ECOLI, None]}, index=pd.Index([5, 4], name="pid"))
    check(lineage(frame), expected_series([ECOLI_LINEAGE, np.nan], index=pd.Index([5, 4], name="pid"),
                                          name="classification"))


def test_series_input_is_not_modified():
    """The input Series, its index and its values are untouched."""
    data = pd.Series([ECOLI, None, ARCHAEA], index=[3, 3, 1], name="c")
    before = data.copy(deep=True)
    lineage(data)
    pd.testing.assert_series_equal(data, before, check_exact=True)
    assert data.index.equals(before.index)


def test_dataframe_input_is_not_modified():
    """The input DataFrame is untouched."""
    frame = pd.DataFrame({"classification": [ECOLI, None], "x": [1, 2]})
    before = frame.copy(deep=True)
    lineage(frame)
    pd.testing.assert_frame_equal(frame, before, check_exact=True)


# --- load_preferred_taxa -----------------------------------------------------

def test_load_preferred_taxa_from_a_file(tmp_path):
    """Comments, blank lines, trailing spaces, CRLF and mixed case are handled."""
    path = tmp_path / "taxa.txt"
    path.write_bytes(b"# a comment\r\nBacteria\r\n\r\n  Archaea  \n   \nPseudomonadota \t\n"
                     b"#Euryarchaeota\nBACTERIA\nBacilli")
    assert load_preferred_taxa(str(path)) == {"bacteria", "archaea", "pseudomonadota", "bacilli"}


def test_load_preferred_taxa_empty_file(tmp_path):
    """A file with only comments and blanks gives an empty set."""
    path = tmp_path / "taxa.txt"
    path.write_text("# nothing\n\n")
    assert load_preferred_taxa(str(path)) == set()


def test_load_preferred_taxa_default():
    """The default list is a set of lowercase names with no blanks or comments."""
    taxa = load_preferred_taxa()
    assert type(taxa) is set
    assert {"bacteria", "archaea", "eukaryota", "enterobacterales"} <= taxa
    assert all(name == name.strip().lower() and name and not name.startswith("#") for name in taxa)


def test_load_preferred_taxa_default_matches_legacy():
    """The new reader gives the old set for the shipped file."""
    assert load_preferred_taxa() == legacy.loadPreferredTaxa()


def test_loadpreferredtaxa_warns_at_the_caller():
    """The deprecated name warns with FutureWarning, pointing at the caller's file."""
    with pytest.warns(FutureWarning, match="load_preferred_taxa") as record:
        result = loadPreferredTaxa()
    assert len(record) == 1
    assert record[0].filename == __file__
    assert result == load_preferred_taxa()
    assert type(result) is set
