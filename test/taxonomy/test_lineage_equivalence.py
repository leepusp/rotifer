#!/usr/bin/env python3
"""Byte-for-byte equivalence of rotifer.taxonomy.utils with the code it replaced.

What has to hold, and why:

  * on every input the old lineage() handled correctly, the current one
    returns the same object down to the bytes: same type, values, dtype,
    name, index class, index dtype and index name. NeighborhoodDF and the
    three taxonomy cursors write this result straight into their tables,
    so any difference, even one that compares equal, reaches users
  * "handled correctly" means no missing values, a unique index sorted
    ascending, preferred_taxa already in the case that `case` returns, and
    no spaces around the input separator. Inputs outside that domain were
    buggy and are pinned with exact values in test_lineage.py instead
  * the default list of preferred taxa read from share/rotifer/taxonomy
    stays the same set

Three checks compare two results, strictest last:

  1. pandas.testing.assert_series_equal with every check on (Series only)
  2. _canonical_bytes, a readable encoding of everything a caller can see,
     so a failure shows which part differs
  3. pickle.dumps of both objects, which must be identical bytes

The reference is _legacy_utils.py, a verbatim copy of the old module, and
data/taxonomy/lineage_golden.tsv, the old module's output on 2,000 real
NCBI classifications. The golden file is compared as raw bytes, so it keeps
guarding the output after the legacy copy is deleted.

Do not edit these tests to make a change pass. A change that breaks them
changes user-visible output and needs its own test in test_lineage.py and
an entry in the PR's "Output that changes" list.

Run under pytest: PYTHONPATH=lib pytest test/taxonomy/test_lineage_equivalence.py
"""

import importlib.util
import pathlib
import pickle
import random
import warnings

import pandas as pd
import pytest

from rotifer.taxonomy import utils as current

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "taxonomy"
FIXTURE = DATA / "ncbi_classifications.tsv"
GOLDEN = DATA / "lineage_golden.tsv"


def _load_legacy():
    """Import _legacy_utils.py by path; test/ is not a package."""
    spec = importlib.util.spec_from_file_location("_legacy_taxonomy_utils", HERE / "_legacy_utils.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = _load_legacy()


def _legacy_taxa():
    """The old default set of preferred taxa (lowercase clade names)."""
    return legacy.loadPreferredTaxa()


PREFERRED = sorted(_legacy_taxa())


# --- comparison --------------------------------------------------------------

def _canonical_bytes(obj):
    """Encode every caller-visible property of a lineage() result as bytes.

    A list becomes its type plus the type and repr of each item. A Series
    adds its dtype (with string storage and NA value), its name, its index
    class, dtype and name, and one 'label<TAB>value' line per row, both as
    repr, so '', nan, None and <NA> never collide.
    """
    parts = [type(obj).__qualname__]
    if isinstance(obj, pd.Series):
        dtype = obj.dtype
        parts.append(f"dtype={dtype!r} storage={getattr(dtype, 'storage', None)!r} "
                     f"na_value={getattr(dtype, 'na_value', None)!r}")
        parts.append(f"name={obj.name!r}")
        index = obj.index
        parts.append(f"index={type(index).__qualname__} dtype={index.dtype!r} name={index.name!r}")
        parts.extend(f"{label!r}\t{type(value).__qualname__}\t{value!r}" for label, value in obj.items())
    else:
        parts.extend(f"{type(value).__qualname__}\t{value!r}" for value in obj)
    return "\n".join(parts).encode("utf-8")


def assert_same_bytes(new, old):
    """Fail unless `new` and `old` are the same object down to the bytes."""
    assert type(new) is type(old)
    if isinstance(old, pd.Series):
        pd.testing.assert_series_equal(
            new, old, check_dtype=True, check_index_type=True, check_series_type=True,
            check_names=True, check_exact=True, check_flags=True, check_freq=True,
        )
    assert _canonical_bytes(new) == _canonical_bytes(old)
    assert pickle.dumps(new, protocol=5) == pickle.dumps(old, protocol=5)


# --- real NCBI classifications -----------------------------------------------

def _fixture():
    """The 2,000 real classifications, as strings, in file order."""
    return pd.read_csv(FIXTURE, sep="\t", comment="#", dtype=str, keep_default_na=False)


def golden_bytes(taxids, lineages):
    """Write (taxid, lineage) pairs the way lineage_golden.tsv was written."""
    lines = ["taxid\tlineage\n"] + [f"{t}\t{x}\n" for t, x in zip(taxids, lineages)]
    return "".join(lines).encode("utf-8")


def test_fixture_is_the_one_the_golden_file_was_made_from():
    """The real-data fixture has 2,000 rows, no missing values and no tabs."""
    table = _fixture()
    assert list(table.columns) == ["taxid", "classification"]
    assert len(table) == 2000
    assert table.taxid.is_unique
    assert (table.classification != "").all()


def test_golden_file_series():
    """lineage() of the real classifications, as a Series, matches the golden file byte for byte."""
    table = _fixture()
    result = current.lineage(table.classification)
    assert golden_bytes(table.taxid, result) == GOLDEN.read_bytes()


def test_golden_file_list():
    """lineage() of the real classifications, as a list, matches the golden file byte for byte."""
    table = _fixture()
    result = current.lineage(table.classification.tolist())
    assert golden_bytes(table.taxid, result) == GOLDEN.read_bytes()


def test_real_data_against_legacy_series():
    """The real classifications give the same Series bytes as the old code."""
    series = _fixture().classification
    assert_same_bytes(current.lineage(series), legacy.lineage(series))


def test_real_data_against_legacy_dataframe():
    """The real fixture, as a DataFrame indexed by taxid, gives the same bytes as the old code."""
    table = _fixture().astype({"taxid": "int64"}).set_index("taxid").sort_index()
    assert_same_bytes(current.lineage(table), legacy.lineage(table))
    renamed = table.rename(columns={"classification": "cls"})
    assert_same_bytes(current.lineage(renamed, column="cls"), legacy.lineage(renamed, column="cls"))


# --- seeded synthetic tables -------------------------------------------------

def _name(rng):
    """One clade name: a preferred taxon in some letter case, or a made-up one."""
    if rng.random() < 0.5:
        name = rng.choice(PREFERRED)
        style = rng.randrange(3)
        return name if style == 0 else name.capitalize() if style == 1 else name.upper()
    return f"Taxon{rng.randrange(10**6)}"


def _classification(rng, sep):
    """A classification of 0 to 12 names; 0 gives the empty string."""
    return sep.join(_name(rng) for _ in range(rng.randrange(13)))


def _values(seed, sep="; "):
    """Rows for table `seed`: 1 to 200 classifications drawn from a small pool."""
    rng = random.Random(seed)
    n = rng.choice([1, 2, 3, 7, 20, 200])
    pool = [_classification(rng, sep) for _ in range(max(1, n // 4))]
    return [rng.choice(pool) for _ in range(n)], rng


def _names_in(values, sep):
    """Every distinct clade name in `values`, in first-seen order."""
    seen = {}
    for value in values:
        for name in value.split(sep):
            seen.setdefault(name, None)
    return list(seen)


def opt_series_rangeindex(values, rng):
    """A str Series with the default RangeIndex."""
    return (pd.Series(values),), {}


def opt_series_object(values, rng):
    """An object-dtype Series."""
    return (pd.Series(values, dtype=object),), {}


def opt_series_sorted_gappy_index(values, rng):
    """A named Series on a sorted, non-contiguous int index."""
    index = sorted(rng.sample(range(10**6), len(values)))
    return (pd.Series(values, index=index, name="classification"),), {}


def opt_series_string_index(values, rng):
    """A Series on a sorted string index, with a name of its own."""
    index = sorted({f"r{rng.randrange(10**9):09d}" for _ in range(len(values) * 3)})[:len(values)]
    return (pd.Series(values, index=pd.Index(index, name="pid"), name="tax"),), {}


def opt_dataframe_column(values, rng):
    """A DataFrame with a named sorted index, read through `column`."""
    frame = pd.DataFrame({"cls": values, "other": range(len(values))},
                         index=pd.Index(range(10, 10 + len(values)), name="row"))
    return (frame,), {"column": "cls"}


def opt_list(values, rng):
    """A plain list."""
    return (list(values),), {}


def opt_tuple(values, rng):
    """A tuple, which takes the same path as a list."""
    return (tuple(values),), {}


def opt_custom_list(values, rng):
    """A Series with a custom lowercase preferred_taxa list."""
    names = [x.lower() for x in _names_in(values, "; ")]
    chosen = rng.sample(names, k=len(names) // 2) + rng.sample(PREFERRED, 50)
    return (pd.Series(values),), {"preferred_taxa": chosen}


def opt_custom_set(values, rng):
    """A list with a custom lowercase preferred_taxa set, never empty.

    An empty set is outside the old domain: the old code read the default
    list for it, which test_lineage.py pins as a deliberate change.
    """
    names = {x.lower() for x in _names_in(values, "; ")}
    chosen = set(rng.sample(sorted(names), k=max(1, len(names) // 2)))
    return (list(values),), {"preferred_taxa": chosen}


def opt_case_none(values, rng):
    """case=None, with a preferred set in the exact case of the input."""
    names = _names_in(values, "; ")
    chosen = set(rng.sample(names, k=(len(names) + 1) // 2))
    return (pd.Series(values),), {"preferred_taxa": chosen, "case": None}


def opt_case_upper(values, rng):
    """case=str.upper, with an upper-case preferred set."""
    chosen = {x.upper() for x in PREFERRED}
    return (pd.Series(values),), {"preferred_taxa": chosen, "case": str.upper}


def opt_separators(values, rng):
    """Custom insep and outsep."""
    return (pd.Series([v.replace("; ", "|") for v in values]),), {"insep": "|", "outsep": " / "}


OPTIONS = [
    opt_series_rangeindex, opt_series_object, opt_series_sorted_gappy_index,
    opt_series_string_index, opt_dataframe_column, opt_list, opt_tuple,
    opt_custom_list, opt_custom_set, opt_case_none, opt_case_upper, opt_separators,
]


@pytest.mark.parametrize("option", OPTIONS, ids=lambda f: f.__name__[4:])
@pytest.mark.parametrize("seed", range(40))
def test_seeded_tables_against_legacy(seed, option):
    """Each of 40 seeded tables, under each option set, gives the old bytes."""
    values, rng = _values(seed)
    args, kwargs = option(values, rng)
    assert_same_bytes(current.lineage(*args, **kwargs), legacy.lineage(*args, **kwargs))


@pytest.mark.parametrize("data", [
    pd.Series([], dtype=str),
    pd.Series([], dtype=object),
    pd.Series([], dtype=str, name="classification", index=pd.Index([], dtype="int64", name="pid")),
    pd.DataFrame({"classification": pd.Series([], dtype=str)}),
], ids=["str", "object", "named", "dataframe"])
def test_empty_series_against_legacy(data):
    """An empty Series or DataFrame gives the old empty Series, dtype included."""
    assert_same_bytes(current.lineage(data), legacy.lineage(data))


def test_empty_list_against_legacy():
    """An empty list gives an empty list."""
    assert_same_bytes(current.lineage([]), legacy.lineage([]))


def test_empty_classification_against_legacy():
    """An empty classification string gives an empty lineage, as before."""
    data = pd.Series(["", "Bacteria", ""])
    assert_same_bytes(current.lineage(data), legacy.lineage(data))


# --- default preferred taxa --------------------------------------------------

def test_default_preferred_taxa_unchanged():
    """loadPreferredTaxa() still returns the old set from the shipped file."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        new = current.loadPreferredTaxa()
    old = legacy.loadPreferredTaxa()
    assert type(new) is type(old)
    assert sorted(new) == sorted(old)
