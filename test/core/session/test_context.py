#!/usr/bin/env python3
"""Tests for rotifer.core.session._context and describe().

Run under pytest, or standalone: python test/core/session/test_context.py
"""

import sys

import numpy as np
import pandas as pd
import pytest

from rotifer.core.session import default_namespace, describe
from rotifer.core.session._context import (
    compare_environment, environment, notebook_context)


def test_default_namespace_outside_ipython():
    assert default_namespace() is sys.modules["__main__"].__dict__


def test_environment_keys():
    env = environment()
    assert set(env) == {"python", "rotifer", "packages", "tools"}
    assert set(env["packages"]) == {
        "pandas", "numpy", "biopython", "pyarrow", "cloudpickle"}
    assert env["python"].startswith(f"{sys.version_info.major}.")
    assert env["packages"]["pandas"] == pd.__version__
    assert env["tools"] == {}


def test_notebook_context_outside_ipython():
    assert notebook_context() == {}


def test_compare_environment():
    cur = environment()
    assert compare_environment(cur, cur) == []
    saved = {**cur, "python": "2.7.0",
             "packages": {**cur["packages"], "pandas": "0.1", "pyarrow": None,
                          "numpy": None}}
    saved["packages"]["cloudpickle"] = "9.9"
    cur = {**cur, "packages": {**cur["packages"], "numpy": "1.0",
                               "cloudpickle": None}}
    text = "\n".join(compare_environment(saved, cur))
    assert "Python version differs: saved with 2.7.0" in text
    assert "pandas version differs: saved with 0.1" in text
    assert "numpy was not installed when saved" in text
    assert "cloudpickle 9.9 was used when saved, but it is not installed" in text


def test_describe_dataframe():
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    d = describe(df)
    assert d["kind"] == "dataframe" and d["shape"] == [2, 2]
    assert d["length"] == 2 and d["columns"] == ["a", "b"]
    assert d["dtypes"]["a"] == "int64" and len(d["preview"]) <= 200


def test_describe_wide_dataframe_truncates_columns():
    df = pd.DataFrame(np.zeros((1, 80)))
    assert len(describe(df)["columns"]) == 50


def test_describe_ndarray():
    d = describe(np.zeros((3, 4)))
    assert d["kind"] == "ndarray" and d["shape"] == [3, 4] and d["length"] == 3


def test_describe_builtins():
    assert describe([1, 2, 3]) == {"kind": "sequence", "length": 3,
                                   "preview": "[1, 2, 3]"}
    assert describe({"a": 1})["kind"] == "mapping"
    assert describe(5) == {"kind": "scalar", "preview": "5"}
    assert len(describe("x" * 1000)["preview"]) <= 200


def test_describe_unknown_object_never_calls_len_or_repr():
    calls = []

    class Thing:
        def __len__(self):
            calls.append("len")
            raise RuntimeError("no len")

        def __repr__(self):
            calls.append("repr")
            raise RuntimeError("no repr")

        __str__ = __repr__

    d = describe(Thing())
    assert calls == []
    assert d["kind"] == "object"
    assert d["type"].endswith("Thing")
    assert d["preview"] == "<%s object>" % d["type"]
    assert "length" not in d


def test_describe_list_subclass_is_not_inspected():
    calls = []

    class Fancy(list):
        def __len__(self):
            calls.append("len")
            return 1

        def __repr__(self):
            calls.append("repr")
            return "x"

    assert describe(Fancy())["kind"] == "object"
    assert calls == []


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
