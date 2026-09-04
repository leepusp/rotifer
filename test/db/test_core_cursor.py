#!/usr/bin/env python3
"""Contract tests for rotifer.db.core.BaseCursor.

These tests define behaviour expected from the base cursor abstraction itself.
They do not access remote services, local biological databases, or institutional
infrastructure.
"""

import pytest

from rotifer.db.core import BaseCursor


# ---------------------------------------------------------------------------
# Identifier normalization
# ---------------------------------------------------------------------------


def test_parse_ids_from_single_string():
    cursor = BaseCursor()

    assert cursor.parse_ids("ABC123") == {"ABC123"}


def test_parse_ids_from_comma_separated_string():
    cursor = BaseCursor()

    assert cursor.parse_ids("ABC123,DEF456") == {
        "ABC123",
        "DEF456",
    }


def test_parse_ids_removes_duplicates():
    cursor = BaseCursor()

    assert cursor.parse_ids(
        ["ABC123", "ABC123", "DEF456"]
    ) == {
        "ABC123",
        "DEF456",
    }


def test_parse_ids_converts_values_to_strings_by_default():
    cursor = BaseCursor()

    assert cursor.parse_ids([1, 2, 2]) == {
        "1",
        "2",
    }


def test_parse_ids_can_preserve_input_types():
    cursor = BaseCursor()

    assert cursor.parse_ids(
        [1, 2, 2],
        as_string=False,
    ) == {
        1,
        2,
    }


# ---------------------------------------------------------------------------
# Missing-accession registry
# ---------------------------------------------------------------------------


def test_missing_registry_starts_empty():
    cursor = BaseCursor()

    assert cursor.missing.empty
    assert cursor.missing_ids() == set()


def test_update_missing_registers_accession():
    cursor = BaseCursor()

    cursor.update_missing(
        ["ABC123"],
        error="temporary error",
        retry=True,
    )

    assert cursor.missing_ids() == {"ABC123"}

    row = cursor.missing.loc["ABC123"]

    assert row["error"] == "temporary error"
    assert bool(row["retry"]) is True


def test_missing_ids_can_filter_by_retry_state():
    cursor = BaseCursor()

    cursor.update_missing(
        ["RETRY"],
        error="temporary error",
        retry=True,
    )

    cursor.update_missing(
        ["STOP"],
        error="permanent error",
        retry=False,
    )

    assert cursor.missing_ids(retry=True) == {
        "RETRY",
    }

    assert cursor.missing_ids(retry=False) == {
        "STOP",
    }


def test_giveup_error_is_not_retryable():
    cursor = BaseCursor()

    cursor.giveup.add("permanent")

    cursor.update_missing(
        ["ABC123"],
        error="permanent database error",
    )

    assert cursor.missing_ids(retry=False) == {
        "ABC123",
    }


def test_remove_missing_removes_selected_accession():
    cursor = BaseCursor()

    cursor.update_missing(
        ["ABC123", "DEF456"],
        error="temporary error",
        retry=True,
    )

    cursor.remove_missing(["ABC123"])

    assert cursor.missing_ids() == {
        "DEF456",
    }


def test_remove_missing_without_arguments_clears_registry():
    cursor = BaseCursor()

    cursor.update_missing(
        ["ABC123", "DEF456"],
        error="temporary error",
        retry=True,
    )

    old = cursor.remove_missing()

    assert cursor.missing_ids() == set()
    assert set(old) == {
        "ABC123",
        "DEF456",
    }


# ---------------------------------------------------------------------------
# Abstract cursor interface
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "method,args",
    [
        ("getids", (object(),)),
        ("__getitem__", ("ABC123",)),
        ("fetchone", (["ABC123"],)),
        ("fetchall", (["ABC123"],)),
    ],
)
def test_base_cursor_abstract_methods_raise_not_implemented(
    method,
    args,
):
    cursor = BaseCursor()

    with pytest.raises(NotImplementedError):
        getattr(cursor, method)(*args)
