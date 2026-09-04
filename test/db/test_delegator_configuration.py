#!/usr/bin/env python3
"""Contract tests for DelegatorCursor backend configuration.

These tests exercise backend alias resolution and module loading without
accessing external databases, the network, or institutional infrastructure.
"""

import importlib
import sys
from types import SimpleNamespace

import pytest

from rotifer.db.delegator import DelegatorCursor


class ConfigurationDelegator(DelegatorCursor):
    """Delegator whose backend construction is disabled for config tests."""

    def reset_cursors(self):
        self.cursors = {}


def test_cursor_modules_loads_configured_readers_and_writers(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.setattr(
        test_module,
        "config",
        {
            "readers": {
                "reader": "fake.reader",
            },
            "writers": {
                "writer": "fake.writer",
            },
        },
        raising=False,
    )

    modules = {
        "fake.reader": SimpleNamespace(
            __name__="fake.reader",
        ),
        "fake.writer": SimpleNamespace(
            __name__="fake.writer",
        ),
    }

    imported = []

    def fake_import_module(name):
        imported.append(name)
        return modules[name]

    monkeypatch.setattr(
        importlib,
        "import_module",
        fake_import_module,
    )

    cursor = ConfigurationDelegator(
        readers=["reader"],
        writers=["writer"],
        progress=False,
    )

    resolved = cursor._cursor_modules

    assert resolved == {
        "reader": modules["fake.reader"],
        "writer": modules["fake.writer"],
    }

    assert set(imported) == {
        "fake.reader",
        "fake.writer",
    }


def test_cursor_modules_requires_module_config(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.delattr(
        test_module,
        "config",
        raising=False,
    )

    cursor = ConfigurationDelegator(
        progress=False,
    )

    with pytest.raises(
        ValueError,
        match='No attribute "config"',
    ):
        cursor._cursor_modules


def test_cursor_modules_requires_readers_mapping(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.setattr(
        test_module,
        "config",
        {
            "writers": {},
        },
        raising=False,
    )

    cursor = ConfigurationDelegator(
        progress=False,
    )

    with pytest.raises(
        ValueError,
        match="Missing dictionary of reader modules",
    ):
        cursor._cursor_modules


def test_cursor_modules_requires_writers_mapping(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.setattr(
        test_module,
        "config",
        {
            "readers": {},
        },
        raising=False,
    )

    cursor = ConfigurationDelegator(
        progress=False,
    )

    with pytest.raises(
        ValueError,
        match="Missing dictionary of writer modules",
    ):
        cursor._cursor_modules


def test_cursor_modules_rejects_unknown_backend_alias(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.setattr(
        test_module,
        "config",
        {
            "readers": {},
            "writers": {},
        },
        raising=False,
    )

    cursor = ConfigurationDelegator(
        readers=["unknown"],
        progress=False,
    )

    with pytest.raises(
        ValueError,
        match='Missing module name "unknown"',
    ):
        cursor._cursor_modules


def test_cursor_modules_reports_backend_import_failure(
    monkeypatch,
):
    test_module = sys.modules[__name__]

    monkeypatch.setattr(
        test_module,
        "config",
        {
            "readers": {
                "reader": "fake.backend",
            },
            "writers": {},
        },
        raising=False,
    )

    def fail_import(name):
        raise ModuleNotFoundError(name)

    monkeypatch.setattr(
        importlib,
        "import_module",
        fail_import,
    )

    cursor = ConfigurationDelegator(
        readers=["reader"],
        progress=False,
    )

    with pytest.raises(
        ImportError,
        match="Unable to load module fake.backend",
    ):
        cursor._cursor_modules
