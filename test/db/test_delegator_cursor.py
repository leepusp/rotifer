#!/usr/bin/env python3
"""Contract tests for Rotifer sequential cursor delegation.

The test backends are entirely in memory. They exercise delegation,
fallback, missing-accession bookkeeping, and writers without external
databases, network access, or institutional configuration.
"""

from rotifer.db.core import BaseCursor
from rotifer.db.delegator import SequentialDelegatorCursor


class MemoryCursor(BaseCursor):
    """Minimal cursor backend used only for delegation contract tests."""

    def __init__(self, available=()):
        super().__init__(progress=False)
        self.available = set(available)
        self.calls = []
        self.inserted = []

    def getids(self, obj):
        if obj is None:
            return set()
        if isinstance(obj, str):
            return {obj}
        return set(obj)

    def __getitem__(self, accessions, *args, **kwargs):
        targets = self.parse_ids(accessions)
        self.calls.append(("getitem", targets.copy()))

        found = targets.intersection(self.available)
        missing = targets - found

        self.remove_missing(found)

        if missing:
            self.update_missing(
                missing,
                error="not found",
                retry=False,
            )

        return sorted(found)

    def fetchone(self, accessions, *args, **kwargs):
        targets = self.parse_ids(accessions)
        self.calls.append(("fetchone", targets.copy()))

        found = targets.intersection(self.available)
        missing = targets - found

        self.remove_missing(found)

        if missing:
            self.update_missing(
                missing,
                error="not found",
                retry=False,
            )

        for accession in sorted(found):
            yield accession

    def insert(self, obj):
        self.inserted.append(obj)


class MemoryDelegator(SequentialDelegatorCursor):
    """Sequential delegator with explicitly supplied in-memory backends."""

    def __init__(
        self,
        backends,
        readers,
        writers=None,
    ):
        self._test_backends = backends

        super().__init__(
            readers=readers,
            writers=[] if writers is None else writers,
            progress=False,
        )

    def reset_cursors(self):
        names = set(self.readers).union(self.writers)

        self.cursors = {
            name: self._test_backends[name]
            for name in names
        }

    def getids(self, obj, *args, **kwargs):
        if obj is None:
            return set()
        if isinstance(obj, str):
            return {obj}
        return set(obj)


def test_getitem_uses_reader_fallback_only_for_unresolved_ids():
    first = MemoryCursor({"A"})
    second = MemoryCursor({"B"})

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = cursor[["A", "B"]]

    assert result == ["A", "B"]

    assert first.calls == [
        ("getitem", {"A", "B"}),
    ]

    assert second.calls == [
        ("getitem", {"B"}),
    ]

    assert cursor.missing_ids() == set()


def test_getitem_records_accession_missing_from_all_readers():
    first = MemoryCursor()
    second = MemoryCursor()

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = cursor["MISSING"]

    assert result == []
    assert cursor.missing_ids() == {"MISSING"}


def test_fetchone_uses_reader_fallback_only_for_unresolved_ids():
    first = MemoryCursor({"A"})
    second = MemoryCursor({"B"})

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = list(cursor.fetchone(["A", "B"]))

    assert result == ["A", "B"]

    assert first.calls == [
        ("fetchone", {"A", "B"}),
    ]

    assert second.calls == [
        ("fetchone", {"B"}),
    ]

    assert cursor.missing_ids() == set()


def test_fetchone_records_accession_missing_from_all_readers():
    first = MemoryCursor()
    second = MemoryCursor()

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = list(cursor.fetchone(["MISSING"]))

    assert result == []
    assert cursor.missing_ids() == {"MISSING"}


def test_fetchone_sends_results_to_configured_writer():
    reader = MemoryCursor({"A"})
    writer = MemoryCursor()

    cursor = MemoryDelegator(
        backends={
            "reader": reader,
            "writer": writer,
        },
        readers=["reader"],
        writers=["writer"],
    )

    result = list(cursor.fetchone(["A"]))

    assert result == ["A"]
    assert writer.inserted == ["A"]


def test_later_reader_clears_missing_state_from_previous_reader():
    first = MemoryCursor()
    second = MemoryCursor({"A"})

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = list(cursor.fetchone(["A"]))

    assert result == ["A"]

    assert first.missing_ids() == set()
    assert second.missing_ids() == set()
    assert cursor.missing_ids() == set()


def test_fetchall_consumes_sequential_fallback():
    first = MemoryCursor({"A"})
    second = MemoryCursor({"B"})

    cursor = MemoryDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    result = cursor.fetchall(["A", "B"])

    assert set(result) == {"A", "B"}
    assert cursor.missing_ids() == set()


# ---------------------------------------------------------------------------
# Shared attribute propagation
# ---------------------------------------------------------------------------


class AttributeCursor(BaseCursor):
    """Backend exposing attributes used to test delegator propagation."""

    def __init__(
        self,
        shared="backend-default",
        nullable="backend-nullable-default",
    ):
        super().__init__(progress=False)
        self.shared = shared
        self.nullable = nullable
        self.local_only = "backend-local"


class AttributeDelegator(MemoryDelegator):
    """In-memory delegator exposing the generic shared-attribute policy."""

    _shared_attributes = [
        "shared",
        "nullable",
    ]

    _nullable_attributes = frozenset({
        "nullable",
    })


def make_attribute_delegator():
    backend = AttributeCursor()

    cursor = AttributeDelegator(
        backends={
            "reader": backend,
        },
        readers=["reader"],
    )

    return cursor, backend


def test_shared_attribute_is_propagated_to_backend():
    cursor, backend = make_attribute_delegator()

    cursor.shared = "updated"

    assert cursor.shared == "updated"
    assert backend.shared == "updated"


def test_non_nullable_none_preserves_backend_value():
    cursor, backend = make_attribute_delegator()

    original = backend.shared

    cursor.shared = None

    assert cursor.shared is None
    assert backend.shared == original


def test_nullable_none_is_propagated_to_backend():
    cursor, backend = make_attribute_delegator()

    cursor.nullable = None

    assert cursor.nullable is None
    assert backend.nullable is None


def test_non_shared_attribute_is_not_propagated():
    cursor, backend = make_attribute_delegator()

    cursor.local_only = "delegator-value"

    assert cursor.local_only == "delegator-value"
    assert backend.local_only == "backend-local"


def test_shared_attribute_propagates_to_all_backends():
    first = AttributeCursor()
    second = AttributeCursor()

    cursor = AttributeDelegator(
        backends={
            "first": first,
            "second": second,
        },
        readers=["first", "second"],
    )

    cursor.shared = "common"

    assert first.shared == "common"
    assert second.shared == "common"


# ---------------------------------------------------------------------------
# Shared attributes during backend construction
# ---------------------------------------------------------------------------


class ConstructionBackend(BaseCursor):
    """Backend recording values received when the delegator constructs it."""

    def __init__(
        self,
        shared="backend-default",
        optional="backend-optional-default",
    ):
        super().__init__(progress=False)
        self.shared = shared
        self.optional = optional


class ConstructionDelegator(SequentialDelegatorCursor):
    """Delegator exercising the real reset_cursors construction path."""

    _shared_attributes = [
        "shared",
        "optional",
    ]

    def __init__(
        self,
        shared="configured",
        optional=None,
    ):
        self.shared = shared
        self.optional = optional

        super().__init__(
            readers=["reader"],
            writers=[],
            progress=False,
        )

    @property
    def _cursor_modules(self):
        class BackendModule:
            ConstructionDelegator = ConstructionBackend

        return {
            "reader": BackendModule,
        }

    def getids(self, obj, *args, **kwargs):
        if obj is None:
            return set()
        if isinstance(obj, str):
            return {obj}
        return set(obj)


def test_reset_cursors_passes_non_none_shared_attributes():
    cursor = ConstructionDelegator(
        shared="configured-value",
    )

    backend = cursor.cursors["reader"]

    assert backend.shared == "configured-value"


def test_reset_cursors_preserves_backend_default_for_none():
    cursor = ConstructionDelegator(
        optional=None,
    )

    backend = cursor.cursors["reader"]

    assert cursor.optional is None
    assert backend.optional == "backend-optional-default"
