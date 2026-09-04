#!/usr/bin/env python3
"""Offline contracts for NCBI gene-neighborhood retrieval.

These tests characterize the behavior required by workflows that resolve
protein identifiers through IPGs and retrieve genomic neighborhoods.
They do not access NCBI, local mirrors, institutional databases, or the
network.
"""

import pandas as pd

from rotifer.db import ncbi
from rotifer.db.core import BaseCursor


TARGET = "TSPA_QUERY"


class OfflineGeneNeighborhoodCursor(ncbi.GeneNeighborhoodCursor):
    """GeneNeighborhoodCursor without construction of external backends."""

    def __init__(
        self,
        readers=(),
        cursors=None,
    ):
        BaseCursor.__init__(
            self,
            progress=False,
        )

        self.readers = list(readers)
        self.writers = []
        self.cursors = dict(cursors or {})

        self.column = "pid"
        self.before = 10
        self.after = 10

        self.save = None
        self.tries = 1


class EmptyMissingBackend(BaseCursor):
    """Backend that records a missing accession but yields no result."""

    def __init__(self, accession):
        super().__init__(
            progress=False,
        )
        self.calls = []
        self._missing = {
            accession: [
                "backend did not resolve accession",
                "EmptyMissingBackend",
                False,
            ],
        }

    def fetchone(
        self,
        accessions,
        ipgs=None,
    ):
        self.calls.append(
            set(accessions)
        )
        return iter(())


class ResolvingBackend(BaseCursor):
    """Backend that resolves requested protein accessions in memory."""

    def __init__(self):
        super().__init__(
            progress=False,
        )
        self.calls = []

    def fetchone(
        self,
        accessions,
        ipgs=None,
    ):
        targets = set(accessions)
        self.calls.append(
            targets.copy()
        )

        for accession in targets:
            yield pd.DataFrame(
                [
                    {
                        "pid": accession,
                    },
                ]
            )


def target_ipgs():
    """IPG fixture containing a query, synonym, and representative."""

    return pd.DataFrame(
        [
            {
                "id": 1,
                "pid": TARGET,
                "representative": "TSPA_REP",
            },
            {
                "id": 1,
                "pid": "TSPA_ALT",
                "representative": "TSPA_REP",
            },
            {
                "id": 2,
                "pid": "UNRELATED",
                "representative": "UNRELATED_REP",
            },
        ]
    )


def empty_ipgs():
    """Structurally valid empty IPG table."""

    return pd.DataFrame(
        columns=[
            "id",
            "pid",
            "representative",
        ]
    )


def test_getids_expands_protein_identifiers_from_same_ipg():
    cursor = OfflineGeneNeighborhoodCursor()

    identifiers = cursor.getids(
        TARGET,
        ipgs=target_ipgs(),
    )

    assert identifiers == {
        TARGET,
        "TSPA_ALT",
        "TSPA_REP",
    }


def test_fetchone_records_target_absent_from_supplied_ipgs():
    cursor = OfflineGeneNeighborhoodCursor()

    results = list(
        cursor.fetchone(
            [TARGET],
            ipgs=empty_ipgs(),
        )
    )

    assert results == []
    assert cursor.missing_ids() == {
        TARGET,
    }
    assert cursor._missing[TARGET][0] == "Not found in IPGs"
    assert cursor._missing[TARGET][2] is False


def test_fetchall_returns_valid_empty_dataframe_when_target_is_unresolved():
    cursor = OfflineGeneNeighborhoodCursor()

    result = cursor.fetchall(
        [TARGET],
        ipgs=empty_ipgs(),
    )

    assert isinstance(
        result,
        pd.DataFrame,
    )
    assert result.empty


def test_fetchone_preserves_missing_state_from_backend_without_results():
    backend = EmptyMissingBackend(
        TARGET,
    )

    cursor = OfflineGeneNeighborhoodCursor(
        readers=["memory"],
        cursors={
            "memory": backend,
        },
    )

    results = list(
        cursor.fetchone(
            [TARGET],
            ipgs=target_ipgs(),
        )
    )

    assert results == []
    assert cursor._missing == backend._missing


def test_fetchone_continues_after_missing_reader_and_clears_resolved_target():
    missing_backend = EmptyMissingBackend(
        TARGET,
    )
    resolving_backend = ResolvingBackend()

    cursor = OfflineGeneNeighborhoodCursor(
        readers=[
            "missing",
            "resolving",
        ],
        cursors={
            "missing": missing_backend,
            "resolving": resolving_backend,
        },
    )

    results = list(
        cursor.fetchone(
            [TARGET],
            ipgs=target_ipgs(),
        )
    )

    assert len(results) == 1

    assert missing_backend.calls == [
        {
            TARGET,
        },
    ]

    assert resolving_backend.calls == [
        {
            TARGET,
        },
    ]

    assert cursor.missing_ids() == set()
    assert missing_backend.missing_ids() == set()
