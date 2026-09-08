#!/usr/bin/env python3
"""Tests for answering a mapping no single backend can answer.

The backends do not hold the same vocabulary. `idmapping.dat` maps UniRef and
knows nothing of AlphaFoldDB; the web service maps AlphaFoldDB and knows
nothing of UniRef. A query between the two therefore finds every backend able
to serve one end and none able to serve both, and before this it came back
empty with nothing in `missing` -- indistinguishable from a mapping that does
not exist.

It is answerable all the same, because every row here is joined through a
UniProtKB accession: resolve the identifiers to accessions with the backend
that knows them, then ask the other backend about those. What these tests pin
down is that the two halves are joined on the accession and reported under the
identifier the caller actually gave, that a backend able to answer in one step
is still used in one step, and that borrowing a backend for the pivot leaves
no trace of accessions nobody asked about in its registry.

Nothing here touches a network or a server.

Run under pytest, or standalone: python test/db/uniprot/test_mapping_bridge.py
"""

import pandas as pd
import pytest

import rotifer.db.core
import rotifer.db.methods
from rotifer.db import uniprot

COLUMNS = ['source', 'source_type', 'accession', 'target', 'target_type']


class Backend(rotifer.db.methods.MappingCursor, rotifer.db.core.BaseCursor):
    """A backend holding a fixed mapping table over a fixed vocabulary."""

    def __init__(self, rows, vocabulary, content=None):
        super().__init__(progress=False)
        self.table = pd.DataFrame(list(rows), columns=COLUMNS)
        self.vocabulary = set(vocabulary)
        self.asked = []
        self._content = content

    def databases(self):
        return self.vocabulary

    def content_id(self):
        return self._content

    def fetchone(self, accessions, source=None, target=None, *args, **kwargs):
        wanted = self.parse_ids(accessions)
        self.asked.append((sorted(wanted), source, target))
        rows = self.table[self.table.source.isin(wanted)]
        if source:
            rows = rows[rows.source_type.isin(source)]
        if target:
            rows = rows[rows.target_type.isin(target)]
        if not rows.empty:
            yield rows.reset_index(drop=True)
        self.update_missing(wanted - set(rows.source), 'Not found.', retry=False)


def delegator(**backends):
    """A MappingCursor over the given backends, asked in the order given."""

    class Cursor(uniprot.MappingCursor):
        def reset_cursors(self):
            self.cursors = dict(getattr(self, '_backends', {}))

    cursor = Cursor.__new__(Cursor)
    cursor._backends = dict(backends)
    uniprot.MappingCursor.__init__(cursor, readers=list(backends), progress=False)
    return cursor


#: What the table has: UniRef in, accession out. No AlphaFoldDB.
def table():
    return Backend(
        [('UniRef50_C7N6L9', 'UniRef50', 'C7N6L9', 'C7N6L9', 'UniProtKB-AC'),
         ('UniRef50_C7N6L9', 'UniRef50', 'C7N6L9', 'X_1', 'RefSeq')],
        {'UniRef50', 'RefSeq'})


#: What the service has: accession in, AlphaFoldDB out. No UniRef.
def service():
    return Backend(
        [('C7N6L9', 'UniProtKB-AC', 'C7N6L9', 'C7N6L9', 'AlphaFoldDB'),
         ('C7N6L9', 'UniProtKB-AC', 'C7N6L9', 'X_1', 'RefSeq')],
        {'AlphaFoldDB', 'RefSeq'})


# ------------------------------------------------ the query that failed

def test_a_query_neither_backend_spans_is_answered():
    cursor = delegator(table=table(), service=service())
    found = cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert len(found) == 1


def test_the_answer_is_reported_under_the_identifier_that_was_given():
    """The caller asked about a UniRef cluster, not about the accession it
    had to be turned into on the way."""
    cursor = delegator(table=table(), service=service())
    row = cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50',
                          target='AlphaFoldDB').iloc[0]
    assert row.source == 'UniRef50_C7N6L9'
    assert row.source_type == 'UniRef50'
    assert row.target_type == 'AlphaFoldDB'


def test_the_accession_that_joined_the_halves_is_shown():
    cursor = delegator(table=table(), service=service())
    row = cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50',
                          target='AlphaFoldDB').iloc[0]
    assert row.accession == 'C7N6L9'


def test_the_second_backend_is_asked_about_accessions():
    cursor = delegator(table=table(), service=service())
    cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert cursor.cursors['service'].asked[-1][0] == ['C7N6L9']


# ------------------------------------------------- when not to bridge

def test_a_backend_that_spans_the_query_answers_it_alone():
    """One step is still one step. A bridge is what happens when no backend
    can do that, not a replacement for it."""
    cursor = delegator(service=service())
    cursor.fetchall(['C7N6L9'], source='UniProtKB-AC', target='RefSeq')
    assert len(cursor.cursors['service'].asked) == 1


def test_an_open_source_is_not_bridged():
    """With the source left open every backend is asked with the identifiers
    as they stand, so a second pass would reach nothing the first did not."""
    cursor = delegator(table=table(), service=service())
    cursor.fetchall(['UniRef50_C7N6L9'], target='AlphaFoldDB')
    assert all(a[0] != ['C7N6L9'] for a in cursor.cursors['service'].asked)


def test_nothing_is_bridged_to_a_database_nobody_has():
    """The pivot costs a query of its own, so it is not paid for an answer
    that cannot exist."""
    cursor = delegator(table=table(), service=service())
    with pytest.raises(ValueError):
        cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='NoSuchDB')


def test_a_backend_repeating_another_is_not_bridged_through_twice():
    first = Backend(table().table.values.tolist(), {'UniRef50','RefSeq'}, content='2026_01')
    second = Backend(table().table.values.tolist(), {'UniRef50','RefSeq'}, content='2026_01')
    cursor = delegator(first=first, second=second, service=service())
    cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert second.asked == []


# ------------------------------------------------------- what it leaves

def test_borrowing_a_backend_leaves_its_registry_as_it_was():
    """The pivot asks about accessions the caller never mentioned, and a
    cluster's other members are exactly the accessions a backend is likely
    not to have. Recording one would report an identifier nobody asked
    about, which is worse than saying nothing."""
    many = Backend(
        [('UniRef50_C7N6L9', 'UniRef50', 'C7N6L9', 'C7N6L9', 'UniProtKB-AC'),
         ('UniRef50_C7N6L9', 'UniRef50', 'UNKNOWN', 'UNKNOWN', 'UniProtKB-AC')],
        {'UniRef50'})
    # It knows C7N6L9 and has never heard of the other member
    partial = Backend(
        [('C7N6L9', 'UniProtKB-AC', 'C7N6L9', 'C7N6L9', 'AlphaFoldDB')],
        {'AlphaFoldDB'})
    cursor = delegator(many=many, partial=partial)
    cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert 'UNKNOWN' in partial.asked[-1][0]      # it really was asked
    assert partial._missing == {}
    assert 'UNKNOWN' not in cursor.missing_ids()


def test_an_answered_identifier_is_not_left_looking_missing():
    cursor = delegator(table=table(), service=service())
    cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert cursor.missing_ids() == set()


def test_an_identifier_nobody_could_look_for_is_reported():
    """Every backend passed over for want of the vocabulary used to leave an
    empty answer with no reason attached, which reads exactly like a mapping
    that does not exist."""
    cursor = delegator(table=table(), service=service())
    found = cursor.fetchall(['UniRef50_NOPE'], source='UniRef50', target='AlphaFoldDB')
    assert found.empty
    assert cursor.missing_ids() == {'UniRef50_NOPE'}


# --------------------------------------------------------- the join

def test_every_accession_of_an_identifier_is_followed():
    """A cluster stands for more than one entry, so the two halves join on
    the accession rather than pairing off in order."""
    many = Backend(
        [('UniRef50_C7N6L9', 'UniRef50', 'C7N6L9', 'C7N6L9', 'UniProtKB-AC'),
         ('UniRef50_C7N6L9', 'UniRef50', 'P00750', 'P00750', 'UniProtKB-AC')],
        {'UniRef50'})
    both = Backend(
        [('C7N6L9', 'UniProtKB-AC', 'C7N6L9', 'C7N6L9', 'AlphaFoldDB'),
         ('P00750', 'UniProtKB-AC', 'P00750', 'P00750', 'AlphaFoldDB')],
        {'AlphaFoldDB'})
    cursor = delegator(many=many, both=both)
    found = cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50', target='AlphaFoldDB')
    assert sorted(found.target) == ['C7N6L9', 'P00750']
    assert set(found.source) == {'UniRef50_C7N6L9'}


def test_a_database_one_backend_already_served_is_not_bridged_for():
    """RefSeq is in both vocabularies, so the first pass answers it. Only
    AlphaFoldDB is still owed, and that is all the second pass asks for."""
    cursor = delegator(table=table(), service=service())
    cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50',
                    target=['RefSeq', 'AlphaFoldDB'])
    assert cursor.cursors['service'].asked[-1][2] == ['AlphaFoldDB']


def test_both_halves_of_a_split_target_come_back():
    cursor = delegator(table=table(), service=service())
    found = cursor.fetchall(['UniRef50_C7N6L9'], source='UniRef50',
                            target=['RefSeq', 'AlphaFoldDB'])
    assert sorted(set(found.target_type)) == ['AlphaFoldDB', 'RefSeq']


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
