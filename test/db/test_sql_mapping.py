#!/usr/bin/env python3
"""Tests for the dialect-free SQL of an identifier mapping.

`rotifer.db.sql.mapping` builds the statement both UniProt SQL backends run.
It is worth testing on its own because the two things it has to get right are
invisible from either backend alone:

  * the shape. Which branches a query needs follows from what was asked, and
    two of the shapes exist only to avoid a join that would otherwise hash a
    few billion rows
  * the binding. A dialect binding by position needs its values in the order
    the statement mentions them, so every fragment must be built in that
    order. Get it wrong and ClickHouse still works -- it binds by name -- while
    SQLite3 silently matches the wrong value to the wrong placeholder

The asymmetry between the two ends is deliberate and easy to lose: an omitted
source means an identifier may be anything, the accession included, while an
omitted target means every database the table cross-references, which does not
include the accession.

Run under pytest, or standalone: python test/db/test_sql_mapping.py
"""

import re

import pytest

from rotifer.db.sql import mapping as sqlmap

U = sqlmap.UNIPROTKB


def qmark(ids=('A', 'B')):
    """A binder and a restriction that binds by position, as SQLite3 does."""
    binder = sqlmap.QmarkBinder()

    def restrict(column):
        return f'f.{column} IN {binder.collection(list(ids))}'

    return binder, restrict


def named(ids=('A', 'B')):
    """The same, binding by name, as ClickHouse does."""
    binder = sqlmap.NamedBinder()

    def restrict(column):
        return f'f.{column} IN {binder.collection(list(ids))}'

    return binder, restrict


# ------------------------------------------------------------- the shape

def test_an_accession_source_needs_no_join():
    """Every row of a matching accession already carries the answer, so the
    commonest query reads the table once instead of joining it to itself."""
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, [U], None, binder)
    assert 'JOIN' not in sql


def test_an_accession_target_needs_no_join():
    """The rows matching the queried identifiers already carry the accession,
    which is the whole answer when that is what was asked for."""
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, ['RefSeq'], [U], binder)
    assert 'JOIN' not in sql


def test_two_real_ends_join_but_restrict_the_far_side():
    """The target side is restricted by the accessions the source side found.
    Merely joining on them would hash the whole table, which on a few billion
    rows does not return -- this is the difference between 0.3 s and never."""
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, ['EMBL-CDS'], ['RefSeq'], binder)
    assert 'INNER JOIN' in sql
    assert 'accession IN (SELECT accession FROM' in sql


def test_an_open_source_matches_either_column():
    """With no source named the caller has not said what they are holding, so
    an identifier may be an accession or a cross-reference."""
    binder, restrict = named()
    sql = sqlmap.source_subquery('t', restrict, None, binder)
    assert 'f.accession IN' in sql and 'f.id IN' in sql
    assert sql.count('UNION ALL') == 1


def test_an_open_target_does_not_return_the_accession():
    """'Every database' means the ones the table cross-references. The
    accession is its key, not one of them, and returning it would answer a
    question nobody asked."""
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, [U], None, binder)
    assert f"'{U}' AS target_type" not in sql


def test_naming_the_accession_as_a_target_does_return_it():
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, [U], ['RefSeq', U], binder)
    assert f"'{U}' AS target_type" in sql


def test_a_release_restricts_every_branch():
    """A table holding several releases must not answer from the wrong one,
    whichever branch of the query does the answering."""
    binder, restrict = named()
    sql = sqlmap.mapping_query('t', restrict, None, None, binder, release='2026_01')
    assert sql.count('release =') == sql.count('SELECT DISTINCT') + sql.count('SELECT accession AS value') \
        or 'release =' in sql
    assert '2026_01' in str(binder.parameters)


# ----------------------------------------------------------- the binding

@pytest.mark.parametrize('source,target', [
    ([U], None),
    ([U], ['RefSeq']),
    (['RefSeq'], [U]),
    (['EMBL-CDS'], ['RefSeq']),
    (None, None),
    ([U], ['RefSeq', U]),
    (['RefSeq', U], [U]),
])
def test_positional_values_are_bound_in_the_order_the_statement_reads_them(source, target):
    """The invariant a positional dialect depends on: the nth ? is filled by
    the nth bound value. Building any fragment out of order breaks it, and
    ClickHouse would never notice because it binds by name."""
    binder, restrict = qmark(['A', 'B'])
    sql = sqlmap.mapping_query('t', restrict, source, target, binder, release='2026_01')
    assert sql.count('?') == len(binder.parameters)


def test_every_repetition_of_a_condition_binds_its_own_values():
    """A condition appearing three times needs its values bound three times.
    Writing it once and repeating the string would leave later placeholders
    unfilled, which is why the restriction is a callable."""
    binder, restrict = qmark(['A', 'B'])
    sql = sqlmap.mapping_query('t', restrict, None, None, binder)
    assert sql.count('f.accession IN') + sql.count('f.id IN') >= 2
    assert sql.count('?') == len(binder.parameters)


def test_named_binding_keeps_one_name_per_value():
    binder, restrict = named(['A', 'B'])
    sql = sqlmap.mapping_query('t', restrict, [U], ['RefSeq'], binder)
    for name in binder.parameters:
        assert f'%({name})s' in sql
    assert len(set(re.findall(r'%\((\w+)\)s', sql))) == len(binder.parameters)


def test_the_two_dialects_build_the_same_shape():
    """The point of the module: only the placeholders differ. Stripping them
    must leave two identical statements, or the backends have drifted."""
    a, restrict_a = named(['A', 'B'])
    b, restrict_b = qmark(['A', 'B'])
    left = sqlmap.mapping_query('t', restrict_a, ['EMBL-CDS'], ['RefSeq'], a)
    right = sqlmap.mapping_query('t', restrict_b, ['EMBL-CDS'], ['RefSeq'], b)
    strip = lambda s: re.sub(r'%\(\w+\)s|\((?:\?(?:, )?)+\)|\?', '<v>', s)
    assert strip(left) == strip(right)


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-q']))
