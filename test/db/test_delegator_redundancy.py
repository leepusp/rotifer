#!/usr/bin/env python3
"""Tests for skipping a delegator backend that repeats another one.

Backends are consulted in order, cheapest first, and a later one normally
covers what the earlier ones lack. When two hold the same data, though, the
later one has nothing to add, and for a backend that scans a 90 GB file that
is the difference between milliseconds and minutes.

What must hold:

  * equal content_id skips, and the skip is positional -- the *later* backend
    is dropped, never the one that already answered
  * None never matches, not even another None, because it means "cannot say",
    not "holds nothing". Two silent backends must both be consulted
  * a backend that raises while identifying itself is consulted, since
    deciding what to skip must not be what breaks a query
  * skipping is about identity of data only, so it must not depend on whether
    the earlier backend actually found anything

That last point is why this is not the `final` flag: `final` claims no backend
will ever find an entry, which would also suppress backends holding entirely
different data, such as a web service that could still resolve it.

Run under pytest, or standalone: python test/db/test_delegator_redundancy.py
"""

import pandas as pd
import pytest

import rotifer.db.core
import rotifer.db.delegator


class Backend(rotifer.db.core.BaseCursor):
    """A backend that reports a fixed content_id and records being asked."""

    def __init__(self, content=None, rows=(), raises=False):
        super().__init__(progress=False)
        self._content = content
        self._rows = list(rows)
        self._raises = raises
        self.asked = 0

    def content_id(self):
        if self._raises:
            raise RuntimeError('cannot identify my data')
        return self._content

    def getids(self, obj, *args, **kwargs):
        if isinstance(obj, pd.DataFrame) and 'id' in obj.columns:
            return set(obj['id'].astype(str))
        return set()

    def fetchone(self, accessions, *args, **kwargs):
        self.asked += 1
        wanted = self.parse_ids(accessions)
        found = [x for x in self._rows if x in wanted]
        if found:
            yield pd.DataFrame({'id': found})

    def __getitem__(self, accessions, *args, **kwargs):
        self.asked += 1
        wanted = self.parse_ids(accessions)
        found = [x for x in self._rows if x in wanted]
        return pd.DataFrame({'id': found}) if found else None


class Delegator(rotifer.db.delegator.SequentialDelegatorCursor):
    """A delegator over Backend instances handed to it directly."""

    def __init__(self, backends, readers):
        self._backends = backends
        super().__init__(readers=readers)

    def reset_cursors(self):
        # Bypass the module lookup: these backends are objects, not
        # classes found in a configured module.
        self.cursors = dict(getattr(self, '_backends', {}))

    def getids(self, obj, *args, **kwargs):
        if isinstance(obj, pd.DataFrame) and 'id' in obj.columns:
            return set(obj['id'].astype(str))
        return set()


def build(*specs):
    """Build a delegator from (name, content_id, rows) triples."""
    backends = {name: Backend(content, rows) for name, content, rows in specs}
    return Delegator(backends, [name for name, _, _ in specs]), backends


def test_same_content_skips_the_later_backend():
    delegator, backends = build(('fast', 'release_1', ['a']),
                                ('slow', 'release_1', ['a', 'b']))
    list(delegator.fetchone(['a', 'b']))
    assert backends['fast'].asked == 1
    assert backends['slow'].asked == 0


def test_the_skip_is_positional():
    """Order decides which one survives: the reader list is the cost order,
    so whichever comes first is the cheap one."""
    delegator, backends = build(('slow', 'release_1', ['a']),
                                ('fast', 'release_1', ['a']))
    list(delegator.fetchone(['a']))
    assert backends['slow'].asked == 1 and backends['fast'].asked == 0


def test_different_content_consults_both():
    delegator, backends = build(('one', 'release_1', ['a']),
                                ('two', 'release_2', ['b']))
    list(delegator.fetchone(['a', 'b']))
    assert backends['one'].asked == 1 and backends['two'].asked == 1


def test_none_never_matches_none():
    """None means 'I cannot say what I hold'. Treating two of those as equal
    would skip a backend on no evidence at all."""
    delegator, backends = build(('one', None, ['a']),
                                ('two', None, ['b']))
    list(delegator.fetchone(['a', 'b']))
    assert backends['one'].asked == 1 and backends['two'].asked == 1


def test_none_is_consulted_even_after_a_known_backend():
    delegator, backends = build(('known', 'release_1', ['a']),
                                ('silent', None, ['b']))
    list(delegator.fetchone(['a', 'b']))
    assert backends['silent'].asked == 1


def test_a_backend_that_cannot_identify_itself_is_consulted():
    """Working out what to skip is an optimisation; it must never be the
    thing that loses data."""
    backends = {'first': Backend('release_1', ['a']),
                'broken': Backend(raises=True, rows=['b'])}
    delegator = Delegator(backends, ['first', 'broken'])
    list(delegator.fetchone(['a', 'b']))
    assert backends['broken'].asked == 1


def test_skipping_does_not_depend_on_the_first_backend_finding_anything():
    """Redundancy is a statement about which data a backend holds, not about
    what a particular query happened to return."""
    delegator, backends = build(('fast', 'release_1', []),
                                ('slow', 'release_1', ['a']))
    results = list(delegator.fetchone(['a']))
    assert backends['slow'].asked == 0
    assert results == []


def test_getitem_skips_too():
    delegator, backends = build(('fast', 'release_1', ['a']),
                                ('slow', 'release_1', ['a']))
    delegator['a']
    assert backends['fast'].asked == 1 and backends['slow'].asked == 0


def test_three_backends_keep_the_distinct_one():
    delegator, backends = build(('fast', 'release_1', ['a']),
                                ('slow', 'release_1', ['a']),
                                ('other', 'release_2', ['c']))
    list(delegator.fetchone(['a', 'c']))
    assert backends['slow'].asked == 0
    assert backends['other'].asked == 1


def test_redundant_reports_who_served_the_data():
    delegator, backends = build(('fast', 'release_1', ['a']),)
    assert delegator.redundant(backends['fast'], {'release_1': 'fast'}) == 'fast'
    assert delegator.redundant(backends['fast'], {}) is None


def test_base_cursor_is_never_skippable_by_default():
    """Every existing cursor inherits this, so the default has to be the safe
    one: opting in is a decision, not something a backend falls into."""
    assert rotifer.db.core.BaseCursor(progress=False).content_id() is None


if __name__ == '__main__':
    import sys
    failures = 0
    for name, test in sorted(globals().items()):
        if not name.startswith('test_') or not callable(test):
            continue
        try:
            test()
            ok, detail = True, ''
        except Exception as exc:
            ok, detail = False, f'{type(exc).__name__}: {exc}'
        failures += not ok
        print('  %-58s %s' % (name, 'ok' if ok else 'FAILED'))
        if detail:
            print('      %s' % detail)
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
