#!/usr/bin/env python3
"""Tests that the UniProt delegator skips a backend repeating another one.

These duplicate, against `rotifer.db.uniprot.BaseUniProtDelegatorCursor`, what
test/db/test_delegator_redundancy.py checks against the generic delegator --
and that duplication is the point.

`BaseUniProtDelegatorCursor` overrides `fetchone` for a reason of its own: the
generic version hands every result to every writer, including results a writer
just returned as a reader. Because it is a separate loop, it did not inherit
the skip when that was added to the generic delegator, and the generic tests
kept passing while a real query still scanned the whole 90 GB mirror. Anything
that reimplements the delegation loop has to be tested on its own loop.

Run under pytest, or standalone: python test/db/uniprot/test_delegator_skip.py
"""

import pandas as pd
import pytest

import rotifer.db.core
import rotifer.db.methods
from rotifer.db import uniprot


class Backend(rotifer.db.methods.MappingCursor, rotifer.db.core.BaseCursor):
    """A mapping backend with a fixed content_id that counts being asked."""

    def __init__(self, content=None, accessions=()):
        super().__init__(progress=False)
        self._content = content
        self._accessions = list(accessions)
        self.asked = 0

    def content_id(self):
        return self._content

    def fetchone(self, accessions, source=None, target=None, *args, **kwargs):
        self.asked += 1
        wanted = self.parse_ids(accessions)
        rows = [{'source': a, 'source_type': self.UNIPROTKB, 'accession': a,
                 'target': f'{a}_ref', 'target_type': 'RefSeq'}
                for a in self._accessions if a in wanted]
        if rows:
            yield pd.DataFrame(rows, columns=self.columns)


class Delegator(uniprot.BaseUniProtDelegatorCursor):
    """A UniProt delegator over backends handed to it directly."""

    def __init__(self, backends, readers):
        self._backends = backends
        super().__init__(readers=readers)

    def reset_cursors(self):
        self.cursors = dict(getattr(self, '_backends', {}))


def build(*specs):
    backends = {name: Backend(content, accs) for name, content, accs in specs}
    return Delegator(backends, [name for name, _, _ in specs]), backends


def test_same_content_skips_the_expensive_backend():
    """The regression: clickhouse and the mirror loaded from one file, so the
    mirror has nothing to add and must not be scanned."""
    delegator, backends = build(('clickhouse', 'idmapping.dat:90:1', ['P00750']),
                                ('mirror', 'idmapping.dat:90:1', ['P00750']))
    delegator.fetchall(['P00750'])
    assert backends['clickhouse'].asked == 1
    assert backends['mirror'].asked == 0


def test_the_skip_holds_when_nothing_was_found():
    """The case that actually cost 78 seconds: an accession absent from both.
    Skipping must follow from the backends holding the same data, not from the
    first one having answered."""
    delegator, backends = build(('clickhouse', 'idmapping.dat:90:1', []),
                                ('mirror', 'idmapping.dat:90:1', []))
    result = delegator.fetchall(['ABSENT'])
    assert backends['mirror'].asked == 0
    assert result.empty and list(result.columns) == [
        'source', 'source_type', 'accession', 'target', 'target_type']


def test_a_differing_release_still_consults_the_mirror():
    """A table loaded from another release is not a substitute: the mirror may
    hold entries it does not."""
    delegator, backends = build(('clickhouse', 'idmapping.dat:90:1', []),
                                ('mirror', 'idmapping.dat:91:2', ['P00750']))
    frame = delegator.fetchall(['P00750'])
    assert backends['mirror'].asked == 1
    assert frame['source'].tolist() == ['P00750']


def test_an_unrecorded_load_consults_the_mirror():
    """content_id is None for a load that was never recorded, e.g. one that was
    interrupted. Nothing may be skipped on the strength of that."""
    delegator, backends = build(('clickhouse', None, []),
                                ('mirror', 'idmapping.dat:90:1', ['P00750']))
    delegator.fetchall(['P00750'])
    assert backends['mirror'].asked == 1


def test_webapi_is_still_consulted_after_both_local_backends():
    """Why this is not the `final` flag: a web service holds different data, so
    it must still be asked for what the local copies lack."""
    delegator, backends = build(('clickhouse', 'idmapping.dat:90:1', []),
                                ('mirror', 'idmapping.dat:90:1', []),
                                ('webapi', None, ['P00750']))
    frame = delegator.fetchall(['P00750'])
    assert backends['mirror'].asked == 0
    assert backends['webapi'].asked == 1
    assert frame['source'].tolist() == ['P00750']


def test_getitem_skips_too():
    delegator, backends = build(('clickhouse', 'idmapping.dat:90:1', ['P00750']),
                                ('mirror', 'idmapping.dat:90:1', ['P00750']))
    delegator['P00750']
    assert backends['mirror'].asked == 0


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
        print('  %-52s %s' % (name, 'ok' if ok else 'FAILED'))
        if detail:
            print('      %s' % detail)
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
