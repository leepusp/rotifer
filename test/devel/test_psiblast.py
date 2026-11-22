#!/usr/bin/env python3
"""Tests for rotifer.devel.beta.blast.psiblast (issue #141).

psiblast() used to write the pairwise text report and parse it with the
blast2table script. It now runs psiblast once with -outfmt 11 and formats the
archive twice (tabular and text) with blast_formatter.

What has to hold, and why:

  * rows, coordinates and order of the table are those blast2table gave on
    the same search (the committed .b2t.tsv); iteration, querycoverage and
    evalue are the correct ones, i.e. what rotifer.io.search makes of the
    committed tabular output
  * the text report is returned too (it has the query id and, with more than
    one round, "Results from round"); report=False returns the same table
    and None, without an archive or a formatter
  * a missing database fails at once and names the database (blast_formatter
    would otherwise try the NCBI servers for seconds); a list of databases,
    which is what config['blastdb'] holds, works
  * a search without hits gives an empty table with BLAST_COLUMNS
  * something that is not a sequence object gives None

The database is built from test/data/search/fixture.fa.gz, with the options
of make_fixtures.sh.

Run under pytest (needs BLAST+ on PATH): PYTHONPATH=lib pytest test/devel/test_psiblast.py
"""

import gzip
import os
import shutil
import subprocess
import time

import pandas as pd
import pytest

if not all(shutil.which(p) for p in ('psiblast', 'blast_formatter', 'makeblastdb', 'blastdbcmd')):
    pytest.skip('BLAST+ is not on PATH', allow_module_level=True)

from rotifer.devel.beta import blast
import rotifer.io.search as st
from rotifer.devel.beta import sequence as rdbs

DATA = os.path.join(os.path.dirname(__file__), '..', 'data', 'search')
CASES = {'src1': ('src', 1), 'src3': ('src', 3), 'arca3': ('arca', 3)}
KEYS = ['hit', 'query', 'hitstart', 'hitend', 'querystart', 'queryend', 'length']
RIGHT = ['iteration', 'querycoverage', 'evalue']


@pytest.fixture(scope='session')
def db(tmp_path_factory):
    d = tmp_path_factory.mktemp('psiblastdb')
    fa = d / 'fixture.fa'
    with gzip.open(os.path.join(DATA, 'fixture.fa.gz'), 'rb') as f:
        fa.write_bytes(f.read())
    subprocess.run(['makeblastdb', '-in', str(fa), '-dbtype', 'prot', '-parse_seqids',
                    '-out', str(d / 'fx')], check=True, stdout=subprocess.DEVNULL)
    return str(d / 'fx')


def query(name):
    return rdbs.sequence.from_file(os.path.join(DATA, f'{name}.fa'))


def run(case, db, **kw):
    q, it = CASES[case]
    return blast.psiblast(query(q), db=db, cpu=2, aln=False, iterations=it, evalue=10, **kw)


@pytest.mark.parametrize('case', list(CASES))
def test_table_matches_blast2table_and_correct_values(case, db):
    t, text = run(case, db)
    old = pd.read_csv(os.path.join(DATA, f'{case}.b2t.tsv'), sep='\t', names=st.BLAST_COLUMNS)
    assert list(t.columns) == st.BLAST_COLUMNS
    pd.testing.assert_frame_equal(t[KEYS].reset_index(drop=True), old[KEYS], check_dtype=False)
    right = st.to_blast(st.collapse_hsps(st.read_blast_tabular(os.path.join(DATA, f'{case}.fmt7.tsv'))))
    pd.testing.assert_frame_equal(t[RIGHT].reset_index(drop=True), right[RIGHT].reset_index(drop=True),
                                  check_dtype=False)


def test_text_report(db):
    t, text = run('src3', db)
    assert 'Results from round' in text
    assert 'P12931' in text
    t1, text1 = run('src1', db)
    assert 'P12931' in text1


def test_report_false(db):
    t, text = run('src3', db)
    t2, none = run('src3', db, report=False)
    assert none is None
    pd.testing.assert_frame_equal(t, t2)


def test_missing_db_fails_fast(tmp_path):
    start = time.time()
    with pytest.raises((FileNotFoundError, ValueError)) as e:
        blast.psiblast(query('src'), db=str(tmp_path / 'nodb'), cpu=2, aln=False)
    assert time.time() - start < 2
    assert 'nodb' in str(e.value)


def test_db_as_list(db):
    t, _ = run('src1', [db])
    t2, _ = run('src1', db)
    pd.testing.assert_frame_equal(t, t2)


def test_no_hits_gives_empty_table(db):
    q = rdbs.sequence.from_string('>nohit\nWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n')
    t, text = blast.psiblast(q, db=db, cpu=2, aln=False)
    assert len(t) == 0
    assert list(t.columns) == st.BLAST_COLUMNS


def test_not_a_sequence_returns_none(db):
    assert blast.psiblast('P12931', db=db) is None


def test_no_leftovers_in_cwd(db, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run('src1', db)
    assert list(tmp_path.iterdir()) == []
