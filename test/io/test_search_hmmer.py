#!/usr/bin/env python3
"""Tests for the HMMER side of the standard search table (issue #141).

rotifer.io.search.from_pyhmmer turns hmmscan and hmmsearch output into the
standard table, and to_pyhmmer turns it back. pyhmmer_to_df also has to
return typed columns when a search finds nothing.

What has to hold, and why:

  * from_pyhmmer(hmmscan(...), 'hmmscan') passes validate and keeps the hits.
    With PYHMMER_COLUMNS, qstart/qend are filled from the model coordinates
  * with the default columns (no aln_hmm_from/aln_hmm_to), qstart/qend are
    missing, but cov is still filled for hmmscan: it only needs the target
  * for hmmscan, cov is the aligned fraction of the sequence, (end - start +
    1) / tlen; for hmmsearch it is the aligned fraction of the model,
    (qend - qstart + 1) / qlen, and qlen is the HMM length
  * to_pyhmmer(from_pyhmmer(df)) gives back df's columns and values, so the
    conversion loses nothing that the native table had
  * an empty search gives a typed table (int64 coordinates, float64 scores),
    and that converts to an empty standard table without errors

The fixture is synthetic: a random protein, mutated copies aligned with gaps,
and an HMM built from that alignment, so nothing depends on Pfam or the
network.

Run under pytest (needs pyhmmer): PYTHONPATH=lib pytest test/io/test_search_hmmer.py
"""

import numpy as np
import pandas as pd
import pytest

ph = pytest.importorskip('pyhmmer')

import rotifer.devel.beta.sequence as rdbs
import rotifer.devel.beta.hmmer as rdbh
import rotifer.io.search as st

AMINO = 'ACDEFGHIKLMNPQRSTVWY'
GAPS = '-._~'
UNGAP = str.maketrans('', '', GAPS)
SCAN = dict(workers=1, cpus_per_worker=1)


def _alignment():
    '''Eight mutated copies of a random protein, with blocks of gaps.'''
    rng = np.random.default_rng(141)
    base = rng.choice(list(AMINO), 120)
    rows = []
    for i in range(8):
        seq = base.copy()
        mutate = rng.random(120) < 0.2
        seq[mutate] = rng.choice(list(AMINO), mutate.sum())
        seq = list(seq)
        for p in range(30 + 4 * i, 36 + 4 * i):
            seq[p] = '-'
        if i % 3 == 0:
            seq[:10] = '-' * 10
        rows.append((f'seq{i}', ''.join(seq)))
    fasta = ''.join(f'>{name}\n{text}\n' for name, text in rows)
    return rdbs.sequence(fasta)


def _ungapped(seqobj):
    result = seqobj.copy()
    result.df['sequence'] = result.df.sequence.str.translate(UNGAP)
    return result


@pytest.fixture(scope='module')
def data(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('search_hmmer')
    aln = _alignment()
    hmm = tmp / 'model.hmm'
    rdbh.hmmbuild(aln, hmm_name='synthetic', save=str(hmm))
    with ph.plan7.HMMFile(str(hmm)) as hmm_file:
        ph.hmmer.hmmpress(hmm_file, str(hmm))
    query = tmp / 'query.fa'
    ungapped = _ungapped(aln).df
    query.write_text(''.join(f'>{i}\n{s}\n' for i, s in zip(ungapped.id, ungapped.sequence)))
    return dict(aln=aln, models=str(hmm), query=str(query))


@pytest.fixture(scope='module')
def hmmscan_full(data):
    return rdbh.hmmscan(_ungapped(data['aln']), models_path=data['models'],
                        columns=st.PYHMMER_COLUMNS, **SCAN)


def test_hmmscan_with_model_coordinates_is_standard(hmmscan_full):
    table = st.from_pyhmmer(hmmscan_full, 'hmmscan')
    st.validate(table)
    assert len(table) >= 8
    assert table.program.eq('hmmscan').all()
    assert table.qstart.notna().all()
    assert table.qend.notna().all()


def test_hmmscan_cov_is_the_aligned_fraction_of_the_sequence(hmmscan_full):
    table = st.from_pyhmmer(hmmscan_full, 'hmmscan')
    expected = (hmmscan_full.aln_target_to - hmmscan_full.aln_target_from + 1) \
        / hmmscan_full.aln_target_length
    assert np.allclose(table['cov'].values, expected.values)
    assert (table.qlen == hmmscan_full.aln_target_length).all()


def test_hmmscan_default_columns_leave_model_coordinates_missing(data):
    default = rdbh.hmmscan(_ungapped(data['aln']), models_path=data['models'], **SCAN)
    table = st.from_pyhmmer(default, 'hmmscan')
    st.validate(table)
    assert len(table) >= 8
    assert table.qstart.isna().all()
    assert table.qend.isna().all()
    assert table['cov'].notna().all()
    expected = (default.aln_target_to - default.aln_target_from + 1) / default.aln_target_length
    assert np.allclose(table['cov'].values, expected.values)


def test_hmmsearch_cov_is_the_aligned_fraction_of_the_model(data):
    raw = rdbh.hmmsearch(data['models'], data['query'], columns=st.PYHMMER_COLUMNS)
    table = st.from_pyhmmer(raw, program='hmmsearch')
    st.validate(table)
    assert len(table) >= 8
    assert table.program.eq('hmmsearch').all()
    assert np.allclose(table['cov'].values, ((table.qend - table.qstart + 1) / table.qlen).values)
    with ph.plan7.HMMFile(data['models']) as hmm_file:
        model_length = next(iter(hmm_file)).M
    assert (table.qlen == model_length).all()


def test_round_trip_gives_back_the_pyhmmer_table(hmmscan_full):
    back = st.to_pyhmmer(st.from_pyhmmer(hmmscan_full, 'hmmscan'))
    columns = sorted(hmmscan_full.columns)
    assert sorted(back.columns) == columns
    pd.testing.assert_frame_equal(back[columns].reset_index(drop=True),
                                  hmmscan_full[columns].reset_index(drop=True),
                                  check_dtype=False)


def test_empty_search_gives_typed_columns():
    empty = rdbh.pyhmmer_to_df([])
    assert len(empty) == 0
    assert list(empty.columns) == st.PYHMMER_COLUMNS[:12]
    assert empty['env_from'].dtype == 'int64'
    assert empty['aln_target_length'].dtype == 'int64'
    assert empty['i_evalue'].dtype == 'float64'
    assert empty['score'].dtype == 'float64'
    assert empty['aln_target_name'].dtype == object


def test_empty_search_converts_to_empty_standard_table():
    empty = rdbh.pyhmmer_to_df([])
    table = st.from_pyhmmer(empty, 'hmmscan')
    st.validate(table)
    assert len(table) == 0
    assert list(table.columns)[:len(st.COLUMNS)] == st.COLUMNS
