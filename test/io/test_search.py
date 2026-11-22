#!/usr/bin/env python3
"""Tests for rotifer.io.search's standard table (issue #141).

PSI-BLAST, PyHMMER and hmmer2table describe matched regions with different
columns. The search functions keep their own layouts; rotifer.io.search converts
them to one standard table and back.

What has to hold, and why:

  * every from_* returns the standard columns, in order, with the standard
    types, so downstream code can rely on them whatever the backend
  * converting a native table to the standard and back gives the same table,
    extra columns included, so nothing is lost by going through the standard
  * coordinates stay 1-based and inclusive, and talilen = end - start + 1
    (hmmer2table writes end - start; to_hmmer2table writes that back)
  * cov and qlen refer to the search query: the sequence for hmmscan, the
    model for hmmsearch, the query protein for PSI-BLAST
  * validate rejects tables that break the standard (missing columns,
    0-based or reversed coordinates, cov outside 0-1)

No search program is needed: the tables are small and written by hand.

Run: PYTHONPATH=lib pytest test/io/test_search.py
"""

import pandas as pd
import pytest

import rotifer.io.search as st


def _blast():
    '''psiblast()'s table: two hits, one with two regions, and an extra column.'''
    return pd.DataFrame({
        'hit':           ['WP_1', 'WP_1', 'WP_2'],
        'query':         ['Q', 'Q', 'Q'],
        'hitstart':      [5, 150, 1],
        'hitend':        [100, 245, 80],
        'evalue':        [1e-20, 1e-8, 0.003],
        'querycoverage': [0.317, 0.317, 0.267],
        'querystart':    [10, 200, 1],
        'queryend':      [104, 294, 80],
        'iteration':     [1, 2, 1],
        'bitscore':      [95.1, 60.2, 40.0],
        'length':        [300, 300, 90],
        'note':          ['a', 'b', 'c'],
    })


def _pyhmmer(renamed=False, model_coordinates=True):
    '''hmmscan's table for two sequences against two models.'''
    df = pd.DataFrame({
        'aln_target_name':   ['p1', 'p1', 'p2'],
        'aln_hmm_name':      ['PF1', 'PF2', 'PF1'],
        'i_evalue':          [1e-10, 2e-5, 3e-3],
        'c_evalue':          [1e-11, 2e-6, 3e-4],
        'score':             [50.0, 20.5, 10.0],
        'env_score':         [51.0, 21.0, 11.0],
        'aln_target_from':   [11, 60, 2],
        'aln_target_to':     [50, 90, 40],
        'aln_target_length': [100, 100, 40],
        'aln_hmm_length':    [45, 35, 45],
        'env_from':          [10, 58, 1],
        'env_to':            [52, 92, 40],
        'aln_hmm_from':      [3, 1, 5],
        'aln_hmm_to':        [42, 31, 43],
        'source':            ['m.hmm'] * 3,
    })
    if not model_coordinates:
        df = df.drop(columns=['aln_hmm_from', 'aln_hmm_to'])
    if renamed:
        df = df.rename(columns=st.PYHMMER_RENAME)
    return df


def _hmmer2table():
    return pd.DataFrame({
        'sequence':  ['p1', 'p2'],
        'model':     ['PF1', 'PF1'],
        'start':     [11, 2],
        'end':       [50, 40],
        'evalue':    [1e-10, 3e-3],
        'cov':       [0.4, 0.97],
        'qstart':    [3, 5],
        'qend':      [42, 43],
        'iteration': [1, 1],
        'score':     [50.0, 10.0],
        'talilen':   [39, 38],
        'qlen':      [100, 40],
        'estart':    [10, 1],
        'eend':      [52, 40],
        'program':   ['hmmscan', 'hmmscan'],
    })


def _same(a, b):
    '''Same columns and values, ignoring column order and dtypes.'''
    assert sorted(a.columns) == sorted(b.columns)
    pd.testing.assert_frame_equal(a[sorted(a.columns)], b[sorted(b.columns)],
                                  check_dtype=False)


# The standard table

def test_empty_has_the_standard_columns_and_types():
    df = st.empty()
    assert list(df.columns) == st.COLUMNS
    assert len(df) == 0
    assert df['start'].dtype == 'int64'
    assert df['qstart'].dtype == 'Int64'
    assert df['evalue'].dtype == 'float64'


@pytest.mark.parametrize('convert', [
    lambda: st.from_blast(_blast()),
    lambda: st.from_pyhmmer(_pyhmmer(), 'hmmscan'),
    lambda: st.from_hmmer2table(_hmmer2table()),
])
def test_conversions_return_the_standard_columns_first(convert):
    df = convert()
    assert list(df.columns[:len(st.COLUMNS)]) == st.COLUMNS
    assert df['start'].dtype == 'int64'
    assert (df['talilen'] == df['end'] - df['start'] + 1).all()


def test_validate_reorders_and_keeps_extra_columns():
    df = st.from_blast(_blast())
    shuffled = df[['note'] + st.COLUMNS[::-1]]
    out = st.validate(shuffled)
    assert list(out.columns) == st.COLUMNS + ['note']


@pytest.mark.parametrize('change, message', [
    (lambda d: d.drop(columns=['cov']), 'missing standard columns: cov'),
    (lambda d: d.assign(start=0), '1-based'),
    (lambda d: d.assign(start=d['end'] + 1), 'start must not be after end'),
    (lambda d: d.assign(cov=31.7), 'cov must be a fraction'),
    (lambda d: d.assign(start=[1, None, 1]), 'column start must be int64'),
    (lambda d: d.assign(evalue='x'), 'column evalue must be float64'),
])
def test_validate_rejects_tables_outside_the_standard(change, message):
    with pytest.raises(ValueError, match=message):
        st.validate(change(st.from_blast(_blast())))


# PSI-BLAST

def test_from_blast_maps_hits_to_sequences_and_the_query_to_the_model():
    df = st.from_blast(_blast())
    row = df.iloc[1]
    assert (row['sequence'], row['model']) == ('WP_1', 'Q')
    assert (row['start'], row['end'], row['estart'], row['eend']) == (150, 245, 150, 245)
    assert (row['qstart'], row['qend']) == (200, 294)
    assert (row['iteration'], row['score'], row['tlen']) == (2, 60.2, 300)
    assert row['talilen'] == 96
    assert pd.isna(row['qlen'])
    assert (df['program'] == 'psiblast').all()


def test_blast_round_trip_is_lossless():
    native = _blast()
    back = st.to_blast(st.from_blast(native))
    assert list(back.columns) == list(native.columns)
    pd.testing.assert_frame_equal(back, native, check_dtype=False)


# PyHMMER

@pytest.mark.parametrize('renamed', [False, True])
def test_from_pyhmmer_reads_renamed_and_raw_columns(renamed):
    df = st.from_pyhmmer(_pyhmmer(renamed=renamed), 'hmmscan')
    row = df.iloc[0]
    assert (row['sequence'], row['model']) == ('p1', 'PF1')
    assert (row['start'], row['end'], row['estart'], row['eend']) == (11, 50, 10, 52)
    assert (row['qstart'], row['qend']) == (3, 42)
    assert row['evalue'] == 1e-10
    assert row['source'] == 'm.hmm'
    assert 'c_evalue' in df.columns and 'env_score' in df.columns


def test_hmmscan_cov_and_qlen_measure_the_sequence():
    row = st.from_pyhmmer(_pyhmmer(), 'hmmscan').iloc[0]
    assert row['qlen'] == 100
    assert row['cov'] == pytest.approx(40 / 100)


def test_hmmsearch_cov_and_qlen_measure_the_model():
    row = st.from_pyhmmer(_pyhmmer(), 'hmmsearch').iloc[0]
    assert row['qlen'] == 45
    assert row['cov'] == pytest.approx(40 / 45)
    assert row['tlen'] == 100


def test_missing_model_coordinates_stay_missing():
    df = st.from_pyhmmer(_pyhmmer(model_coordinates=False), 'hmmsearch')
    assert df['qstart'].isna().all() and df['qend'].isna().all()
    assert df['cov'].isna().all()
    scan = st.from_pyhmmer(_pyhmmer(model_coordinates=False), 'hmmscan')
    assert scan['cov'].notna().all()


def test_from_pyhmmer_needs_to_know_the_program():
    with pytest.raises(ValueError, match='program must be one of'):
        st.from_pyhmmer(_pyhmmer(), 'phmmer')


@pytest.mark.parametrize('renamed', [False, True])
@pytest.mark.parametrize('model_coordinates', [False, True])
@pytest.mark.parametrize('program', st.PYHMMER_PROGRAMS)
def test_pyhmmer_round_trip_is_lossless(renamed, model_coordinates, program):
    native = _pyhmmer(renamed=renamed, model_coordinates=model_coordinates)
    _same(st.to_pyhmmer(st.from_pyhmmer(native, program)), native)


def test_empty_pyhmmer_result_converts():
    native = _pyhmmer().iloc[0:0]
    df = st.from_pyhmmer(native, 'hmmscan')
    assert list(df.columns[:len(st.COLUMNS)]) == st.COLUMNS
    assert len(df) == 0


# hmmer2table

def test_from_hmmer2table_corrects_talilen():
    df = st.from_hmmer2table(_hmmer2table())
    assert df['talilen'].tolist() == [40, 39]
    assert df['tlen'].isna().all()


def test_hmmer2table_round_trip_restores_its_talilen():
    native = _hmmer2table()
    back = st.to_hmmer2table(st.from_hmmer2table(native))
    pd.testing.assert_frame_equal(back, native, check_dtype=False)


def test_tables_from_different_backends_concatenate():
    df = pd.concat([st.from_blast(_blast()), st.from_pyhmmer(_pyhmmer(), 'hmmscan'),
                    st.from_hmmer2table(_hmmer2table())], ignore_index=True)
    out = st.validate(df)
    assert sorted(out['program'].unique()) == ['hmmscan', 'psiblast']
    assert len(out) == 8
