#!/usr/bin/env python3
"""Tests for BLAST tabular input and HSP merging in rotifer.io.search (issue #141).

psiblast() reads the table that bin/blast2table makes from the text report.
read_blast_tabular() reads BLAST's own tabular output instead, and
collapse_hsps() merges its HSPs the way blast2table does by default, so the
standard table has the same rows blast2table gives.

What has to hold, and why:

  * read_blast_tabular gives one standard row per HSP, in report order, with
    the iteration of the "# Iteration: N" comment above it (1 without one,
    and 1 again for the next query's report), cov = aligned fraction of the
    query, and an empty standard table when there are no hits; a path and
    the text itself are read alike
  * two HSPs of the same (model, sequence) pair are the same region only if
    they overlap on BOTH the sequence and the model axes; overlapping on one
    is two regions (a repeat, or a domain shared by different parts)
  * a float overlap is a fraction of both HSPs' lengths, an int is a number
    of residues (blast2table's same_hsp)
  * a region already below the inclusion threshold keeps its E-value, score
    and iteration and only widens (each axis on its own, start and end
    together); one above it is replaced by an overlapping HSP with a lower
    E-value, whatever the iteration; either way it keeps its position
  * E-values are rounded as blast2table prints them before comparing, so
    a value that prints as the threshold is not below it
  * derived columns follow the final coordinates (talilen, estart/eend, cov)
  * no input order, grouping or pandas version detail leaks into the result:
    rows come out in order of first appearance, pairs interleaved as in the
    report
  * PARITY: on reports made by psiblast against a fixture database (SRC in
    1 and up to 3 rounds, ArcA in 3 rounds), to_blast(collapse_hsps(
    read_blast_tabular(fmt 7))) has exactly the rows, order and coordinates
    of blast2table's table for the same search. Differences that are
    blast2table's own are accounted for below: its E-values and bit scores
    come from the text report, which prints fewer digits (E-values to 1 or 2
    significant digits), its iteration is always 1 with the BioPerl text
    parser, and a region widened in a later round keeps the coverage of its
    first HSP

The parity tests read committed files (test/data/search, made by
make_fixtures.sh), so they don't need BLAST. A last test runs blast2table on
the committed text reports when perl and BioPerl are there.

Run: PYTHONPATH=lib pytest test/io/test_search_collapse_hsps.py
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import rotifer.io.search as st

DATA = Path(__file__).resolve().parents[1] / 'data' / 'search'
ROOT = Path(__file__).resolve().parents[2]
CASES = ['src1', 'src3', 'arca3']


def hsps(rows, qlen=100):
    '''
    Standard table of HSPs from tuples
    (model, sequence, start, end, qstart, qend, evalue, iteration, score).
    '''
    out = pd.DataFrame(rows, columns=['model', 'sequence', 'start', 'end', 'qstart', 'qend',
                                      'evalue', 'iteration', 'score'])
    out['cov'] = (out['qend'] - out['qstart'] + 1) / qlen
    out['talilen'] = out['end'] - out['start'] + 1
    out['qlen'] = qlen
    out['estart'], out['eend'] = out['start'], out['end']
    out['program'] = 'psiblast'
    out['tlen'] = 500
    out['source'] = None
    return st.validate(out)


# collapse_hsps

def test_overlap_on_one_axis_only_is_not_merged():
    # Same part of the sequence matched by two parts of the query, and the
    # reverse: two regions each time
    same_hit = hsps([('Q', 'S', 10, 60, 1, 50, 1e-10, 1, 80.0),
                     ('Q', 'S', 12, 58, 60, 100, 1e-8, 1, 70.0)])
    same_query = hsps([('Q', 'S', 10, 60, 1, 50, 1e-10, 1, 80.0),
                       ('Q', 'S', 200, 250, 2, 49, 1e-8, 1, 70.0)])
    assert len(st.collapse_hsps(same_hit)) == 2
    assert len(st.collapse_hsps(same_query)) == 2


def test_overlap_on_both_axes_is_merged():
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 60, 1, 50, 1e-10, 1, 80.0),
                                 ('Q', 'S', 12, 58, 3, 49, 1e-8, 1, 70.0)]))
    assert len(out) == 1
    assert out.loc[0, 'evalue'] == 1e-10 and out.loc[0, 'score'] == 80.0


def test_lower_evalue_replaces_a_region_that_is_not_significant_yet():
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 60, 1, 50, 0.5, 1, 20.0),
                                 ('Q', 'X', 1, 40, 1, 40, 1e-9, 1, 60.0),
                                 ('Q', 'S', 15, 70, 5, 55, 0.2, 2, 25.0)]))
    assert list(out['sequence']) == ['S', 'X']          # S keeps its place
    row = out.iloc[0]
    assert (row['start'], row['end'], row['qstart'], row['qend']) == (15, 70, 5, 55)
    assert (row['evalue'], row['iteration'], row['score']) == (0.2, 2, 25.0)


def test_higher_evalue_does_not_replace_before_threshold():
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 60, 1, 50, 0.2, 1, 25.0),
                                 ('Q', 'S', 15, 70, 5, 55, 0.5, 2, 20.0)]))
    assert len(out) == 1
    assert (out.loc[0, 'start'], out.loc[0, 'evalue'], out.loc[0, 'iteration']) == (10, 0.2, 1)


def test_significant_region_keeps_first_evalue_and_iteration_but_widens():
    # Round 1 finds the region significantly; round 2 finds it better and
    # wider on the sequence, narrower on the query
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 60, 1e-5, 1, 40.0),
                                 ('Q', 'S', 5, 60, 5, 50, 1e-9, 2, 60.0)]))
    assert len(out) == 1
    row = out.iloc[0]
    assert (row['start'], row['end']) == (5, 60)        # widened
    assert (row['qstart'], row['qend']) == (1, 60)      # query axis unchanged
    assert (row['evalue'], row['iteration'], row['score']) == (1e-5, 1, 40.0)
    assert (row['talilen'], row['estart'], row['eend']) == (56, 5, 60)
    assert row['cov'] == pytest.approx(0.60)


def test_threshold_crossed_between_rounds():
    # Not significant in round 1, replaced in round 2 by a significant HSP,
    # which then only widens in round 3
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 0.3, 1, 22.0),
                                 ('Q', 'S', 10, 52, 1, 42, 1e-4, 2, 45.0),
                                 ('Q', 'S', 10, 80, 1, 42, 1e-12, 3, 70.0)]))
    assert len(out) == 1
    row = out.iloc[0]
    assert (row['evalue'], row['iteration'], row['score']) == (1e-4, 2, 45.0)
    assert (row['start'], row['end']) == (10, 80)


def test_widening_is_decided_per_axis_and_moves_both_ends():
    # The new HSP is longer on the sequence (30 > 20) but shifted: both ends
    # are taken together, not the smallest start and the largest end
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 30, 1, 30, 1e-5, 1, 40.0),
                                 ('Q', 'S', 15, 45, 10, 35, 1e-6, 1, 41.0)]))
    assert (out.loc[0, 'start'], out.loc[0, 'end']) == (15, 45)
    assert (out.loc[0, 'qstart'], out.loc[0, 'qend']) == (1, 30)


def test_only_the_first_matching_region_is_considered():
    # The third HSP overlaps both kept regions; only the first one is widened
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 1e-5, 1, 40.0),
                                 ('Q', 'S', 45, 100, 35, 80, 1e-5, 1, 40.0),
                                 ('Q', 'S', 8, 100, 1, 80, 1e-8, 2, 50.0)]))
    assert len(out) == 2
    assert (out.loc[0, 'start'], out.loc[0, 'end'], out.loc[0, 'qend']) == (8, 100, 80)
    assert (out.loc[1, 'start'], out.loc[1, 'end'], out.loc[1, 'qend']) == (45, 100, 80)


def test_fraction_overlap_must_exceed_both_lengths():
    # 50 residues shared by a 100 and a 60 residue span on both axes
    rows = [('Q', 'S', 1, 100, 1, 100, 1e-5, 1, 40.0),
            ('Q', 'S', 51, 110, 51, 110, 1e-6, 1, 41.0)]
    assert len(st.collapse_hsps(hsps(rows), overlap=0.4)) == 1    # 50 > 40 and 50 > 24
    assert len(st.collapse_hsps(hsps(rows), overlap=0.5)) == 2    # 50 > 50 is false
    assert len(st.collapse_hsps(hsps(rows), overlap=0.8)) == 2


def test_integer_overlap_counts_residues():
    rows = [('Q', 'S', 1, 100, 1, 100, 1e-5, 1, 40.0),
            ('Q', 'S', 98, 200, 98, 200, 1e-6, 1, 41.0)]         # 3 residues shared
    table = hsps(rows, qlen=300)
    assert len(st.collapse_hsps(table, overlap=2)) == 1
    assert len(st.collapse_hsps(table, overlap=3)) == 2           # 3 > 3 is false
    assert len(st.collapse_hsps(table, overlap=0)) == 1
    assert len(st.collapse_hsps(table, overlap=0.4)) == 2         # a float is a fraction


def test_evalues_are_rounded_like_blast2table_before_comparing():
    # 0.0099996 prints as 1.000e-02: not below the threshold, so a lower
    # E-value replaces it instead of just widening it
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 0.0099996, 1, 30.0),
                                 ('Q', 'S', 10, 60, 1, 50, 0.0099, 2, 31.0)]))
    assert (out.loc[0, 'evalue'], out.loc[0, 'iteration'], out.loc[0, 'end']) == (0.0099, 2, 60)
    # the value returned is the unrounded one of the HSP kept
    kept = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 0.123456, 1, 30.04)]))
    assert (kept.loc[0, 'evalue'], kept.loc[0, 'score']) == (0.123456, 30.04)
    assert list(st.blast2table_evalue([0, 12.34, 0.0099996, 1.0])) == [0, 12.3, 0.01, 1.0]
    # a tie is not "lower"
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 0.5, 1, 30.0),
                                 ('Q', 'S', 10, 50, 1, 40, 0.5, 2, 31.0)]))
    assert out.loc[0, 'iteration'] == 1


def test_cov_follows_the_widened_query_span():
    out = st.collapse_hsps(hsps([('Q', 'S', 10, 50, 1, 40, 1e-5, 1, 40.0),
                                 ('Q', 'S', 10, 60, 1, 90, 1e-6, 2, 41.0)], qlen=120))
    assert (out.loc[0, 'qstart'], out.loc[0, 'qend']) == (1, 90)
    assert out.loc[0, 'cov'] == pytest.approx(0.75)
    assert out['cov'].between(0, 1).all()


def test_pairs_are_independent_and_rows_keep_first_seen_order():
    out = st.collapse_hsps(hsps([('Q1', 'A', 1, 50, 1, 50, 1e-5, 1, 40.0),
                                 ('Q2', 'A', 1, 50, 1, 50, 1e-5, 1, 40.0),
                                 ('Q1', 'B', 1, 50, 1, 50, 1e-5, 1, 40.0),
                                 ('Q1', 'A', 1, 50, 1, 50, 1e-9, 2, 60.0),
                                 ('Q2', 'A', 90, 140, 1, 50, 1e-5, 1, 40.0)]))
    assert list(zip(out['model'], out['sequence'], out['start'])) == [
        ('Q1', 'A', 1), ('Q2', 'A', 1), ('Q1', 'B', 1), ('Q2', 'A', 90)]


def test_empty_and_untouched_tables():
    out = st.collapse_hsps(st.empty())
    assert list(out.columns) == st.COLUMNS and out.empty
    one = hsps([('Q', 'S', 10, 50, 1, 40, 1e-5, 1, 40.0)])
    pd.testing.assert_frame_equal(st.collapse_hsps(one), one)
    assert len(one) == 1                                  # the input is not modified


def test_extra_columns_follow_the_row_that_gave_the_evalue():
    df = hsps([('Q', 'S', 10, 50, 1, 40, 0.5, 1, 20.0), ('Q', 'S', 12, 52, 2, 42, 0.1, 2, 25.0)])
    df['note'] = ['first', 'second']
    out = st.collapse_hsps(df)
    assert list(out['note']) == ['second'] and list(out.columns[-1:]) == ['note']


def test_missing_query_coordinates_are_refused():
    df = hsps([('Q', 'S', 10, 50, 1, 40, 0.5, 1, 20.0)])
    df['qstart'] = pd.array([pd.NA], dtype='Int64')
    with pytest.raises(ValueError, match='qstart'):
        st.collapse_hsps(df)


# read_blast_tabular

REPORT = '''# PSIBLAST 2.17.0+
# Iteration: 1
# Query: Q1
# Database: db
# Fields: query acc.ver, subject acc.ver, s. start, s. end, evalue, q. start, q. end, bit score, subject length, query length
# 2 hits found
Q1\tA\t1\t50\t1e-10\t1\t50\t80.5\t300\t100
Q1\tB\t5\t40\t0.003\t10\t45\t30.1\t200\t100
# PSIBLAST 2.17.0+
# Iteration: 2
# Query: Q1
# 1 hits found
Q1\tA\t1\t60\t1e-20\t1\t60\t99.9\t300\t100
# BLAST processed 1 queries
'''


def test_read_blast_tabular_text_and_iterations():
    df = st.read_blast_tabular(REPORT, source='db')
    assert list(df.columns) == st.COLUMNS
    assert list(df['iteration']) == [1, 1, 2]
    assert list(df['sequence']) == ['A', 'B', 'A'] and set(df['model']) == {'Q1'}
    first = df.iloc[0]
    assert (first['start'], first['end'], first['estart'], first['eend']) == (1, 50, 1, 50)
    assert (first['qstart'], first['qend'], first['qlen'], first['tlen']) == (1, 50, 100, 300)
    assert first['evalue'] == 1e-10 and first['score'] == 80.5 and first['talilen'] == 50
    assert first['cov'] == pytest.approx(0.5)
    assert (first['program'], first['source']) == ('psiblast', 'db')
    pd.testing.assert_frame_equal(df, st.validate(df))


def test_read_blast_tabular_without_iteration_comments_and_multiple_queries():
    # blastp style rows have no "# Iteration"; a new header starts again at 1
    text = ('# BLASTP 2.17.0+\n# Query: Q1\nQ1\tA\t1\t50\t1e-10\t1\t50\t80.5\t300\t100\n'
            '# BLASTP 2.17.0+\n# Query: Q2\nQ2\tA\t1\t20\t1e-3\t1\t20\t30.5\t300\t50\n')
    df = st.read_blast_tabular(text)
    assert list(df['model']) == ['Q1', 'Q2'] and list(df['iteration']) == [1, 1]
    assert df.loc[1, 'cov'] == pytest.approx(0.4)
    # psiblast with several queries: rounds restart for each query
    text = REPORT + REPORT.replace('Q1', 'Q2').replace('# Iteration: 2', '# Iteration: 3')
    assert list(st.read_blast_tabular(text)['iteration']) == [1, 1, 2, 1, 1, 3]


def test_read_blast_tabular_zero_hits():
    text = '# PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: Q1\n# 0 hits found\n# BLAST processed 1 queries\n'
    df = st.read_blast_tabular(text)
    assert df.empty and list(df.columns) == st.COLUMNS
    pd.testing.assert_frame_equal(df, st.empty())
    assert st.collapse_hsps(df).empty
    assert st.read_blast_tabular('').empty


def test_read_blast_tabular_path_and_text_are_the_same(tmp_path):
    path = tmp_path / 'report.tsv'
    path.write_text(REPORT)
    expected = st.read_blast_tabular(REPORT)
    pd.testing.assert_frame_equal(st.read_blast_tabular(path), expected)
    pd.testing.assert_frame_equal(st.read_blast_tabular(str(path)), expected)
    gz = tmp_path / 'report.tsv.gz'
    with gzip.open(gz, 'wt') as handle:
        handle.write(REPORT)
    pd.testing.assert_frame_equal(st.read_blast_tabular(gz), expected)
    with open(path) as handle:
        pd.testing.assert_frame_equal(st.read_blast_tabular(handle), expected)
    # a report with a single row and no comments is still text
    pd.testing.assert_frame_equal(st.read_blast_tabular('Q1\tA\t1\t50\t1e-10\t1\t50\t80.5\t300\t100'),
                                  expected.iloc[[0]].reset_index(drop=True).assign(iteration=1))


def test_read_blast_tabular_rejects_other_layouts():
    with pytest.raises(ValueError, match='fields'):
        st.read_blast_tabular('Q1\tA\t1\t50\t1e-10\n')


def test_blast_outfmt_names_the_fields_read():
    assert st.BLAST_OUTFMT == '7 ' + st.BLAST_TABULAR_FIELDS
    assert len(st.BLAST_TABULAR_FIELDS.split()) == 10


# Parity with blast2table, on reports made by psiblast against the fixtures

def _fixture(case):
    hsp = st.read_blast_tabular(DATA / f'{case}.fmt7.tsv', source='fixture')
    expected = pd.read_csv(DATA / f'{case}.b2t.tsv', sep='\t', names=st.BLAST_COLUMNS)
    return hsp, st.collapse_hsps(hsp), expected


@pytest.mark.parametrize('case,hsp_rows,rows', [('src1', 453, 453), ('src3', 913, 465),
                                                ('arca3', 1071, 366)])
def test_parity_rows_order_and_coordinates(case, hsp_rows, rows):
    hsp, merged, expected = _fixture(case)
    assert (len(hsp), len(merged), len(expected)) == (hsp_rows, rows, rows)
    got = st.to_blast(merged)
    for column in ['hit', 'query', 'hitstart', 'hitend', 'querystart', 'queryend', 'length']:
        assert list(got[column]) == list(expected[column]), column


@pytest.mark.parametrize('case', CASES)
def test_parity_evalues_and_bit_scores(case):
    _, merged, expected = _fixture(case)
    got = st.to_blast(merged)
    mine, theirs = got['evalue'].to_numpy(), expected['evalue'].to_numpy()
    assert ((mine == 0) == (theirs == 0)).all()
    # blast2table reads the text report: E-values below 1 have 1 or 2
    # significant digits (so differ from the tabular ones by up to half a
    # unit of the last printed digit), above 1 one decimal
    big = mine > 1
    assert (mine[big] == theirs[big]).all()
    small = (mine > 0) & ~big
    unit = 10.0 ** np.floor(np.log10(mine[small]))
    assert (np.abs(mine[small] - theirs[small]) <= 0.5 * unit * (1 + 1e-9)).all()
    # bit scores: blast2table prints no decimals above 99
    rounded = np.where(got['bitscore'] > 99, got['bitscore'].round(0), got['bitscore'].round(1))
    assert np.allclose(rounded, expected['bitscore'], atol=1e-9)


@pytest.mark.parametrize('case', CASES)
def test_parity_coverage(case):
    hsp, merged, expected = _fixture(case)
    qlen = hsp['qlen'].iloc[0]
    got = st.to_blast(merged)
    own = (expected['queryend'] - expected['querystart'] + 1) / qlen
    consistent = (expected['querycoverage'] - own).abs() <= 1e-3 + 1e-9
    # always the aligned fraction of the final coordinates
    assert ((got['querycoverage'] - own).abs() < 1e-12).all()
    # equal to 3 decimals wherever blast2table's own coverage matches its coordinates
    assert ((got['querycoverage'] - expected['querycoverage'])[consistent].abs() <= 1e-3 + 1e-9).all()
    # elsewhere blast2table kept the coverage of the region's first HSP, from
    # before it was widened by a later round
    assert (expected['querycoverage'][~consistent] < got['querycoverage'][~consistent]).all()
    assert case == 'src1' or (~consistent).any()
    assert case != 'src1' or consistent.all()


@pytest.mark.parametrize('case', CASES)
def test_parity_iterations(case):
    # blast2table (BioPerl's text parser) labels every row 1, so rows are
    # checked against the HSPs they came from instead
    hsp, merged, expected = _fixture(case)
    got = st.to_blast(merged)
    key = ['model', 'sequence', 'score', 'iteration']
    assert len(merged.merge(hsp[key], on=key)) >= len(merged)
    if case == 'src1':
        assert (got['iteration'] == expected['iteration']).all() and (got['iteration'] == 1).all()
    else:
        assert merged['iteration'].max() == hsp['iteration'].max() >= 2
        assert (merged['iteration'] >= 1).all()


@pytest.mark.parametrize('case', ['src3', 'arca3'])
def test_fixtures_exercise_the_threshold_and_widening(case):
    hsp, merged, _ = _fixture(case)
    first = hsp.drop_duplicates(['model', 'sequence'])
    # regions first seen above the threshold and settled in a later round
    late = merged[merged['iteration'] > 1]
    assert (late['evalue'] >= 0.01).any() and (late['evalue'] < 0.01).any()
    # regions widened after the first round
    widened = merged.merge(first, on=['model', 'sequence'], suffixes=('', '_first'))
    assert ((widened['end'] - widened['start']) > (widened['end_first'] - widened['start_first'])).any()


@pytest.mark.parametrize('case', CASES)
def test_collapsed_fixtures_follow_the_standard(case):
    hsp, merged, _ = _fixture(case)
    pd.testing.assert_frame_equal(merged, st.validate(merged))
    assert merged['cov'].between(0, 1).all()
    assert (merged['talilen'] == merged['end'] - merged['start'] + 1).all()
    # no two regions of a pair overlap on both axes like HSPs did
    rows = hsp.drop_duplicates(['model', 'sequence', 'start', 'end', 'qstart', 'qend'])
    assert len(merged) <= len(rows)


def _blast2table_works():
    perl = shutil.which('perl')
    if not perl or not (ROOT / 'bin' / 'blast2table').exists():
        return False
    return subprocess.run([perl, '-MBio::SearchIO', '-e', '1'], capture_output=True).returncode == 0


@pytest.mark.skipif(not _blast2table_works(), reason='needs perl with BioPerl')
@pytest.mark.parametrize('case', CASES)
def test_committed_blast2table_tables_match_the_text_reports(case, tmp_path):
    report = tmp_path / f'{case}.txt'
    with gzip.open(DATA / f'{case}.txt.gz', 'rb') as src:
        report.write_bytes(src.read())
    done = subprocess.run([shutil.which('perl'), str(ROOT / 'bin' / 'blast2table'), str(report)],
                          stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300)
    assert done.returncode == 0, done.stderr
    assert done.stdout == (DATA / f'{case}.b2t.tsv').read_text()


def _merge(rows, evalues, threshold=0.01, overlap=0.4, fraction=True):
    """Run _merge_hsps on one group; rows are (hit_start, hit_end, query_start, query_end)."""
    cols = list(zip(*rows))
    return st._merge_hsps([range(len(rows))], *cols, evalues, threshold, overlap, fraction)


def test_merge_hsps_widens_a_significant_region():
    """A region below the threshold keeps its source and widens on the axes where the new HSP is longer."""
    out = _merge([(10, 20, 1, 10), (5, 30, 3, 8)], [1e-9, 1e-20])
    assert out == [[0, 0, 5, 30, 1, 10]]


def test_merge_hsps_replaces_a_weak_region_with_a_better_hsp():
    """Above the threshold, an overlapping HSP with a lower E-value replaces the region entirely."""
    out = _merge([(10, 20, 1, 10), (12, 22, 2, 11)], [0.5, 0.2])
    assert out == [[0, 1, 12, 22, 2, 11]]


def test_merge_hsps_drops_a_worse_overlapping_hsp():
    """Above the threshold, an overlapping HSP that is not better (ties included) is dropped."""
    assert _merge([(10, 20, 1, 10), (12, 22, 2, 11)], [0.2, 0.5]) == [[0, 0, 10, 20, 1, 10]]
    assert _merge([(10, 20, 1, 10), (12, 22, 2, 11)], [0.2, 0.2]) == [[0, 0, 10, 20, 1, 10]]


def test_merge_hsps_fraction_versus_residues():
    """The same pair overlaps by 3 residues: a fraction of 0.4 (of 11) rejects it, 2 residues accepts it."""
    rows = [(1, 11, 1, 11), (9, 19, 9, 19)]
    assert len(_merge(rows, [1e-9, 1e-9], overlap=0.4, fraction=True)) == 2
    assert len(_merge(rows, [1e-9, 1e-9], overlap=2, fraction=False)) == 1
    assert len(_merge(rows, [1e-9, 1e-9], overlap=3, fraction=False)) == 2


def test_merge_hsps_sorts_regions_by_first_position_across_groups():
    """Groups are independent and the output is ordered by the position of each region's first HSP."""
    rows = [(1, 10, 1, 10), (1, 10, 1, 10), (50, 60, 1, 10)]
    out = st._merge_hsps([[1], [0, 2]], *zip(*rows), [1e-9] * 3, 0.01, 0.4, True)
    assert [region[0] for region in out] == [0, 1, 2]
