#!/usr/bin/env python3
"""Tests for the region slices cut by search2aln and split_by_model (issue #141).

ROTIFER coordinates are 1-based and inclusive: the region s-e of a protein
is seq[s - 1:e]. Both functions sliced the hit region as seq[s:e], which
dropped the first residue while the id still said s-e.

What has to hold, and why:

  * search2aln keeps the hit's first residue: a hit at 11-20 is exactly
    residues 11 to 20, and its id '<pid>/11-20' describes what was cut
  * split_by_model (hmmer) cuts the envelope estart-eend, both ends included
  * split_by_model (blast table) cuts hitstart-hitend, both ends included
  * a region starting at residue 1 is not empty and keeps residue 1
  * split_by_model with a blast table no longer fails on 'if df:' (the truth
    value of a DataFrame is ambiguous) when a table is given

The fixtures are synthetic: full-length sequences are runs of distinct
letters, so every residue's position can be read off directly. The database
fetch is patched and the aligner is stubbed, so nothing external is needed.

Run under pytest: PYTHONPATH=lib pytest test/devel/test_search_slices.py
"""

import pandas as pd
import pytest

import rotifer.devel.beta.sequence as rdbs
import rotifer.devel.beta.blast as rdbb

# Residue i (1-based) of FULL is FULL[i - 1]: A=1 ... Z=26, a=27 ... z=52
FULL = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'


@pytest.fixture
def database(monkeypatch):
    """Serve full-length sequences from a dict and stub the aligner.

    sequence(list_of_ids) fetches from the database. The patched class
    looks the ids up in the returned dict instead (the local database
    path is ignored). align() returns just the
    ids and sequences, unchanged, so no aligner is needed. Fill the dict in the test.
    """
    proteins = {}
    original = rdbs.sequence

    class served(original):
        def __init__(self, input=None, *args, **kwargs):
            if isinstance(input, list):
                kwargs.pop('local_database_path', None)
                found = [i for i in input if i in proteins]
                input = pd.DataFrame({'id': found, 'sequence': [proteins[i] for i in found]}) if found else None
            super().__init__(input, *args, **kwargs)

        def align(self, *args, **kwargs):
            # like the real aligner, keep only the ids and the aligned sequences
            return served(self.df[['id', 'sequence']].copy())

    monkeypatch.setattr(rdbs, 'sequence', served)
    return proteins


def blast_table(start, end, hit='p1'):
    return pd.DataFrame({'hit': [hit], 'evalue': [1e-10], 'querycoverage': [0.9],
                         'hitstart': [start], 'hitend': [end]})


def hmmer_table(start, end, protein='p1'):
    return pd.DataFrame({'sequence': [protein], 'evalue': [1e-10], 'cov': [0.9],
                         'estart': [start], 'eend': [end]})


class TestSearch2aln:
    def test_first_residue_of_the_hit_is_kept(self, database):
        database['p1'] = FULL
        s = rdbb.search2aln(blast_table(11, 20), id_with_coord=True)
        assert s.df.sequence.tolist() == ['KLMNOPQRST']   # residues 11..20
        assert s.df.id.tolist() == ['p1/11-20']

    def test_id_and_coordinate_columns_still_split(self, database):
        database['p1'] = FULL
        s = rdbb.search2aln(blast_table(11, 20))
        assert s.df.id.tolist() == ['p1']
        assert (s.df.hit_start.tolist(), s.df.hit_end.tolist()) == (['11'], ['20'])

    def test_region_starting_at_residue_1(self, database):
        database['p1'] = FULL
        s = rdbb.search2aln(blast_table(1, 10), id_with_coord=True)
        assert s.df.sequence.tolist() == ['ABCDEFGHIJ']
        assert s.df.id.tolist() == ['p1/1-10']

    def test_whole_protein_hit(self, database):
        database['p1'] = FULL
        s = rdbb.search2aln(blast_table(1, len(FULL)), id_with_coord=True)
        assert s.df.sequence.tolist() == [FULL]


class TestSplitByModelHmmer:
    def test_envelope_keeps_both_ends(self, database):
        database['p1'] = FULL
        s = rdbb.split_by_model(df=hmmer_table(11, 20), method='hmmer')
        assert s.df.sequence.tolist() == ['KLMNOPQRST']
        assert s.df.id.tolist() == ['p1/11-20']

    def test_region_starting_at_residue_1(self, database):
        database['p1'] = FULL
        s = rdbb.split_by_model(df=hmmer_table(1, 10), method='hmmer')
        assert s.df.sequence.tolist() == ['ABCDEFGHIJ']

    def test_two_domains_in_one_protein(self, database):
        database['p1'] = FULL
        df = pd.concat([hmmer_table(1, 5), hmmer_table(27, 31)], ignore_index=True)
        s = rdbb.split_by_model(df=df, method='hmmer')
        assert s.df.sequence.tolist() == ['ABCDE', 'abcde']


class TestSplitByModelBlast:
    def split(self, database, start, end):
        database['p1'] = FULL
        full = rdbs.sequence(['p1'])
        return rdbb.split_by_model(seqobj=full, df=blast_table(start, end), method='blast')

    def test_hit_keeps_both_ends(self, database):
        s = self.split(database, 11, 20)
        assert s.df.sequence.tolist() == ['KLMNOPQRST']

    def test_region_starting_at_residue_1(self, database):
        s = self.split(database, 1, 10)
        assert s.df.sequence.tolist() == ['ABCDEFGHIJ']
