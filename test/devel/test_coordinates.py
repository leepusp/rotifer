#!/usr/bin/env python3
"""Tests for residue coordinates in sequence and envelope helpers (issue #142).

ROTIFER coordinates are 1-based and inclusive: residue i is seq[i - 1] and
the region s-e is seq[s - 1:e]. Several helpers mixed this with Python's
0-based offsets and returned wrong coordinates, lengths and ids without
any warning.

What has to hold, and why:

  * fragment_coordinates returns 1-based inclusive start/end, and missing
    values plus a warning when a fragment can't be located, instead of
    find()'s -1
  * gaps_to_numbers(adjust_coordinates=True), used by the HTML and LaTeX
    alignment views, numbers each row from its first to its last residue
    in the full-length sequence (the trailing number was one past the end)
  * add_cordinates_to_aln's start is 1-based like its end (it was find()'s
    0-based offset), so extend_aln adds exactly n_terminal residues
  * extend_aln's length column is the extended sequence's length, not the
    full-length protein's
  * extend_aln only warns about missing coordinates for aligned proteins,
    not for every other protein in full_length_seqobj

The fixtures are synthetic: full-length sequences are runs of distinct
letters, so every residue's position can be read off directly, and the
full-length sequences are passed in (or the database fetch is patched),
so nothing needs a sequence database.

Run under pytest: PYTHONPATH=lib pytest test/devel/test_coordinates.py
"""

import logging

import pandas as pd
import pytest

import rotifer.devel.beta.sequence as rdbs
import rotifer.devel.alpha.gian_func as rdag

# Residue i (1-based) of FULL is FULL[i - 1]: A=1 ... Z=26, a=27 ... z=52
FULL = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'


def seqobj(ids, sequences):
    """A sequence object built from plain lists, without any database."""
    return rdbs.sequence(pd.DataFrame({'id': ids, 'sequence': sequences}))


@pytest.fixture
def database(monkeypatch):
    """Serve full-length sequences from a dict instead of a database.

    sequence(list_of_ids) fetches from the database. The patched class
    looks the ids up in the returned dict instead (unknown ids are left
    out, as a failed fetch would). Fill the dict in the test.
    """
    proteins = {}
    original = rdbs.sequence

    class served(original):
        def __init__(self, input=None, *args, **kwargs):
            if isinstance(input, list):
                found = [i for i in input if i in proteins]
                input = pd.DataFrame({'id': found, 'sequence': [proteins[i] for i in found]}) if found else None
            super().__init__(input, *args, **kwargs)

    monkeypatch.setattr(rdbs, 'sequence', served)
    return proteins


# --- fragment_coordinates ----------------------------------------------------

class TestFragmentCoordinates:
    def test_coordinates_are_1_based_and_inclusive(self):
        c = rdbs.fragment_coordinates(pd.Series(['KLMNOPQRST', 'ABCDEFGHIJ']), pd.Series([FULL, FULL]))
        assert c.start.tolist() == [11, 1]
        assert c.end.tolist() == [20, 10]

    def test_gaps_are_ignored(self):
        c = rdbs.fragment_coordinates(pd.Series(['--KLM-NO--']), pd.Series([FULL]))
        assert (c.start[0], c.end[0]) == (11, 15)

    def test_slice_recovers_the_fragment(self):
        c = rdbs.fragment_coordinates(pd.Series(['xyz']), pd.Series([FULL]))
        assert FULL[c.start[0] - 1:c.end[0]] == 'xyz'

    def test_fragment_not_found_is_missing_and_warned(self, caplog):
        caplog.set_level(logging.WARNING)
        c = rdbs.fragment_coordinates(pd.Series(['XXXX'], index=['p1']), pd.Series([FULL], index=['p1']))
        assert c.start.isna().all() and c.end.isna().all()
        assert 'p1' in caplog.text

    def test_missing_full_sequence_is_missing_and_warned(self, caplog):
        caplog.set_level(logging.WARNING)
        c = rdbs.fragment_coordinates(pd.Series(['ABC'], index=['p9']), pd.Series([None], index=['p9']))
        assert c.start.isna().all()
        assert 'p9' in caplog.text

    def test_no_warning_when_everything_is_found(self, caplog):
        caplog.set_level(logging.WARNING)
        rdbs.fragment_coordinates(pd.Series(['ABC']), pd.Series([FULL]))
        assert caplog.text == ''


# --- gaps_to_numbers (#7) ----------------------------------------------------

# m1 is residues 11-16 of FULL and is the model whose gap columns are hidden;
# p1 is residues 5-12 of its own protein and has a 2-residue insertion (QQ)
GAPPED = (['m1', 'p1'], ['KLM--NOP', 'KLMQQNOP'])
P1_FULL = 'AAAAKLMQQNOPCCCC'


class TestGapsToNumbers:
    def test_rows_are_numbered_from_first_to_last_residue(self, database):
        database.update(m1=FULL, p1=P1_FULL)
        _, numbers, _, _ = seqobj(*GAPPED).gaps_to_numbers('m1', adjust_coordinates=True)
        assert numbers.iloc[:, 0].tolist() == [11, 5]   # first residue
        assert numbers.iloc[:, 1].tolist() == [0, 2]    # hidden insertion
        assert numbers.iloc[:, -1].tolist() == [16, 12] # last residue

    def test_compact_view_shows_the_last_residue(self, database):
        database.update(m1=FULL, p1=P1_FULL)
        view = seqobj(*GAPPED).compact_residue_df(50, remove_gaps='m1', adjust_coordinates=True)
        assert view.loc['m1'].iloc[[0, -1]].str.strip().tolist() == ['11', '16']
        assert view.loc['p1'].iloc[[0, -1]].str.strip().tolist() == ['5', '12']

    def test_without_full_length_sequences_rows_are_numbered_on_their_own(self, database):
        _, numbers, _, _ = seqobj(*GAPPED).gaps_to_numbers('m1', adjust_coordinates=True)
        assert numbers.iloc[:, 0].tolist() == [1, 1]
        assert numbers.iloc[:, -1].tolist() == [6, 8]


# --- add_cordinates_to_aln and extend_aln (#1, #2) ---------------------------

def ungapped(aln):
    """Map each id to its sequence without gaps."""
    return dict(zip(aln.df.id, aln.df.sequence.str.replace('-', '')))


@pytest.fixture
def aligned():
    """p1 is residues 11-20 and p2 residues 2-10 of FULL; p3 is a bystander."""
    full = seqobj(['p1', 'p2', 'p3'], [FULL] * 3)
    aln = seqobj(['p1', 'p2'], ['KLMNOPQRST', '-BCDEFGHIJ'])
    return aln, full


class TestAddCoordinates:
    def test_start_and_end_are_1_based_and_inclusive(self, aligned):
        aln, full = aligned
        c = rdag.add_cordinates_to_aln(aln, full_length_seqobj=full).df.set_index('id')
        assert c.start.tolist() == [11, 2]
        assert c.end.tolist() == [20, 10]
        assert c.C_term.tolist() == [32, 42]
        assert c.full_length.tolist() == [52, 52]

    def test_fragment_not_found_is_missing_not_minus_one(self, aligned, caplog):
        caplog.set_level(logging.WARNING)
        _, full = aligned
        c = rdag.add_cordinates_to_aln(seqobj(['p1'], ['XXXX']), full_length_seqobj=full)
        assert c.df.start.isna().all()
        assert 'p1' in caplog.text


class TestExtendAln:
    def test_extends_by_exactly_the_requested_residues(self, aligned):
        aln, full = aligned
        e = ungapped(rdag.extend_aln(aln, n_terminal=3, c_terminal=3, full_length_seqobj=full))
        assert e['p1'] == FULL[11 - 3 - 1:20 + 3]   # HIJ ... TUVW
        assert e['p2'] == FULL[0:10 + 3]            # capped at residue 1

    def test_sequences_not_extended_keep_their_start(self, aligned):
        aln, full = aligned
        e = ungapped(rdag.extend_aln(aln, n_terminal=3, c_terminal=0, n_term_list=['p1'], full_length_seqobj=full))
        assert e['p2'] == FULL[1:10]                # still starts at residue 2 (B)

    def test_length_is_the_extended_length(self, aligned):
        aln, full = aligned
        e = rdag.extend_aln(aln, n_terminal=3, c_terminal=3, full_length_seqobj=full).df.set_index('id')
        assert e.length.to_dict() == {'p1': 16, 'p2': 13}

    def test_no_warning_for_proteins_outside_the_alignment(self, aligned, caplog, capsys):
        caplog.set_level(logging.WARNING)
        aln, full = aligned                         # p3 is only in full
        e = rdag.extend_aln(aln, n_terminal=3, c_terminal=3, full_length_seqobj=full)
        assert sorted(e.df.id) == ['p1', 'p2']
        reported = caplog.text + capsys.readouterr().out   # it used to be printed
        assert 'missing the coordinates' not in reported

    def test_aligned_protein_without_coordinates_is_warned_and_dropped(self, aligned, caplog):
        caplog.set_level(logging.WARNING)
        _, full = aligned
        aln = seqobj(['p1', 'p2'], ['KLMNOPQRST', 'XXXXXXXXXX'])
        e = rdag.extend_aln(aln, n_terminal=3, c_terminal=3, full_length_seqobj=full)
        assert e.df.id.tolist() == ['p1']
        assert 'missing the coordinates  p2' in caplog.text

    def test_two_fragments_of_one_protein_are_both_kept(self):
        full = seqobj(['p1'], [FULL])
        aln = seqobj(['p1', 'p1'], ['KLMNOPQRST', 'efghijklmn'])
        e = rdag.extend_aln(aln, n_terminal=0, c_terminal=0, full_length_seqobj=full)
        assert sorted(e.df.sequence.str.replace('-', '')) == ['KLMNOPQRST', 'efghijklmn']
