#!/usr/bin/env python3
"""Tests for rotifer.io.dali, which turns a DALI server TSV into a star alignment.

What has to hold, and why:

  * dali_tsv lays every hit on the query's columns: residues a hit inserts
    between two aligned blocks open gap columns in all other rows, so every
    row (sequences and DSSP annotations) has the same length
  * ``unaligned`` decides how inserted residues are written: lowercase
    (``lower``), as they are (``keep``), as gaps (``gap``), or dropped with
    their columns (``remove``)
  * DSSP rows (``ss_from:<id>``) come first and follow their sequence's
    columns; the query's DSSP row is missing when the file has no
    ``# DSSP:`` line
  * an invalid ``unaligned`` and a table left empty by ``max_hits`` raise
    ValueError
  * importing rotifer.io.dali loads no rotifer.devel or rotifer.db module and
    prints nothing; rotifer.devel.beta.sequence is imported when dali_tsv runs
  * _make_star_row writes aligned residues uppercase in the query's columns,
    uncovered query positions as '-', and insertions formatted per mode and
    padded to insert_sizes, so the row is query_len + sum(insert_sizes) long

data/toy_dali.tsv is a hand-written table: a 9-residue query, a hit s1
with two blocks and a 5-residue insertion between them plus one trailing
residue, and a hit s2 with one block. data/toy_dali.<mode>.tsv hold the
alignment written for each ``unaligned`` mode.

Run under pytest: PYTHONPATH=lib pytest test/io/test_dali.py
"""

import os
import subprocess
import sys

import pandas as pd
import pytest

DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
TOY = os.path.join(DATA, 'toy_dali.tsv')
LIB = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'lib'))
COLUMNS = ['id', 'sequence', 'type', 'description']


def _dali():
    """Import rotifer.io.dali (skipped when its dependencies are missing)."""
    return pytest.importorskip('rotifer.io.dali')


@pytest.mark.parametrize('mode', ['lower', 'gap', 'remove', 'keep'])
def test_alignment_matches_pinned_output(mode):
    """Each unaligned mode writes the rows in data/toy_dali.<mode>.tsv."""
    aln = _dali().dali_tsv(TOY, unaligned=mode)
    expected = pd.read_csv(os.path.join(DATA, f'toy_dali.{mode}.tsv'), sep='\t', dtype=str)
    got = aln.df[COLUMNS].reset_index(drop=True).astype(str)
    pd.testing.assert_frame_equal(got, expected, check_dtype=False)


def test_rows_have_the_same_length():
    """Sequences and DSSP rows share the alignment's columns."""
    aln = _dali().dali_tsv(TOY)
    assert aln.df.sequence.str.len().nunique() == 1


def test_query_dssp_row_needs_a_dssp_line(tmp_path):
    """Without '# DSSP:', only the hits get an ss_from row."""
    path = tmp_path / 'no_dssp.tsv'
    with open(TOY) as handle:
        path.write_text(''.join(line for line in handle if not line.startswith('# DSSP:')))
    aln = _dali().dali_tsv(str(path))
    assert list(aln.df.id) == ['ss_from:s1', 'ss_from:s2', 'q1', 's1', 's2']


def test_invalid_mode_and_empty_table_raise():
    """A bad unaligned mode and max_hits=0 are errors."""
    dali = _dali()
    with pytest.raises(ValueError):
        dali.dali_tsv(TOY, unaligned='bogus')
    with pytest.raises(ValueError):
        dali.dali_tsv(TOY, max_hits=0)


def test_import_loads_no_devel_or_db_module(tmp_path):
    """Importing rotifer.io.dali pulls in no rotifer.devel or rotifer.db module and prints nothing."""
    environ = {k: v for k, v in os.environ.items() if not k.startswith('ROTIFER_')}
    environ['PYTHONPATH'] = LIB
    environ['ROTIFER_USER'] = os.path.join(str(tmp_path), 'user')
    # Rotifer logs to stdout: create what it would announce creating
    os.makedirs(os.path.join(environ['ROTIFER_USER'], 'etc'), exist_ok=True)
    code = (
        'import sys, rotifer.io.dali\n'
        'print(sorted(m for m in sys.modules if m.startswith(("rotifer.devel", "rotifer.db"))))'
    )
    out = subprocess.run(
        [sys.executable, '-c', code],
        env=environ, capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == '[]'
    assert out.stderr == ''


def test_make_star_row_on_a_hand_built_row():
    """Blocks go in uppercase, gaps where no block covers, insertions per mode and padded."""
    dali = _dali()
    # query of 6; blocks q2-3 <- s2-3 and q5-6 <- s6-7 (0-based starts 1,4 / 1,5);
    # the hit has 'a' before the first block, 'XY' between the blocks and 'FGz' after
    row = pd.Series({'qstarts': '1,4', 'sstarts': '1,5', 'lengths': '2,2'})
    residues = 'aBcXYdeFGz'
    sizes = [0, 1, 0, 0, 3, 0, 4]
    expected = {
        'keep': '-aBC-XY-DEFGz-',
        'lower': '-aBC-xy-DEfgz-',
        'gap': '--BC----DE----',
    }
    for mode, want in expected.items():
        got = dali._make_star_row(row, residues, 6, sizes, unaligned=mode)
        assert got == want
        assert len(got) == 6 + sum(sizes)
    assert dali._make_star_row(row, residues, 6, [0] * 7, unaligned='remove') == '-BC-DE'
    with pytest.raises(ValueError):
        dali._make_star_row(row, residues, 6, sizes, unaligned='bogus')
