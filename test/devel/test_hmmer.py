#!/usr/bin/env python3
"""Tests for rotifer.devel.beta.hmmer's hmmscan and hmmscan_linear.

What has to hold, and why:

  * unpressed HMM files work: HMMFile.is_pressed is a method in recent
    pyhmmer and a property in older releases, and testing the method without
    calling it was always true
  * hmmscan_linear behaves like hmmscan in all of the above

The fixture is synthetic: a random protein, mutated copies aligned with
blocks of gaps, and an HMM built from that alignment, so nothing depends on
Pfam or the network.

Run under pytest (needs pyhmmer): PYTHONPATH=lib pytest test/devel/test_hmmer.py
"""

import shutil

import numpy as np
import pytest

ph = pytest.importorskip('pyhmmer')

import rotifer.devel.beta.sequence as rdbs
import rotifer.devel.beta.hmmer as rdbh

AMINO = 'ACDEFGHIKLMNPQRSTVWY'
GAPS = '-._~'
UNGAP = str.maketrans('', '', GAPS)
SCAN = dict(workers=1, cpus_per_worker=1)


def _alignment():
    '''Eight mutated copies of a random protein, aligned with gaps.

    Leading and trailing gaps, deletions (internal gap blocks), an insertion
    present in only half of the rows, and one row using '.' for its gaps.
    '''
    rng = np.random.default_rng(137)
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
        if i % 3 == 1:
            seq[-12:] = '-' * 12
        insertion = list(rng.choice(list(AMINO), 6)) if i < 4 else ['-'] * 6
        seq = seq[:80] + insertion + seq[80:]
        text = ''.join(seq)
        if i == 7:
            text = text.replace('-', '.')
        rows.append((f'seq{i}', text))
    fasta = ''.join(f'>{name}\n{text}\n' for name, text in rows)
    return rdbs.sequence(fasta), fasta


def _ungapped(seqobj):
    result = seqobj.copy()
    result.df['sequence'] = result.df.sequence.str.translate(UNGAP)
    return result


@pytest.fixture(scope='module')
def data(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('hmmer')
    aln, fasta = _alignment()
    unpressed = tmp / 'plain' / 'model.hmm'
    pressed = tmp / 'pressed' / 'model.hmm'
    unpressed.parent.mkdir()
    pressed.parent.mkdir()
    rdbh.hmmbuild(aln, hmm_name='synthetic', save=str(unpressed))
    shutil.copy(unpressed, pressed)
    with ph.plan7.HMMFile(str(pressed)) as hmm_file:
        ph.hmmer.hmmpress(hmm_file, str(pressed))
    aligned_fasta = tmp / 'aligned.fa'
    aligned_fasta.write_text(fasta)
    return dict(aln=aln, unpressed=str(unpressed), pressed=str(pressed),
                aligned_fasta=str(aligned_fasta))


def _scans():
    '''Both entry points, called the same way.'''
    return {
        'hmmscan': lambda seqs=None, **kw: rdbh.hmmscan(seqs, **SCAN, **kw),
        'hmmscan_linear': lambda seqs=None, **kw: rdbh.hmmscan_linear(seqs, cpus=1, **kw),
    }


@pytest.fixture(params=['hmmscan', 'hmmscan_linear'])
def scan(request):
    return _scans()[request.param]


@pytest.mark.parametrize('database', ['unpressed', 'pressed'])
def test_pressed_and_unpressed_databases(data, scan, database):
    result = scan(_ungapped(data['aln']), models_path=data[database])
    assert len(result) >= 8


class _Method:
    def __init__(self, value):
        self.value = value

    def is_pressed(self):
        return self.value


class _Property:
    def __init__(self, value):
        self.is_pressed = value


@pytest.mark.parametrize('value', [True, False])
def test_is_pressed_supports_method_and_property(value):
    assert rdbh._is_pressed(_Method(value)) is value
    assert rdbh._is_pressed(_Property(value)) is value
