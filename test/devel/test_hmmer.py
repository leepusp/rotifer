#!/usr/bin/env python3
"""Tests for rotifer.devel.beta.hmmer's hmmscan on aligned input (issue #137).

Easel reads '-', '.' and '_' as gaps and '~' as missing data. Left in a
sequence to be searched, HMMER scores them as residues: on a real alignment
of 41 TraG proteins scanned against Pfam, P-loop_TraG went from 41 hits to
none, and every hit reported the alignment width as the target length.

What has to hold, and why:

  * an aligned sequence object is scanned exactly like the same sequences
    without gaps: same hits, coordinates, E-values and scores. Removing gaps
    changes nothing for unaligned input, so it is always done
  * a warning says when gaps were removed, since coordinates then refer to
    the ungapped sequences rather than to alignment columns; unaligned input
    stays silent
  * rows that are not sequences, like add_consensus's, are not scanned (they
    used to crash the digitizer with "invalid chars in sequence")
  * digitalize_seqobj(..., alignment=True) keeps the gaps: hmmbuild needs
    the columns
  * an aligned FASTA file given as file= works too (Easel's parser rejects
    gap characters outright)
  * unpressed HMM files work: HMMFile.is_pressed is a method in recent
    pyhmmer and a property in older releases, and testing the method without
    calling it was always true
  * hmmscan_linear behaves like hmmscan in all of the above

The fixture is synthetic: a random protein, mutated copies aligned with
blocks of gaps, and an HMM built from that alignment, so nothing depends on
Pfam or the network.

Run under pytest (needs pyhmmer): PYTHONPATH=lib pytest test/devel/test_hmmer.py
"""

import logging
import shutil

import numpy as np
import pytest

ph = pytest.importorskip('pyhmmer')

import rotifer.devel.beta.sequence as rdbs
import rotifer.devel.beta.hmmer as rdbh

AMINO = 'ACDEFGHIKLMNPQRSTVWY'
GAPS = '-._~'
UNGAP = str.maketrans('', '', GAPS)
KEY = ['sequence', 'model', 'estart', 'eend', 'aln_target_from', 'aln_target_to',
       'aln_target_length', 'evalue', 'score']
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


def _normalized(df):
    return df[KEY].sort_values(KEY[:4]).reset_index(drop=True)


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


@pytest.fixture
def warnings_logged():
    '''Collect warnings from the hmmer module's logger.'''
    records = []

    class Collect(logging.Handler):
        def emit(self, record):
            records.append(record)

    handler = Collect(level=logging.WARNING)
    log = logging.getLogger(rdbh.__name__)
    log.addHandler(handler)
    yield records
    log.removeHandler(handler)


def _scans():
    '''Both entry points, called the same way.'''
    return {
        'hmmscan': lambda seqs=None, **kw: rdbh.hmmscan(seqs, **SCAN, **kw),
        'hmmscan_linear': lambda seqs=None, **kw: rdbh.hmmscan_linear(seqs, cpus=1, **kw),
    }


@pytest.fixture(params=['hmmscan', 'hmmscan_linear'])
def scan(request):
    return _scans()[request.param]


def test_aligned_input_scans_like_ungapped_input(data, scan):
    aligned = scan(data['aln'], models_path=data['pressed'])
    ungapped = scan(_ungapped(data['aln']), models_path=data['pressed'])
    assert len(ungapped) >= 8, 'every sequence should hit its own model'
    assert _normalized(aligned).equals(_normalized(ungapped))
    lengths = _ungapped(data['aln']).df.set_index('id').sequence.str.len()
    assert (aligned.aln_target_length.values == lengths[aligned.sequence].values).all()


def test_removing_gaps_is_reported_once(data, scan, warnings_logged):
    scan(data['aln'], models_path=data['pressed'])
    assert len(warnings_logged) == 1
    assert 'gaps were removed' in warnings_logged[0].getMessage()


def test_unaligned_input_is_not_reported(data, scan, warnings_logged):
    scan(_ungapped(data['aln']), models_path=data['pressed'])
    assert warnings_logged == []


def test_consensus_rows_are_not_scanned(data, scan):
    with_consensus = data['aln'].add_consensus()
    assert len(with_consensus.df) > len(with_consensus)
    result = scan(with_consensus, models_path=data['pressed'])
    assert set(result.sequence) <= set(data['aln'].df.id)
    expected = scan(data['aln'], models_path=data['pressed'])
    assert _normalized(result).equals(_normalized(expected))


def test_building_an_hmm_keeps_the_gaps(data):
    msa = rdbh.digitalize_seqobj(data['aln'], alignment=True)
    width = data['aln'].get_alignment_length()
    assert all(len(row) == width for row in msa.textize().alignment)
    hmm = rdbh.hmmbuild(data['aln'], save=None)[0]
    assert hmm.M > 0


def test_digitalize_removes_gaps_for_searching(data):
    digital = rdbh.digitalize_seqobj(data['aln'])
    lengths = _ungapped(data['aln']).df.sequence.str.len().tolist()
    assert [len(d) for d in digital] == lengths


def test_aligned_fasta_file(data, scan):
    from_file = scan(file=data['aligned_fasta'], models_path=data['pressed'])
    from_seqobj = scan(data['aln'], models_path=data['pressed'])
    assert _normalized(from_file).equals(_normalized(from_seqobj))


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
