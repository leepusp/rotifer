#!/usr/bin/env python3
"""Tests for the HH-suite parser in rotifer.io.hhsuite (issue #123).

The parser used to live in bin/hhsuite2table, and rotifer.io.hhsuite found
it by cutting its own path at the first "rotifer" and appending
rotifer/bin/hhsuite2table. That broke when the checkout had another name,
when a parent folder contained "rotifer", and in every installed copy,
which has no bin/ at all. The library now owns the parser and the script
imports it.

What has to hold, and why:

  * rotifer.io.hhsuite imports from a copy of the package with no bin/
    beside it, under a path containing "rotifer" more than once
  * read_hhr and parse_hhr return the hits of an .hhr file, from a file,
    a list of files or a directory
  * bin/hhsuite2table prints the library's table, with score to one decimal
    and probability as an integer, as it always has; the library itself
    returns both as float64
  * hhsuite2pandas gives the documented dtypes, on empty reports too, and
    rejects an unknown column before reading any file
  * a single path is one file, and a directory is read in name order
  * query_ali and template_ali are str for every hit

The sample file, data/sample.hhr, holds two hits. The second, last one
used to get the first hit's p-value; both are checked.

Run under pytest, or standalone: python test/io/test_hhsuite.py
"""

import io
import os
import shutil
import subprocess
import sys
import tempfile

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
LIB = os.path.join(REPO, 'lib')
SAMPLE = os.path.join(REPO, 'test', 'io', 'data', 'sample.hhr')

def _environ(pythonpath, tmp_path):
    environ = { k: v for k, v in os.environ.items() if not k.startswith('ROTIFER_') }
    environ['PYTHONPATH'] = pythonpath
    environ['ROTIFER_USER'] = os.path.join(str(tmp_path), 'user')
    # Rotifer logs to stdout: create what it would announce creating
    os.makedirs(os.path.join(environ['ROTIFER_USER'], 'etc'), exist_ok=True)
    return environ

def _hhsuite():
    if LIB not in sys.path:
        sys.path.insert(0, LIB)
    from rotifer.io import hhsuite
    return hhsuite

def test_import_from_installed_copy(tmp_path):
    site = os.path.join(str(tmp_path), 'rotifer-lab', 'rotifer', 'site-packages')
    shutil.copytree(os.path.join(LIB, 'rotifer'), os.path.join(site, 'rotifer'),
                    ignore=shutil.ignore_patterns('__pycache__'))
    out = subprocess.run(
        [sys.executable, '-c',
         'import sys\n'
         'from rotifer.io import hhsuite\n'
         f'print(len(hhsuite.read_hhr({SAMPLE!r})), "hhsuite2table" in sys.modules)'],
        env=_environ(site, tmp_path), capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.split() == ['2', 'False']

def test_read_hhr_returns_hits():
    df = _hhsuite().read_hhr(SAMPLE)
    assert list(df.columns) == _hhsuite().default_columns
    assert list(df.hit) == ['PF00001.1', 'PF00002.2']
    assert set(df['query']) == {'sp|P00001|TEST_HUMAN'}
    assert list(df.qstart) == [3, 10] and list(df.qend) == [47, 29]
    assert list(df.hstart) == [1, 5] and list(df.hend) == [45, 24]
    assert list(df.hlen) == [60, 30]
    assert list(df.probability) == [99.9, 50.0]
    assert list(df.evalue) == [1.1e-20, 0.5]
    assert list(df.identity) == [0.4, 0.25]
    assert list(df.pvalue) == [2.2e-24, 1e-4]

def test_read_hhr_accepts_directory_and_list(tmp_path):
    for name in ('a.hhr', 'b.hhr'):
        shutil.copy(SAMPLE, os.path.join(str(tmp_path), name))
    shutil.copy(SAMPLE, os.path.join(str(tmp_path), 'ignored.txt'))
    hhsuite = _hhsuite()
    assert len(hhsuite.read_hhr(str(tmp_path))) == 4
    assert len(hhsuite.read_hhr([SAMPLE, SAMPLE])) == 4

def test_parse_hhr_yields_one_table_per_file(tmp_path):
    for name in ('a.hhr', 'b.hhr'):
        shutil.copy(SAMPLE, os.path.join(str(tmp_path), name))
    tables = list(_hhsuite().parse_hhr(str(tmp_path)))
    assert [len(t) for t in tables] == [2, 2]

def test_cli_prints_formatted_library_output(tmp_path):
    """The command line equals the library table once score and probability are formatted.

    The library returns them as float64; the script prints the score with
    one decimal and the probability as an integer. Replaces a test that
    compared the CLI to a plain to_csv of the library, which no longer
    holds by design.
    """
    out = subprocess.run(
        [sys.executable, os.path.join(REPO, 'bin', 'hhsuite2table'), SAMPLE],
        env=_environ(LIB, tmp_path), capture_output=True, text=True, cwd=str(tmp_path),
        stdin=subprocess.DEVNULL,
    )
    assert out.returncode == 0, out.stderr
    df = _hhsuite().hhsuite2pandas(SAMPLE)
    df['score'] = df['score'].map(lambda x: '%-.1f' % x)
    df['probability'] = df['probability'].astype(int)
    expected = io.StringIO()
    df.to_csv(expected, sep="\t", index=None, float_format='%-.3g')
    assert out.stdout == expected.getvalue()

def test_cli_output_is_pinned(tmp_path):
    """bin/hhsuite2table prints the bytes in data/sample.hhsuite2table.tsv."""
    out = subprocess.run(
        [sys.executable, os.path.join(REPO, 'bin', 'hhsuite2table'), SAMPLE],
        env=_environ(LIB, tmp_path), capture_output=True, text=True, cwd=str(tmp_path),
        stdin=subprocess.DEVNULL,
    )
    assert out.returncode == 0, out.stderr
    with open(os.path.join(REPO, 'test', 'io', 'data', 'sample.hhsuite2table.tsv')) as handle:
        assert out.stdout == handle.read()

def test_hhsuite2pandas_values_are_pinned():
    """The library's score, probability and pvalue columns, now as floats.

    Changed on purpose: up to version 0.002 score was a string ('120.5')
    and probability an int (99); both are float64 now. bin/hhsuite2table
    keeps printing the old form.
    """
    df = _hhsuite().hhsuite2pandas(SAMPLE)
    assert list(df.score) == [120.5, 20.1]
    assert list(df.probability) == [99.9, 50.0]
    assert str(df.score.dtype) == 'float64' and str(df.probability.dtype) == 'float64'
    assert list(df.pvalue) == [2.2e-24, 1e-4]

def test_parse_hhr_single_path_is_one_table():
    """A single file name is one report, not a list of characters."""
    tables = list(_hhsuite().parse_hhr(SAMPLE))
    assert len(tables) == 1 and len(tables[0]) == 2

def test_unknown_column_raises_before_reading():
    """An unknown column raises ValueError naming it, even if the file does not exist."""
    with pytest.raises(ValueError, match='nonsense'):
        _hhsuite().hhsuite2pandas('/no/such/file.hhr', columns=['hit', 'nonsense'])

def test_empty_report_gives_typed_empty_frame(tmp_path):
    """A report with a header and no hits gives zero rows and the documented dtypes."""
    with open(SAMPLE) as handle:
        header = handle.read().split(' No Hit')[0]
    path = os.path.join(str(tmp_path), 'empty.hhr')
    with open(path, 'w') as handle:
        handle.write(header)
    df = _hhsuite().hhsuite2pandas(path)
    assert len(df) == 0
    assert list(df.columns) == _hhsuite().default_columns
    expected = {'query': 'str', 'hit': 'str', 'qstart': 'int64', 'qend': 'int64', 'hstart': 'int64',
                'hend': 'int64', 'qalilen': 'int64', 'qlen': 'int64', 'hlen': 'int64',
                'aligned_cols': 'int64', 'qcov': 'float64', 'hcov': 'float64',
                'probability': 'float64', 'evalue': 'float64', 'pvalue': 'float64',
                'score': 'float64', 'similarity': 'float64', 'sum_probs': 'float64',
                'identity': 'float64'}
    assert {c: str(df[c].dtype) for c in df.columns} == expected

def test_alignments_are_str():
    """query_ali and template_ali are str for every hit, the last one included."""
    df = _hhsuite().hhsuite2pandas(SAMPLE, columns=['hit', 'query_ali', 'template_ali'])
    assert len(df) == 2
    assert all(isinstance(v, str) for v in df.query_ali)
    assert all(isinstance(v, str) for v in df.template_ali)

def test_read_hhr_directory_is_in_name_order(tmp_path):
    """Files of a directory are read in name order, whatever order they were written in."""
    with open(SAMPLE) as handle:
        text = handle.read()
    for name, query in (('b.hhr', 'QUERY_B'), ('a.hhr', 'QUERY_A')):
        with open(os.path.join(str(tmp_path), name), 'w') as handle:
            handle.write(text.replace('sp|P00001|TEST_HUMAN', query))
    df = _hhsuite().read_hhr(str(tmp_path))
    assert list(df['query']) == ['QUERY_A', 'QUERY_A', 'QUERY_B', 'QUERY_B']


if __name__ == '__main__':
    failures = 0
    for name, test in sorted(globals().items()):
        if not name.startswith('test_') or not callable(test):
            continue
        try:
            if 'tmp_path' in test.__code__.co_varnames[:test.__code__.co_argcount]:
                with tempfile.TemporaryDirectory() as tmp:
                    test(tmp)
            else:
                test()
            ok, detail = True, ''
        except Exception as exc:
            ok, detail = False, f'{type(exc).__name__}: {exc}'
        failures += not ok
        print('  %-56s %s' % (name, 'ok' if ok else 'FAILED'))
        if detail:
            print('      %s' % detail)
    print('ALL PASS' if not failures else '%d FAILURES' % failures)
    sys.exit(1 if failures else 0)
