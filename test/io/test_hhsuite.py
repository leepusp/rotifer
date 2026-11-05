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
  * bin/hhsuite2table prints what the library returns, so the command line
    and the library cannot drift apart

The sample file, data/sample.hhr, holds two hits. Only the first hit's
p-value is checked: the parser assigns the last hit of a file the p-value
of the hit before it, which is a separate bug.

Run under pytest, or standalone: python test/io/test_hhsuite.py
"""

import io
import os
import shutil
import subprocess
import sys
import tempfile

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
    assert list(df.probability) == [99, 50]
    assert list(df.evalue) == [1.1e-20, 0.5]
    assert list(df.identity) == [0.4, 0.25]
    assert df.pvalue[0] == 2.2e-24

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

def test_cli_prints_library_output(tmp_path):
    out = subprocess.run(
        [sys.executable, os.path.join(REPO, 'bin', 'hhsuite2table'), SAMPLE],
        env=_environ(LIB, tmp_path), capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert out.returncode == 0, out.stderr
    expected = io.StringIO()
    _hhsuite().hhsuite2pandas(SAMPLE).to_csv(expected, sep="\t", index=None, float_format='%-.3g')
    assert out.stdout == expected.getvalue()


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
