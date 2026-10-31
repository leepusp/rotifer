#!/usr/bin/env python3
"""Tests for Rotifer's logger (rotifer.core.logger) and loadConfig's errors.

Rotifer's command line tools write their data to stdout and are chained
in pipes, so anything else written there corrupts the data downstream.

What has to hold, and why:

  * log messages go to stderr, never stdout
  * a message is written once, however many modules, or calls to
    get_logger, are involved: the handler sits on the "rotifer" logger
    and module loggers reach it by propagation
  * loggers outside Rotifer's namespace, like a script's "__main__",
    still get a handler, once
  * rotifer.logger.setLevel still governs every module's logger, as
    rotifer.db.ncbi relies on
  * loadConfig logs a YAML file it cannot parse and goes on, instead of
    raising TypeError (logger.error(..., file=...)) or NameError

Each case runs a fresh interpreter, with HOME pointing to an empty
directory so that importing rotifer has something to log.

Run under pytest, or standalone: python test/core/test_logger.py
"""

import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
LIB = os.path.join(REPO, 'lib')

def _run(code, tmp_path):
    '''Run code in a fresh interpreter and return its stdout and stderr.'''
    environ = { k: v for k, v in os.environ.items() if not k.startswith('ROTIFER_') }
    environ['PYTHONPATH'] = LIB
    environ['HOME'] = os.path.join(str(tmp_path), 'home')
    os.makedirs(environ['HOME'], exist_ok=True)
    out = subprocess.run([sys.executable, '-c', code], env=environ,
                         capture_output=True, text=True, cwd=str(tmp_path))
    assert out.returncode == 0, out.stderr
    return out.stdout, out.stderr

def test_import_warnings_go_to_stderr(tmp_path):
    stdout, stderr = _run('import rotifer; print("DATA")', tmp_path)
    assert stdout == 'DATA\n'
    assert 'Creating directory' in stderr

def test_module_message_written_once(tmp_path):
    stdout, stderr = _run(
        'import rotifer.core.functions as rcf\n'
        'rcf.logger.error("logged once")', tmp_path)
    assert stdout == ''
    assert stderr.count('logged once') == 1

def test_repeated_get_logger_adds_no_handlers(tmp_path):
    stdout, stderr = _run(
        'from rotifer.core import logger as rcl\n'
        'for _ in range(3): log = rcl.get_logger("rotifer.demo")\n'
        'log.warning("said once")', tmp_path)
    assert stdout == ''
    assert stderr.count('said once') == 1

def test_foreign_logger_gets_one_handler(tmp_path):
    stdout, stderr = _run(
        'from rotifer.core import logger as rcl\n'
        'for _ in range(3): log = rcl.get_logger("__main__")\n'
        'log.warning("from a script")', tmp_path)
    assert stdout == ''
    assert stderr.count('from a script') == 1
    assert 'WARNING __main__ from a script' in stderr

def test_rotifer_level_governs_modules(tmp_path):
    stdout, stderr = _run(
        'import rotifer\n'
        'log = rotifer.logging.getLogger("rotifer.db.demo")\n'
        'log.info("hidden")\n'
        'rotifer.logger.setLevel(rotifer.logging.INFO)\n'
        'log.info("shown")', tmp_path)
    assert stdout == ''
    assert 'hidden' not in stderr
    assert stderr.count('shown') == 1

def test_loadconfig_logs_unparsable_yaml(tmp_path):
    conf = os.path.join(str(tmp_path), 'conf')
    os.makedirs(conf)
    with open(os.path.join(conf, 'multi.yml'), 'w') as fh:
        fh.write('a: 1\n---\nb: 2\n')
    stdout, stderr = _run(
        'from rotifer.core.functions import loadConfig\n'
        f'print(loadConfig(":multi", user_path={conf!r}, system_path={conf!r}, defaults={{"x": 1}}))',
        tmp_path)
    assert stdout == "{'x': 1}\n"
    assert stderr.count('Error while parsing file') == 2  # system and user paths
    assert 'expected a single document' in stderr

def test_loadconfig_logs_unparsable_yaml_file_path(tmp_path):
    path = os.path.join(str(tmp_path), 'multi.yml')
    with open(path, 'w') as fh:
        fh.write('a: 1\n---\nb: 2\n')
    stdout, stderr = _run(
        'from rotifer.core.functions import loadConfig\n'
        f'print(loadConfig({path!r}, defaults={{"x": 1}}))', tmp_path)
    assert stdout == "{'x': 1}\n"
    assert f'Error while parsing file {path}' in stderr


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
