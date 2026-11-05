#!/usr/bin/env python3
"""Tests for the path resolver in rotifer.core (issue #122).

rotifer.core.config tells the rest of Rotifer where to find configuration
and data files and where the user directory is. It used to assume a git
checkout, which broke pip/conda installs, containers and shared servers.

What has to hold, and why:

  * in a checkout, 'base' is the repository root and the bundled files are
    the ones under lib/rotifer/_data, which etc/rotifer and share/rotifer
    point to, so scripts using the old paths keep working
  * in an installed copy, with no checkout around it, the bundled files are
    found inside the package itself
  * ROTIFER_BASE overrides both, for sites keeping their own configuration
  * an unwritable user directory (read-only $HOME on HPC or in containers)
    falls back to a private temporary directory instead of failing
  * ~/.rotifer/lib/python is only added to sys.path if it exists
  * loadConfig resolves its default locations when called, not when the
    module is imported

The configuration is computed at import time, so each case runs a fresh
interpreter.

Run under pytest, or standalone: python test/core/test_paths.py
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
LIB = os.path.join(REPO, 'lib')

def _probe(pythonpath, env={}, code='print(json.dumps(rotifer.core.config))'):
    '''Import rotifer in a fresh interpreter and return what code prints, parsed.'''
    environ = { k: v for k, v in os.environ.items() if not k.startswith('ROTIFER_') }
    environ.update(env)
    environ['PYTHONPATH'] = pythonpath
    out = subprocess.run(
        [sys.executable, '-c', 'import json, sys, rotifer, rotifer.core\n' + code],
        env=environ, capture_output=True, text=True, cwd=tempfile.gettempdir(),
    )
    assert out.returncode == 0, out.stderr
    return json.loads(out.stdout.strip().splitlines()[-1])

def _installed_copy(tmp_path):
    '''Copy the package somewhere with no checkout around it.'''
    site = os.path.join(str(tmp_path), 'site-packages')
    shutil.copytree(os.path.join(LIB, 'rotifer'), os.path.join(site, 'rotifer'),
                    ignore=shutil.ignore_patterns('__pycache__'))
    return site

def test_checkout_layout(tmp_path):
    config = _probe(LIB, {'ROTIFER_USER': os.path.join(str(tmp_path), 'user')})
    assert config['base'] == REPO
    assert config['baseConfig'] == os.path.join(LIB, 'rotifer', '_data', 'etc')
    assert config['baseDataDirectory'] == os.path.join(LIB, 'rotifer', '_data', 'share')
    assert os.path.realpath(os.path.join(REPO, 'etc', 'rotifer')) == config['baseConfig']
    assert os.path.realpath(os.path.join(REPO, 'share', 'rotifer')) == config['baseDataDirectory']

def test_installed_layout(tmp_path):
    site = _installed_copy(tmp_path)
    config = _probe(site, {'ROTIFER_USER': os.path.join(str(tmp_path), 'user')})
    package = os.path.join(os.path.realpath(site), 'rotifer')
    assert config['package'] == package
    assert config['baseConfig'] == os.path.join(package, '_data', 'etc')
    assert config['baseDataDirectory'] == os.path.join(package, '_data', 'share')
    assert os.path.isfile(os.path.join(config['baseConfig'], 'cli', 'core.yaml'))
    assert os.path.isdir(os.path.join(config['baseDataDirectory'], 'taxonomy'))

def test_installed_layout_loads_bundled_config(tmp_path):
    site = _installed_copy(tmp_path)
    loaded = _probe(site, {'ROTIFER_USER': os.path.join(str(tmp_path), 'user')},
        'from rotifer.core.functions import loadConfig\n'
        'print(json.dumps(sorted(loadConfig(":cli.core"))))')
    assert loaded, 'cli/core.yaml was not found in the installed package'

def test_base_override(tmp_path):
    base = os.path.join(str(tmp_path), 'site')
    os.makedirs(os.path.join(base, 'etc', 'rotifer'))
    os.makedirs(os.path.join(base, 'share', 'rotifer'))
    config = _probe(LIB, {'ROTIFER_BASE': base, 'ROTIFER_USER': os.path.join(str(tmp_path), 'user')})
    assert config['base'] == os.path.realpath(base)
    assert config['baseConfig'] == os.path.join(os.path.realpath(base), 'etc', 'rotifer')
    assert config['baseDataDirectory'] == os.path.join(os.path.realpath(base), 'share', 'rotifer')

def test_user_override(tmp_path):
    user = os.path.join(str(tmp_path), 'user')
    config = _probe(LIB, {'ROTIFER_USER': user})
    assert config['user'] == os.path.realpath(user)
    assert config['userConfig'] == os.path.join(os.path.realpath(user), 'etc')
    assert config['userDataDirectory'] == os.path.join(os.path.realpath(user), 'share')

def test_unwritable_user_falls_back(tmp_path):
    if hasattr(os, 'geteuid') and os.geteuid() == 0:
        return  # root can write anywhere
    readonly = os.path.join(str(tmp_path), 'readonly')
    os.makedirs(readonly)
    os.chmod(readonly, 0o500)
    try:
        config = _probe(LIB, {'ROTIFER_USER': os.path.join(readonly, '.rotifer')})
    finally:
        os.chmod(readonly, 0o700)
    assert not config['user'].startswith(os.path.realpath(readonly))
    assert config['user'].startswith(tempfile.gettempdir())
    assert not os.path.exists(os.path.join(readonly, '.rotifer'))

def test_missing_user_lib_not_on_path(tmp_path):
    user = os.path.join(str(tmp_path), 'user')
    libdir = os.path.join(os.path.realpath(user), 'lib', 'python')
    assert not _probe(LIB, {'ROTIFER_USER': user}, f'print(json.dumps({libdir!r} in sys.path))')
    os.makedirs(libdir)
    assert _probe(LIB, {'ROTIFER_USER': user}, f'print(json.dumps({libdir!r} in sys.path))')

def test_cache_override(tmp_path):
    cache = os.path.join(str(tmp_path), 'cache')
    os.makedirs(cache)
    found = _probe(LIB, {'ROTIFER_CACHE': cache, 'ROTIFER_USER': os.path.join(str(tmp_path), 'user')},
        'print(json.dumps(rotifer.config["cache"]))')
    assert found == cache

def test_loadconfig_defaults_resolved_at_call_time(tmp_path):
    user = os.path.join(str(tmp_path), 'later')
    os.makedirs(user)
    with open(os.path.join(user, 'probe.yml'), 'w') as fh:
        fh.write('answer: 42\n')
    found = _probe(LIB, {'ROTIFER_USER': os.path.join(str(tmp_path), 'user')},
        'from rotifer.core.functions import loadConfig\n'
        f'rotifer.core.config["userConfig"] = {user!r}\n'
        'print(json.dumps(loadConfig(":probe")))')
    assert found == {'answer': 42}


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
