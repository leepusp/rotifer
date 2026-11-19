"""
Helpers for golden tests: compare a command's output with the output the
former Perl tool produced for the same command line.

Golden data lives in test/data/golden (see its README.md). For a command
``<command>``, each ``<command>/<case>.cmd`` holds a command line and
``<case>.out`` the expected standard output. Command lines run from
test/data/golden with bin/ first in PATH, so they exercise the executable
users and scripts call, not only the library behind it.
"""

import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
LIB = os.path.join(REPO, "lib")
GOLDEN = os.path.join(REPO, "test", "data", "golden")
INPUTS = os.path.join(GOLDEN, "inputs")

if LIB not in sys.path:
    sys.path.insert(0, LIB)

def cases(command):
    """
    List the golden cases of a command as (name, command line) pairs.
    """
    directory = os.path.join(GOLDEN, command)
    return [
        (name[:-4], read(os.path.join(command, name)).strip())
        for name in sorted(os.listdir(directory)) if name.endswith(".cmd")
    ]

def read(path):
    """
    Read a file under test/data/golden.
    """
    with open(os.path.join(GOLDEN, path)) as handle:
        return handle.read()

def expected(command, case):
    """
    The recorded output of a golden case.
    """
    return read(os.path.join(command, f"{case}.out"))

def environ(tmp_path):
    """
    Environment for running commands: this checkout's bin/ and the current
    Python first in PATH, and a throw-away rotifer user directory.
    """
    env = {k: v for k, v in os.environ.items() if not k.startswith("ROTIFER_")}
    env["PATH"] = os.pathsep.join([
        os.path.join(REPO, "bin"), os.path.dirname(sys.executable), env.get("PATH", ""),
    ])
    env["ROTIFER_USER"] = os.path.join(str(tmp_path), "user")
    os.makedirs(os.path.join(env["ROTIFER_USER"], "etc"), exist_ok=True)
    env.pop("PYTHONPATH", None)
    return env

def run(cmdline, tmp_path):
    """
    Run a golden command line and return its standard output.
    """
    result = subprocess.run(
        ["bash", "-c", cmdline], cwd=GOLDEN, env=environ(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
    )
    assert result.returncode == 0, f"{cmdline}\n{result.stderr}"
    return result.stdout

def check(command, case, cmdline, tmp_path):
    """
    Assert that a golden command line still prints the recorded output.
    """
    assert run(cmdline, tmp_path) == expected(command, case), cmdline
