__doc__ = """
Parsers for the output of external programs.

Submodules are imported on demand, e.g. ``from rotifer.io import hhsuite``.
"""

import os
import re
import sys
from contextlib import contextmanager

STDIN = "-"

@contextmanager
def open_input(path):
    """
    Open a text input for reading.

    Parameters
    ----------
    path : str
        File name, or ``-`` for standard input. Names ending in ``.gz``
        or ``.bz2`` are decompressed on the fly (see
        ``rotifer.core.functions.open_compressed``).
    """
    if path == STDIN:
        yield sys.stdin
        return
    from rotifer.core.functions import open_compressed
    handle = open_compressed(path)
    try:
        yield handle
    finally:
        handle.close()

def each_source(source):
    """
    Yield ``(name, handle)`` for each report in ``source``: a file name,
    a list of file names, or an already open file.
    """
    if hasattr(source, "read"):
        yield getattr(source, "name", STDIN), source
        return
    for path in [source] if isinstance(source, str) else source:
        with open_input(path) as handle:
            yield path, handle

def split_fields(text, pattern):
    """
    Split ``text`` on a regular expression, dropping trailing empty
    fields (like Perl's ``split``, which the former tools relied on).
    """
    fields = re.split(pattern, text)
    while fields and fields[-1] == "":
        fields.pop()
    return fields

def written_as_fraction(value):
    """
    True if a number is written as a fraction: with a decimal point
    (``0.4``) or in scientific notation (``1e-2``). The former Perl tools
    used this to tell a fractional overlap from a number of residues.
    """
    return bool(re.search(r"\d*\.\d+", str(value)) or re.search(r"\d+e-\d+", str(value)))

def strip_suffix(path, suffixes):
    """
    Base name of ``path`` without the first matching suffix (regular
    expressions, as in Perl's ``File::Basename::fileparse``).
    """
    name = os.path.basename(path)
    for suffix in suffixes:
        match = re.search(f"(?:{suffix})$", name)
        if match:
            return name[:match.start()]
    return name
