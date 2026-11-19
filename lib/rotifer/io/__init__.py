__doc__ = """
Parsers for the output of external programs.

Each submodule reads one kind of report and is imported on demand, for
example ``from rotifer.io import hmmer``:

================================ ===========================================
:mod:`rotifer.io.blast`          BLAST+ pairwise text reports (``blast2table``)
:mod:`rotifer.io.hmmer`          HMMER 3 reports (``hmmer2table``)
:mod:`rotifer.io.phobius`        Phobius predictions (``phobius2table``)
:mod:`rotifer.io.hhsuite`        HH-suite ``.hhr`` reports (``hhsuite2table``)
:mod:`rotifer.io.dali`           DALI server TSV files, as alignments
:mod:`rotifer.io.dalidat`        DaliLite ``.dat`` files (``dat2fa``)
:mod:`rotifer.io.search`         one standard table for search results
================================ ===========================================

The table parsers share a shape: ``parse_<tool>(handle)`` (or
``iter_blast``) yields rows of strings, which the command line tools
write as they come, and ``read_<tool>(source)`` loads the same rows into
a pandas DataFrame.

This module holds the helpers they share: :func:`open_input` and
:func:`each_source` open files, standard input and compressed files;
:func:`split_fields`, :func:`written_as_fraction` and :func:`strip_suffix`
reproduce Perl behaviour the former tools relied on.
"""

import os
import re
import sys
from contextlib import contextmanager

STDIN = "-"
"""File name that stands for standard input."""

@contextmanager
def open_input(path):
    """
    Open a text input for reading, as a context manager.

    Parameters
    ----------
    path : str or os.PathLike
        File name, or ``-`` (:data:`STDIN`) for standard input. Names ending
        in ``.gz`` or ``.bz2`` are decompressed on the fly (see
        ``rotifer.core.functions.open_compressed``).

    Yields
    ------
    file object
        An open text handle. It is closed when the ``with`` block ends,
        except standard input, which is left open.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.

    Examples
    --------
    >>> from rotifer.io import open_input
    >>> with open_input("test/data/golden/inputs/phobius.out") as handle:
    ...     handle.readline().split()
    ['ID', 'NP_414542.1']
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
    Open each report in ``source`` in turn.

    Parameters
    ----------
    source : str, os.PathLike, list of them, or file object
        A file name, a list of file names, or an already open file. A single
        string or path is one file, not a list of characters. Names are
        opened with :func:`open_input`, so ``-`` is standard input and
        compressed files are read as text.

    Yields
    ------
    tuple of (str, file object)
        The report's name and an open handle. Files this function opens are
        closed when the next one is requested; an open file passed in is
        yielded once, under its ``name`` attribute (``-`` if it has none),
        and is not closed.

    Raises
    ------
    FileNotFoundError
        If a file does not exist (raised when that file is reached).

    Examples
    --------
    >>> from rotifer.io import each_source
    >>> [name for name, handle in each_source(["test/io/data/sample.hhr"])]
    ['test/io/data/sample.hhr']
    """
    if hasattr(source, "read"):
        yield getattr(source, "name", STDIN), source
        return
    for path in [source] if isinstance(source, (str, os.PathLike)) else source:
        path = os.fspath(path)
        with open_input(path) as handle:
            yield path, handle

def split_fields(text, pattern):
    """
    Split text on a regular expression and drop trailing empty fields.

    This is what Perl's ``split`` does, and the former Perl tools relied on
    it: ``"a\\t\\tb\\t\\t"`` gives three fields, not five.

    Parameters
    ----------
    text : str
        The text to split.
    pattern : str
        Regular expression that separates fields.

    Returns
    -------
    list of str
        The fields, without the empty ones at the end. Empty fields in the
        middle are kept. An empty text gives an empty list.

    Examples
    --------
    >>> from rotifer.io import split_fields
    >>> split_fields("a\\t\\tb\\t\\t", "\\t")
    ['a', '', 'b']
    """
    fields = re.split(pattern, text)
    while fields and fields[-1] == "":
        fields.pop()
    return fields

def written_as_fraction(value):
    """
    Tell whether a number is written as a fraction.

    The former Perl tools read an overlap written with a decimal point
    (``0.4``, ``.5``) or as ``<digits>e-<digits>`` (``1e-2``) as a fraction
    of a length, and any other number as a count of residues. The test is on
    the text of the number, so ``1`` and ``1.0`` mean different things.

    Parameters
    ----------
    value : str, int or float
        The number, as typed by the user or as Python writes it (a float is
        converted with ``str``, so ``0.4`` is a fraction and ``1.0`` too).

    Returns
    -------
    bool
        True for a fraction, False for a residue count.

    Notes
    -----
    The rules are Perl's, kept for identical output: ``5.`` (no digit
    after the point), ``1E-2`` (capital E) and ``1e+3`` are residue counts.

    Examples
    --------
    >>> from rotifer.io import written_as_fraction
    >>> [written_as_fraction(v) for v in ["0.4", "1", 1.0, "1e-2", "5."]]
    [True, False, True, True, False]
    """
    return bool(re.search(r"\d*\.\d+", str(value)) or re.search(r"\d+e-\d+", str(value)))

def strip_suffix(path, suffixes):
    """
    Return a file's base name without the first suffix that matches.

    Parameters
    ----------
    path : str
        File name; only its base name is used.
    suffixes : list of str
        Regular expressions, tried in order, each anchored at the end of
        the name (as in Perl's ``File::Basename::fileparse``).

    Returns
    -------
    str
        The base name, without the first matching suffix, or unchanged if
        none matches.

    Examples
    --------
    >>> from rotifer.io import strip_suffix
    >>> strip_suffix("/data/PF00001.hmmsearch", [r"\\.hmm\\w+"])
    'PF00001'
    """
    name = os.path.basename(path)
    for suffix in suffixes:
        match = re.search(f"(?:{suffix})$", name)
        if match:
            return name[:match.start()]
    return name
