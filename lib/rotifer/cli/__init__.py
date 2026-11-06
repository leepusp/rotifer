__doc__ = """
Command line entry points for ROTIFER tools.

Every tool exists twice: as a library function (for example
``rotifer.io.hmmer.read_hmmer``) and as a shell command of the same name
(``hmmer2table``). This package holds the command half. Each module
``rotifer.cli.<command>`` defines:

``add_arguments(parser)``
    declare the command's options on an ``argparse.ArgumentParser``
``run(args)``
    do the work, write the result to standard output, return an exit code
``main(argv=None)``
    build the parser, parse ``argv`` and call ``run``

Commands are kept thin on purpose: they parse options and call the
library, so the command line framework can change without touching the
code that does the work. To keep start-up fast, command modules import
only the standard library at load time and never import pandas.

The executables in ``bin/`` are symbolic links to ``bin/_rotifer_cli``,
which runs ``rotifer.cli.<name of the link>``. Each module can also be run
with ``python -m rotifer.cli.<command>``.
"""

import sys
from rotifer.io import STDIN, open_input

def input_files(files):
    """
    Return the list of inputs to process: the files given or, when
    none was given, standard input.
    """
    return list(files) if files else [STDIN]

def write_rows(rows, header=None, out=None, sep="\t"):
    """
    Write rows (sequences of strings) as a delimited text table.
    """
    out = out or sys.stdout
    if header:
        out.write(sep.join(header) + "\n")
    for row in rows:
        out.write(sep.join(row) + "\n")

def write_padded(rows, out=None, sep="\t"):
    """
    Write rows as a delimited table whose columns are padded with spaces
    to a common width, for reading on screen.
    """
    out = out or sys.stdout
    rows = [list(row) for row in rows]
    widths = [max(len(row[i]) for row in rows if i < len(row))
              for i in range(max((len(row) for row in rows), default=0))]
    for row in rows:
        cells = row + [""] * (len(widths) - len(row))
        out.write(sep.join(cell.ljust(width) for cell, width in zip(cells, widths)) + "\n")

def key_value(text):
    """
    ``argparse`` type for ``key=value`` options.
    """
    key, sep, value = text.partition("=")
    if not sep or not key:
        raise ValueError(f"expected key=value, got {text!r}")
    return key, value
