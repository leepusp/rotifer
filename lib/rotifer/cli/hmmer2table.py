__doc__ = """
hmmer2table: convert HMMER 3 reports to a domain table.

Library equivalent: ``rotifer.io.hmmer.read_hmmer``.
"""

import argparse
import glob
import os
import sys

from rotifer.cli import input_files, key_value, write_rows
from rotifer.io import open_input

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="HMMER report(s) or directories of reports (default: standard input)")
    parser.add_argument("-y", "-a", "--header", "--add_header", dest="header", action="store_true",
                        help="add column names to the output")
    parser.add_argument("-e", "--evalue", "--pcut", dest="evalue", type=float,
                        help="maximum domain E-value")
    parser.add_argument("-c", "--clean", dest="clean", type=key_value, action="append", default=[],
                        metavar="TARGET=CLEANER",
                        help="clean identifiers: TARGET is model or sequence, CLEANER is version "
                             "(drop a trailing .<number>) or profiledb; e.g. -c model=version")
    parser.add_argument("-b", "--force_basename", dest="force_basename", action="append", default=[],
                        metavar="SUFFIX",
                        help="use the report's file name, without SUFFIX, as model name")
    parser.add_argument("-x", "--extension", dest="extension", default=".hmmsearch",
                        help="suffix of the reports to read from directories (default: %(default)s)")

def _expand(paths, extension):
    """
    Replace directories by the reports they contain.
    """
    files = []
    for path in paths:
        if os.path.isdir(path):
            files.extend(sorted(glob.glob(os.path.join(path, "*" + extension))))
        elif os.path.isfile(path) or path == "-":
            files.append(path)
        else:
            raise SystemExit(f"hmmer2table: {path}: not a file or directory")
    return files

def run(args):
    from rotifer.io.hmmer import COLUMNS, parse_hmmer
    clean = dict(args.clean)
    write_rows([], header=COLUMNS if args.header else None)
    for path in _expand(input_files(args.files), args.extension):
        with open_input(path) as handle:
            write_rows(parse_hmmer(handle, evalue=args.evalue, clean=clean,
                                   force_basename=args.force_basename, filename=path))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="hmmer2table", description="Convert HMMER 3 reports to a domain table.")
    add_arguments(parser)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except ValueError as error:
        parser.error(str(error))

if __name__ == "__main__":
    sys.exit(main())
