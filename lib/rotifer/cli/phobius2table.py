__doc__ = """
phobius2table: convert Phobius predictions to a domain table.

Library equivalent: ``rotifer.io.phobius.read_phobius``.
"""

import argparse
import sys

from rotifer.cli import input_files, write_rows
from rotifer.io import open_input

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="Phobius report(s) in the default format (default: standard input)")
    parser.add_argument("-e", "--evalue", dest="evalue", default="0",
                        help="value of the E-value column (default: %(default)s)")

def run(args):
    from rotifer.io.phobius import parse_phobius
    for path in input_files(args.files):
        with open_input(path) as handle:
            write_rows(parse_phobius(handle, evalue=args.evalue))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="phobius2table", description="Convert Phobius predictions to a domain table.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
