__doc__ = """
architecture2table: expand an architecture table (domain2architecture's
output) into one row per region: ID, domain, start, end, evalue.

Library equivalent: ``rotifer.interval.architecture.architecture2table``.
"""

import argparse
import sys

from rotifer.cli import input_files, write_rows
from rotifer.io import open_input

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="architecture table(s) (default: standard input)")
    parser.add_argument("-d", "--domain_delimiter", dest="domain_delimiter", default=None,
                        help="regular expression separating domain names in architectures "
                             r"(default: [~+\[\]]+)")
    parser.add_argument("-s", "--input_delimiter", dest="delimiter", default="\t",
                        help="input column delimiter (default: TAB)")
    parser.add_argument("-r", "--output_delimiter", dest="output_delimiter", default="\t",
                        help="output column delimiter (default: TAB)")

def _lines(paths):
    """
    All input files as a single stream of lines (only the first header
    line is recognized as a header).
    """
    for path in paths:
        with open_input(path) as handle:
            yield from handle

def run(args):
    from rotifer.interval.architecture import DOMAIN_DELIMITER, iter_architecture_table
    write_rows(iter_architecture_table(_lines(input_files(args.files)), delimiter=args.delimiter,
                                       domain_delimiter=args.domain_delimiter or DOMAIN_DELIMITER),
               sep=args.output_delimiter)
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="architecture2table",
                                     description="Expand an architecture table into one row per region.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
