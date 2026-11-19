__doc__ = """
distinct_dom_per_prot: number of distinct domains in each protein.

Input: an architecture table, as written by domain2architecture.

Library equivalent: ``rotifer.tools.domains.distinct_dom_per_prot``.
"""

import argparse
import sys

from rotifer.cli import write_rows

def add_arguments(parser):
    parser.add_argument("architectures", metavar="FILE", help="architecture table (- for standard input)")

def run(args):
    from rotifer.tools.domains import distinct_dom_per_prot
    table = distinct_dom_per_prot(sys.stdin if args.architectures == "-" else args.architectures)
    write_rows(table.astype(str).values.tolist(), header=list(table.columns))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="distinct_dom_per_prot", description="Count distinct domains in each protein.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
