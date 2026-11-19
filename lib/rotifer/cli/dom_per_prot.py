__doc__ = """
dom_per_prot: number of regions in each protein.

Input: an architecture table, as written by domain2architecture.

Library equivalent: ``rotifer.tools.domains.dom_per_prot``.
"""

import argparse
import sys

from rotifer.cli import write_rows

def add_arguments(parser):
    parser.add_argument("architectures", metavar="FILE", help="architecture table (- for standard input)")

def run(args):
    from rotifer.tools.domains import dom_per_prot
    table = dom_per_prot(sys.stdin if args.architectures == "-" else args.architectures)
    write_rows(table.astype(str).values.tolist(), header=list(table.columns))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="dom_per_prot", description="Count regions in each protein.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
