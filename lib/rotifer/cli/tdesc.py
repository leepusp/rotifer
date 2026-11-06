__doc__ = """
tdesc: number the columns of a text table (index and name of each
column in the first line).

Library equivalent: ``rotifer.tools.tables.tdesc``.
"""

import argparse
import sys

from rotifer.cli import write_rows

def add_arguments(parser):
    parser.add_argument("table", metavar="FILE", help="table (- for standard input)")
    parser.add_argument("sep", nargs="?", default="\t", metavar="SEP",
                        help="column delimiter, a regular expression (default: TAB)")

def run(args):
    from rotifer.tools.tables import tdesc
    write_rows([str(index), name] for index, name in tdesc(args.table, args.sep or "\t"))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="tdesc", description="Number the columns of a text table.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
