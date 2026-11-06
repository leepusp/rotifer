__doc__ = """
blast2table: convert BLAST+ pairwise text reports to a table of HSPs.

Library equivalent: ``rotifer.io.blast.read_blast``.
"""

import argparse
import sys

from rotifer.cli import input_files, write_rows
from rotifer.io import open_input

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="BLAST+ text report(s) (default: standard input)")
    parser.add_argument("-s", "--swap", dest="swap", action="store_true",
                        help="swap query and hit columns (e.g. RPS-BLAST to domain table)")
    parser.add_argument("-c", "--clean", dest="clean", action="append", default=[], metavar="CLEANER",
                        help="clean identifiers: cdd (domain names from RPS-BLAST descriptions), "
                             "genbank2gi, genbank2acc, profiledb or none; may be repeated")
    parser.add_argument("-e", "--evalue", "--pcut", dest="evalue", type=float, default=10,
                        help="maximum HSP E-value (default: %(default)s)")
    parser.add_argument("-t", "--threshold", dest="threshold", type=float, default=0.01,
                        help="PSI-BLAST inclusion threshold (default: %(default)s)")
    parser.add_argument("-o", "--overlap", dest="overlap", default="0.4",
                        help="overlap above which HSPs are redundant: a fraction or a number of "
                             "residues (default: %(default)s)")
    parser.add_argument("--nonr", "--no_non_redundant_hsps", dest="nonredundant", action="store_false",
                        help="keep redundant (overlapping) HSPs")
    parser.add_argument("-k", "--column", dest="columns", action="append", default=[], metavar="COLUMN",
                        help="output column; may be repeated (default: the standard domain table)")
    parser.add_argument("-y", "--header", dest="header", action="store_true",
                        help="add column names to the output")
    parser.add_argument("-r", "--rename_query", dest="rename_query", default=None,
                        help="new query name")
    parser.add_argument("-b", "--force_basename", dest="force_basename", action="append", default=[],
                        metavar="SUFFIX", help="use the report's file name, without SUFFIX, as query name")
    parser.add_argument("--out", "--outfile", "--output", dest="output", default=None,
                        help="output file (default: standard output)")

def run(args):
    from rotifer.io.blast import AVAILABLE_COLUMNS, COLUMNS, iter_blast
    columns = args.columns or COLUMNS
    unknown = [column for column in columns if column not in AVAILABLE_COLUMNS]
    if unknown:
        raise SystemExit(f"blast2table: unknown column(s) {', '.join(unknown)}; "
                         f"available: {', '.join(AVAILABLE_COLUMNS)}")
    out = open(args.output, "w") if args.output else sys.stdout
    try:
        write_rows([], header=columns if args.header else None, out=out)
        for path in input_files(args.files):
            with open_input(path) as handle:
                write_rows(iter_blast(handle, columns=columns, swap=args.swap, evalue=args.evalue,
                                      threshold=args.threshold, overlap=args.overlap,
                                      nonredundant=args.nonredundant, clean=args.clean,
                                      rename_query=args.rename_query, force_basename=args.force_basename,
                                      filename=path), out=out)
    finally:
        if args.output:
            out.close()
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="blast2table", description="Convert BLAST+ text reports to a table of HSPs.")
    add_arguments(parser)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except ValueError as error:
        parser.error(str(error))

if __name__ == "__main__":
    sys.exit(main())
