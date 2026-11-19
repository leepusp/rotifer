__doc__ = """
dat2fa: extract sequences from DaliLite .dat files as FASTA.

Library equivalent: ``rotifer.io.dalidat.parse_dat``.
"""

import argparse
import sys

from rotifer.cli import input_files
from rotifer.io import open_input

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="DaliLite .dat file(s) (default: standard input)")

def run(args):
    from rotifer.io.dalidat import dat2fasta
    for path in input_files(args.files):
        with open_input(path) as handle:
            sys.stdout.writelines(dat2fasta(handle))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="dat2fa", description="Extract sequences from DaliLite .dat files as FASTA.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
