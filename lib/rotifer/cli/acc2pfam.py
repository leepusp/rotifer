__doc__ = """
acc2pfam: Pfam domains, signal peptides and transmembrane segments of
proteins, as a padded table (ID, domain, start, end, evalue).

Runs hmmscan against Pfam and phobius. Unless a FASTA file is given,
sequences are fetched from the local sequence database or, failing that,
from NCBI (``rotifer.db.proteins``).

Library equivalent: ``rotifer.tools.domains.acc2pfam``.
"""

import argparse
import sys

from rotifer.cli import write_padded

def add_arguments(parser):
    parser.add_argument("accessions", nargs="*", metavar="ACCESSION", help="protein accession(s)")
    parser.add_argument("-f", "--fasta", dest="fasta", default=None,
                        help="read sequences from this FASTA file (- for standard input) instead of fetching them")
    parser.add_argument("-d", "--database", dest="database", default=None,
                        help="Pfam HMM library (default: $ROTIFER_DATA/pfam/Pfam-A.hmm)")
    parser.add_argument("--cpu", dest="cpu", type=int, default=4, help="hmmscan threads (default: %(default)s)")

def run(args):
    from rotifer.io import open_input
    from rotifer.tools.domains import pfam_search_rows
    fasta = None
    if args.fasta:
        with open_input(args.fasta) as handle:
            fasta = handle.read()
    elif not args.accessions:
        raise SystemExit("acc2pfam: give protein accessions or --fasta")
    write_padded(pfam_search_rows(args.accessions, fasta, args.database, args.cpu))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="acc2pfam", description="Pfam domains, signal peptides and "
                                     "transmembrane segments of proteins.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
