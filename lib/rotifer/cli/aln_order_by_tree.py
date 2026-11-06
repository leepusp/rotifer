__doc__ = """
aln_order_by_tree: reorder the sequences of an alignment to follow the
leaves of a tree, writing FASTA (one line per sequence).

Library equivalent: ``rotifer.tools.alignments.aln_order_by_tree``.
"""

import argparse
import sys

def add_arguments(parser):
    parser.add_argument("alignment", metavar="ALIGNMENT", help="alignment file")
    parser.add_argument("tree", metavar="TREE", help="tree file")
    parser.add_argument("--alignment_format", default="fasta", help="alignment format (default: %(default)s)")
    parser.add_argument("--tree_format", default="newick", help="tree format (default: %(default)s)")

def run(args):
    from rotifer.tools.alignments import aln_order_by_tree
    for identifier, sequence in aln_order_by_tree(args.alignment, args.tree, args.alignment_format, args.tree_format):
        sys.stdout.write(f">{identifier}\n{sequence}\n")
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="aln_order_by_tree",
                                     description="Reorder an alignment to follow the leaves of a tree.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
