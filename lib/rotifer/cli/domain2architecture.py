__doc__ = """
domain2architecture: build domain architectures from tables of regions.

Input tables are TAB-separated, without header, with columns sequence,
domain, start, end and E-value (as written by hmmer2table, blast2table -s
and phobius2table); extra columns are ignored.

Library equivalent: ``rotifer.interval.architecture.domain2architecture``.
"""

import argparse
import os
import sys

from rotifer.cli import input_files, key_value, write_rows
from rotifer.io import open_input

def _words(values):
    """
    Expand list options: a value naming an existing file is replaced by
    the words in that file.
    """
    words = []
    for value in values:
        if os.path.isfile(value):
            with open(value) as handle:
                words.extend(handle.read().split())
        else:
            words.append(value)
    return words

def _mapping(pairs):
    """
    Build a dict from KEY=VALUE options; ``file=PATH`` loads a two-column
    table.
    """
    mapping = {}
    for key, value in pairs:
        if key == "file":
            with open(value) as handle:
                mapping.update(line.split()[:2] for line in handle if len(line.split()) >= 2)
        else:
            mapping[key] = value
    return mapping

def add_arguments(parser):
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="domain table(s) (default: standard input)")
    parser.add_argument("-e", "--evalue", "--pcut", dest="evalue", type=float, default=0.01,
                        help="maximum E-value of accepted regions (default: %(default)s)")
    parser.add_argument("-o", "--overlap", dest="overlap", default="0.4",
                        help="maximum overlap between conflicting regions: a fraction (e.g. 0.4) "
                             "or a number of residues (e.g. 10) (default: %(default)s)")
    parser.add_argument("-l", "--length", dest="by_length", action="store_true",
                        help="prefer longer regions instead of lower E-values")
    parser.add_argument("-i", "--ignore_overlap", dest="ignore_overlap", action="append", default=[],
                        metavar="DOMAIN", help="never discard DOMAIN for overlaps (all: any domain); "
                                               "a file name loads a list")
    parser.add_argument("-x", "--exclude", dest="exclude", action="append", default=[],
                        metavar="DOMAIN", help="ignore DOMAIN; a file name loads a list")
    parser.add_argument("-m", "--rename", dest="rename", type=key_value, action="append", default=[],
                        metavar="OLD=NEW", help="rename domains (file=PATH loads a two-column table)")
    parser.add_argument("-tc", "--tc", "--target_cutoff", dest="target_cutoff", type=key_value,
                        action="append", default=[], metavar="DOMAIN=EVALUE",
                        help="domain-specific E-value cutoff (file=PATH loads a two-column table)")
    parser.add_argument("-s", "--delimiter", "--input_delimiter", dest="delimiter", default="\t",
                        help="column delimiter of input tables (regular expression; default: TAB)")
    parser.add_argument("-r", "--output", "--output_delimiter", dest="output_delimiter", default="\t",
                        help="column delimiter of the output header (default: TAB)")
    parser.add_argument("-ia", "--ia", "--input_args", dest="input_args", type=key_value, action="append",
                        default=[], metavar="KEY=VALUE",
                        help="input options: columns=<comma-separated names>, header=1")
    parser.add_argument("-oa", "--oa", "--output_args", dest="output_args", type=key_value, action="append",
                        default=[], metavar="KEY=VALUE",
                        help="output options: header=0, connector=<string>, columns=<comma-separated names>")

def run(args):
    from rotifer.interval.architecture import INPUT_COLUMNS, iter_architectures, parse_regions
    input_args = {"columns": ",".join(INPUT_COLUMNS), **dict(args.input_args)}
    output_args = {"header": "1", "connector": "+", "columns": "ID,architecture,regions,evalue",
                   **dict(args.output_args)}
    data = {}
    for path in input_files(args.files):
        with open_input(path) as handle:
            parse_regions(handle, data, evalue=args.evalue, target_cutoff=_mapping(args.target_cutoff),
                          exclude=set(_words(args.exclude)), rename=_mapping(args.rename),
                          columns=[c.strip() for c in input_args["columns"].split(",")],
                          header=input_args.get("header", "0") not in ("", "0"),
                          delimiter=args.delimiter)
    if output_args["header"] not in ("", "0"):
        sys.stdout.write(args.output_delimiter.join(c.strip() for c in output_args["columns"].split(",")) + "\n")
    write_rows(iter_architectures(data, overlap=args.overlap, by_length=args.by_length,
                                  ignore_overlap=set(_words(args.ignore_overlap)),
                                  connector=output_args["connector"]))
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="domain2architecture",
                                     description="Build domain architectures from tables of sequence regions.")
    add_arguments(parser)
    return run(parser.parse_args(argv))

if __name__ == "__main__":
    sys.exit(main())
