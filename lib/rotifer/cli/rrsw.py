__doc__ = """
rrsw: download and unpack mirrors of public databases.

Library equivalent: ``rotifer.tools.rrsw.rrsw`` (see that module for the
configuration files).
"""

import argparse
import json
import sys

def add_arguments(parser):
    parser.add_argument("targets", nargs="*", metavar="TARGET",
                        help="target name(s) from the table of targets, or all")
    parser.add_argument("-l", "--list", dest="list", action="store_true",
                        help="show the targets in the table and exit")
    parser.add_argument("-n", "--test", dest="test", action="store_true",
                        help="show what would be done, without changing anything")
    parser.add_argument("-t", "--table", dest="table", default=None,
                        help="table of targets (default: rrsw.<host>.yml or rrsw.yml in etc/rotifer/rrsw)")
    parser.add_argument("-s", "--source", dest="source", default=None,
                        help="root of downloaded files (default: <data>/src)")
    parser.add_argument("-p", "--processed", dest="processed", default=None,
                        help="root of processed files (default: rotifer's data directory)")
    parser.add_argument("-x", "--exclude", dest="exclude", action="append", default=[], metavar="TARGET",
                        help="skip this target; may be repeated")
    parser.add_argument("--nodownload", "--no-download", dest="download", action="store_false",
                        help="do not download, only process")
    parser.add_argument("--noprocess", "--no-process", dest="process", action="store_false",
                        help="do not process downloaded files")

def _format_table(rows):
    """
    Render the table of targets as a text grid.
    """
    cells = [["" if value is None else json.dumps(value, separators=(",", ":")) if isinstance(value, (dict, list))
              else str(value) for value in row] for row in rows]
    widths = [max(len(row[i]) for row in cells) for i in range(len(cells[0]))] if cells else []
    rule = "+" + "+".join("-" * (width + 2) for width in widths) + "+"
    lines = [rule]
    for i, row in enumerate(cells):
        lines.append("| " + " | ".join(cell.ljust(width) for cell, width in zip(row, widths)) + " |")
        if i == 0:
            lines.append(rule)
    lines.append(rule)
    return "\n".join(lines) + "\n"

def run(args):
    import logging
    from rotifer.tools.rrsw import load_table, rrsw, select_targets
    # Report each download (and, with -n, what would be downloaded), as before
    logging.getLogger("rotifer.tools.rrsw").setLevel(logging.INFO)
    if args.list:
        sys.stdout.write(_format_table(select_targets(load_table(args.table), args.targets, args.exclude)))
        return 0
    rrsw(args.targets, table=args.table, source=args.source, processed=args.processed,
         download=args.download, process=args.process, exclude=args.exclude, test=args.test)
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(prog="rrsw", description="Download and unpack mirrors of public databases.")
    add_arguments(parser)
    args = parser.parse_args(argv)
    if not args.list and not args.targets:
        parser.print_help()
        return 1
    return run(args)

if __name__ == "__main__":
    sys.exit(main())
