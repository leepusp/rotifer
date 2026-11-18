"""Command line: python -m rotifer.core.session schema"""

import json
import sys

from rotifer.core.session._manifest import manifest_schema

USAGE = "usage: python -m rotifer.core.session schema\n"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv == ["schema"]:
        sys.stdout.write(json.dumps(manifest_schema(), indent=2) + "\n")
        return 0
    sys.stderr.write(USAGE)
    return 2


if __name__ == "__main__":
    sys.exit(main())
