"""
Print the session manifest schema from the command line.

Run as `python -m rotifer.core.session schema`. Any other argument prints the
usage line to stderr.
"""

import json
import sys

from rotifer.core.session._manifest import manifest_schema

USAGE = "usage: python -m rotifer.core.session schema\n"


def main(argv=None):
    """
    Run the command line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments, without the program name. Default:
        `sys.argv[1:]`.

    Returns
    -------
    int
        Exit status: 0 if the schema was written, 2 otherwise.

    See Also
    --------
    manifest_schema : The schema this prints.

    Notes
    -----
    `schema` writes the JSON Schema of manifest v1 to stdout. Any other
    argument writes the usage line to stderr.

    Examples
    --------
    >>> import contextlib
    >>> import io
    >>> from rotifer.core.session.__main__ import main
    >>> buffer = io.StringIO()
    >>> with contextlib.redirect_stdout(buffer):
    ...     status = main(["schema"])
    >>> status
    0
    >>> buffer.getvalue().lstrip().startswith("{")
    True
    """
    argv = sys.argv[1:] if argv is None else argv
    if argv == ["schema"]:
        sys.stdout.write(json.dumps(manifest_schema(), indent=2) + "\n")
        return 0
    sys.stderr.write(USAGE)
    return 2


if __name__ == "__main__":
    sys.exit(main())
