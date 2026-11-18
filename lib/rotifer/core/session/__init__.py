"""
rotifer.core.session: save and restore the variables of a working session.

A session is a directory with one or more files per variable and a JSON manifest
(`manifest.json`) that records what each file holds, how it was written and
the versions of Python and the main packages. The manifest structure is given
by `manifest_schema()`; `python -m rotifer.core.session schema` prints it.

This package does not print to stdout and does not ask questions. Messages go
to the `rotifer` logger.
"""

from rotifer.core.session._context import default_namespace
from rotifer.core.session._describe import describe
from rotifer.core.session._manifest import manifest_schema

from rotifer.core.session._session import DiffResult, Report, Session, open_session

__all__ = ["DiffResult", "Report", "Session", "default_namespace", "describe",
           "manifest_schema", "open_session"]
