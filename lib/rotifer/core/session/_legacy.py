"""Reader for old ROTIFER shelve session files.

This module will be removed on 2026-11-30. It exists only to read and convert
sessions saved by the old shelve-based ``save_session``. Nothing in ROTIFER
writes shelve files any more.

The reader opens the file once, read-only, and never changes, locks for
writing, or deletes the old files. Convert an old session to the new format
with ``rotifer.core.session.Session.convert(path)``.
"""

import dbm
import glob
import shelve
import warnings
from pathlib import Path

from rotifer.core.logger import get_logger

logger = get_logger(__name__)

REMOVAL_DATE = "2026-11-30"


def is_legacy_session(path):
    """True if `path` is an existing dbm/shelve file. Never raises for ordinary paths."""
    try:
        p = Path(str(path))
        if p.is_dir():
            return False
        return bool(dbm.whichdb(str(p)))
    except (OSError, ValueError, TypeError):
        return False


def legacy_target(path):
    """Directory a legacy session at `path` is converted to: ``<path>.session``."""
    return Path(str(path) + ".session")


def read_legacy(path):
    """Read an old shelve session without touching it.

    Returns ``(values, errors)``: `values` maps variable names to unpickled
    objects, `errors` maps names that could not be unpickled to
    ``"ExcType: message"``.
    """
    if not is_legacy_session(path):
        p = str(path)
        if Path(p).is_dir():
            raise ValueError(f"{p} is a directory, not a shelve session file")
        if not Path(p).exists() and not glob.glob(glob.escape(p) + ".*"):
            raise FileNotFoundError(f"No session file found at {p}")
        raise ValueError(f"{p} is not a shelve session file")
    warnings.warn(
        f"Shelve session files are deprecated and this reader will be removed on "
        f"{REMOVAL_DATE}. Convert {path} with rotifer.core.session.Session.convert(path), "
        f"which writes {legacy_target(path)}/.",
        FutureWarning, stacklevel=2)
    values, errors = {}, {}
    with shelve.open(str(path), flag="r") as shelf:
        for key in list(shelf.keys()):
            try:
                values[key] = shelf[key]
            except Exception as exc:
                errors[key] = f"{type(exc).__name__}: {exc}"
                logger.warning("Could not load %r from %s: %s", key, path, errors[key])
    return values, errors
