"""Safe file writes and the session lock.

Writers replace files atomically, so a reader sees either the old or the new
content and never a partial file. Only writers take the lock; readers do not.
"""

import contextlib
import os
import tempfile
from pathlib import Path

import filelock

LOCK_NAME = ".lock"


class SessionLocked(RuntimeError):
    """
    Another process holds the session lock.

    See Also
    --------
    session_lock : Hold the exclusive lock of a session.
    """


@contextlib.contextmanager
def atomic_writer(path, mode="wb"):
    """
    Open a temporary file next to `path` and replace `path` with it on success.

    The data is flushed and fsynced before `os.replace`. If the block raises,
    the temporary file is removed and `path` keeps its old content.

    Parameters
    ----------
    path : str or pathlib.Path
        File to write. Its parent directory is created if missing.
    mode : str, default "wb"
        Mode the temporary file is opened with.

    Yields
    ------
    file object
        The open temporary file, to write to inside the `with` block.

    See Also
    --------
    atomic_write_bytes : Write bytes to a file atomically.
    atomic_write_text : Write text to a file atomically.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, mode) as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def atomic_write_bytes(path, data):
    """
    Write `data` (bytes) to `path` atomically.

    Parameters
    ----------
    path : str or pathlib.Path
        File to write.
    data : bytes
        Content of the file.

    See Also
    --------
    atomic_write_text : Write text to a file atomically.
    atomic_writer : Context manager for writing a file atomically.
    """
    with atomic_writer(path, "wb") as handle:
        handle.write(data)


def atomic_write_text(path, text, encoding="utf-8"):
    """
    Write `text` to `path` atomically.

    Parameters
    ----------
    path : str or pathlib.Path
        File to write.
    text : str
        Content of the file.
    encoding : str, default "utf-8"
        Encoding used to turn `text` into bytes.

    See Also
    --------
    atomic_write_bytes : Write bytes to a file atomically.
    atomic_writer : Context manager for writing a file atomically.
    """
    atomic_write_bytes(path, text.encode(encoding))


@contextlib.contextmanager
def session_lock(directory, timeout=10.0):
    """
    Hold the exclusive lock of the session in `directory`.

    The lock file is `<directory>/.lock`.

    Parameters
    ----------
    directory : str or pathlib.Path
        The session directory. It is created if it does not exist.
    timeout : float, default 10.0
        Seconds to wait for the lock.

    Yields
    ------
    None
        Nothing; the lock is held inside the `with` block and released on
        exit.

    Raises
    ------
    SessionLocked
        If the lock is not obtained within `timeout` seconds.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    lock = filelock.FileLock(str(directory / LOCK_NAME))
    try:
        lock.acquire(timeout=timeout)
    except filelock.Timeout as exc:
        raise SessionLocked(
            f"session {directory} is locked by another process "
            f"(waited {timeout} s)"
        ) from exc
    try:
        yield
    finally:
        lock.release()
