"""
Save and load the variables of a working session.

A session is a directory with a `manifest.json` and an `objects/` folder that
holds one set of files per variable. `Session` is the only entry point. It
never prints and never asks questions: results come back as reports, and
messages go to the logger.

Only writers take the lock (`save`, `remove`). Readers (`load`, `get`, `info`,
`ls`, `diff`) read the manifest, which is always replaced atomically.
"""
from __future__ import annotations

import fnmatch
from glob import escape as glob_escape
import hashlib
import json
import os
import re
import secrets
import shutil
import sys
import types
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from rotifer.core.logger import get_logger
from rotifer.core.session import _context, codecs
from rotifer.core.session._describe import describe
from rotifer.core.session._files import session_lock
from rotifer.core.session._fingerprint import fingerprint, hash_order_dependent
from rotifer.core.session._manifest import (MANIFEST_NAME, Manifest, VariableEntry,
                                       read_manifest, write_manifest)

logger = get_logger(__name__)

OBJECTS_DIR = "objects"
DIFF_LEVELS = ("names", "summary", "content")
IPYTHON_NAMES = frozenset({"In", "Out", "exit", "quit", "get_ipython", "open"})
_SUMMARY_KEYS = ("kind", "length", "shape", "columns", "dtypes")


def _legacy_module():
    """The legacy reader, or None if it has been removed."""
    try:
        from rotifer.core.session import _legacy
    except ImportError:
        return None
    return _legacy


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Report:
    """
    Record what a `Session` call did.

    `Session.save`, `Session.load` and `Session.remove` return one, and a
    conversion from an old shelve session stores one in
    `Session.conversion_report`. `to_dict` gives a plain dict that can be
    written as JSON.

    Parameters
    ----------
    action : str
        What the call did: "save", "load", "remove" or "convert".
    path : str
        Path of the session directory.
    saved : list of str, optional
        Names of the variables this call saved.
    loaded : list of str, optional
        Names of the variables this call loaded.
    skipped : dict of str to str, optional
        Names this call left out, mapped to the reason they were left out.
    failed : dict of str to str, optional
        Names that could not be saved or loaded, mapped to
        "ExcType: message".
    removed : list of str, optional
        Names of the variables this call removed.
    warnings : list of str, optional
        Messages produced by the call, such as a saved environment that
        differs from the running one.
    converted_from : str or None, optional
        Path of the old shelve session this session was converted from, or
        None if it was not converted.

    See Also
    --------
    DiffResult : What `Session.diff` returns instead.

    Examples
    --------
    >>> import os, tempfile
    >>> from rotifer.core.session import Session
    >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
    >>> report = session.save(ns={"x": 1, "y": 2})
    >>> report.action, sorted(report.saved), report.failed
    ('save', ['x', 'y'], {})
    """

    action: str
    path: str
    saved: list = field(default_factory=list)
    loaded: list = field(default_factory=list)
    skipped: dict = field(default_factory=dict)
    failed: dict = field(default_factory=dict)
    removed: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    converted_from: str | None = None

    def to_dict(self):
        """
        Return the report as a JSON-serializable dict.

        Returns
        -------
        dict
            One key per field: "action", "path", "saved", "loaded",
            "skipped", "failed", "removed", "warnings" and
            "converted_from".

        Examples
        --------
        >>> from rotifer.core.session import Report
        >>> report = Report(action="save", path="s", saved=["x"])
        >>> report.to_dict()["saved"]
        ['x']
        >>> report.to_dict()["converted_from"] is None
        True
        """
        return {
            "action": self.action,
            "path": self.path,
            "saved": list(self.saved),
            "loaded": list(self.loaded),
            "skipped": dict(self.skipped),
            "failed": dict(self.failed),
            "removed": list(self.removed),
            "warnings": list(self.warnings),
            "converted_from": self.converted_from,
        }


@dataclass
class DiffResult:
    """
    Record the result of `Session.diff`.

    Parameters
    ----------
    level : str
        Comparison level: "names", "summary" or "content".
    path : str
        Path of the session directory.
    added : list of str, optional
        Names in the namespace but not in the session.
    removed : list of str, optional
        Names in the session but not in the namespace.
    changed : list of str, optional
        Names in both, different at the chosen level.
    unchanged : list of str, optional
        Names in both with no difference found at the chosen level. At level
        "names" it means only "present in both".
    unknown : list of str, optional
        Names in both that could not be compared: the live object is larger
        than `max_bytes`, or it cannot be fingerprinted.

    See Also
    --------
    Session.diff : The method that returns it.

    Examples
    --------
    >>> import os, tempfile
    >>> from rotifer.core.session import Session
    >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
    >>> _ = session.save(ns={"x": 1, "y": 2})
    >>> result = session.diff(ns={"x": 1, "z": 3})
    >>> result.added, result.removed, result.unchanged
    (['z'], ['y'], ['x'])
    """

    level: str
    path: str
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    changed: list = field(default_factory=list)
    unchanged: list = field(default_factory=list)
    unknown: list = field(default_factory=list)

    def to_dict(self):
        """
        Return the result as a JSON-serializable dict.

        Returns
        -------
        dict
            One key per field: "level", "path", "added", "removed",
            "changed", "unchanged" and "unknown".

        Examples
        --------
        >>> from rotifer.core.session import DiffResult
        >>> result = DiffResult(level="names", path="s", added=["y"])
        >>> result.to_dict()["added"]
        ['y']
        >>> result.to_dict()["level"]
        'names'
        """
        return {
            "level": self.level,
            "path": self.path,
            "added": list(self.added),
            "removed": list(self.removed),
            "changed": list(self.changed),
            "unchanged": list(self.unchanged),
            "unknown": list(self.unknown),
        }


def _matches(name, patterns):
    return any(fnmatch.fnmatchcase(name, p) for p in patterns)


def _as_list(value):
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    return list(value)


def _is_glob(pattern):
    return any(c in pattern for c in "*?[")


def _skip_reason(name, value):
    """Why the default filters leave this variable out, or None to keep it."""
    if not isinstance(name, str):
        return "name is not a string"
    if name.startswith("_"):
        return "name starts with underscore"
    if name in IPYTHON_NAMES:
        return "IPython name"
    if isinstance(value, types.ModuleType):
        return "module"
    func = value
    if isinstance(value, types.MethodType):
        func = value.__func__
    if isinstance(func, (types.FunctionType, types.BuiltinFunctionType, type)):
        if getattr(func, "__module__", None) != "__main__":
            return "imported function or class (comes back by import)"
    return None


def _file_prefix(name):
    digest = hashlib.blake2b(name.encode("utf-8"), digest_size=4).hexdigest()
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)[:80]
    return f"{safe}-{digest}"


def _type_name(value):
    cls = type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def _error_text(exc):
    return f"{type(exc).__name__}: {exc}"


def _estimated_size(value):
    try:
        import pandas as pd
        if isinstance(value, (pd.DataFrame, pd.Series)):
            usage = value.memory_usage(deep=False)
            return int(usage.sum() if hasattr(usage, "sum") else usage)
    except ImportError:
        pass
    try:
        import numpy as np
        if isinstance(value, np.ndarray):
            return int(value.nbytes)
    except ImportError:
        pass
    return sys.getsizeof(value)


def _safe_describe(value):
    try:
        return describe(value)
    except Exception:
        return {"kind": "object"}


def _comparable_summary(summary):
    kept = {k: summary.get(k) for k in _SUMMARY_KEYS if summary.get(k) is not None}
    return json.loads(json.dumps(kept))


class Session:
    """
    Open a session directory.

    A session holds the variables of a working session as one set of files
    per variable, described by a manifest. The session never prints and never
    asks questions: results come back as reports.

    Parameters
    ----------
    path : str or pathlib.Path
        Path of the session directory. If it is an old shelve session, it is
        converted to `<path>.session/` and that directory is used.
    create : bool, default True
        Create the session at `path` if it does not exist. If False, opening
        a path without a session raises `FileNotFoundError`.
    lock_timeout : float, default 10.0
        Seconds a writer waits for the session lock before giving up.

    Attributes
    ----------
    path : pathlib.Path
        Directory of the session. For a converted legacy session this is
        `<path>.session`.
    lock_timeout : float
        Seconds a writer waits for the session lock before giving up.
    converted_from : str or None
        Path of the old shelve session this session was converted from, or
        None if it was not converted.
    conversion_report : Report or None
        Report of the conversion (`action` is "convert"), or None if this
        session was not converted.

    Raises
    ------
    FileNotFoundError
        If `create` is False and there is no session at `path`.
    ValueError
        If `path` exists and is not a session directory.

    See Also
    --------
    open_session : Open an existing session without creating it.

    Notes
    -----
    Only writers (`save`, `remove`, and the conversion done by the
    constructor) take the session lock; readers (`load`, `get`, `info`, `ls`,
    `diff`) only read the manifest, which is always replaced atomically. A
    writer raises `rotifer.core.session._files.SessionLocked` if the lock is
    not free within `lock_timeout` seconds.

    When `path` is an old shelve session, the old files are never changed:
    the conversion is written to a temporary directory and renamed to
    `<path>.session` when it is complete. `converted_from` then holds the old
    path and `conversion_report` records what was converted.

    Examples
    --------
    >>> import os
    >>> import tempfile
    >>> from rotifer.core.session import Session
    >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
    >>> report = session.save({"x": 1})
    >>> report.saved
    ['x']
    >>> session.get("x")
    1
    """

    def __init__(self, path, create=True, lock_timeout=10.0):
        self.lock_timeout = lock_timeout
        self.converted_from = None
        self.conversion_report = None
        path = Path(path)
        _legacy = _legacy_module()
        if _legacy is not None and _legacy.is_legacy_session(path):
            target = _legacy.legacy_target(path)
            if self._is_complete_conversion(target):
                logger.warning("%s is an old shelve session; it was already "
                               "converted to %s, opening that", path, target)
                self.path = target
            else:
                self._convert_atomically(path, target)
            return
        self.path = path
        if (path / MANIFEST_NAME).exists():
            return
        if not create:
            raise FileNotFoundError(f"no session at {path}")
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"{path} exists and is not a session directory")
        self._init_dir()

    @classmethod
    def convert(cls, old, new=None, lock_timeout=10.0):
        """
        Convert an old shelve session into a session directory.

        The old files are read, never changed.

        Parameters
        ----------
        old : str or pathlib.Path
            Path of the old shelve session.
        new : str or pathlib.Path, optional
            Directory to write. Default `<old>.session`, the target `Session`
            uses when it converts a legacy session.
        lock_timeout : float, default 10.0
            Seconds a writer waits for the session lock before giving up.

        Returns
        -------
        Session
            The new session. Its `conversion_report` records what was
            converted.

        Raises
        ------
        ValueError
            If `new` is already a session or exists and is not a session
            directory, or if `old` is not a shelve session.
        FileNotFoundError
            If there is no shelve file at `old`.

        Warns
        -----
        FutureWarning
            If the old shelve files are read. The legacy reader is removed on
            2026-11-30.

        See Also
        --------
        Session : Opening a legacy session converts it automatically.

        Examples
        --------
        >>> import os
        >>> import shelve
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> old = os.path.join(tempfile.mkdtemp(), "old")
        >>> with shelve.open(old) as db:
        ...     db["n"] = 5
        >>> session = Session.convert(old)
        >>> session.path.name
        'old.session'
        >>> session.get("n")
        5
        """
        from rotifer.core.session import _legacy
        old = Path(old)
        new = Path(new) if new is not None else _legacy.legacy_target(old)
        self = cls.__new__(cls)
        self.lock_timeout = lock_timeout
        self.converted_from = None
        self.conversion_report = None
        self.path = new
        if (new / MANIFEST_NAME).exists():
            raise ValueError(f"{new} is already a session")
        self._convert_atomically(old, new)
        return self

    # ------------------------------------------------------------ internals

    @staticmethod
    def _is_complete_conversion(target):
        """True if `target` holds a finished conversion.

        A manifest at revision 0 with no variables and no `converted_from` is
        what an interrupted conversion by older code left behind: it is
        removed so the conversion can be redone.
        """
        if not (target / MANIFEST_NAME).exists():
            return False
        try:
            manifest = read_manifest(target)
        except Exception:
            return True   # unreadable: let the caller see the real error
        if (manifest.revision == 0 and not manifest.variables
                and "converted_from" not in manifest.context):
            shutil.rmtree(target, ignore_errors=True)
            return False
        return True

    def _convert_atomically(self, old, target):
        """Convert `old` into a temporary sibling of `target`, then rename it.

        `target` only ever appears complete. On failure the temporary
        directory is removed.
        """
        target = Path(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and (not target.is_dir() or any(target.iterdir())):
            raise ValueError(f"{target} exists and is not a session directory")
        tmp = target.parent / f"{target.name}.tmp-{secrets.token_hex(4)}"
        self.path = tmp
        try:
            self._init_dir()
            self._convert_from(old)
            if target.is_dir():
                target.rmdir()      # empty (checked above)
            os.replace(tmp, target)
        except BaseException:
            shutil.rmtree(tmp, ignore_errors=True)
            self.converted_from = None
            self.conversion_report = None
            raise
        self.path = target
        if self.conversion_report is not None:
            self.conversion_report.path = str(target)
        logger.warning("converted old shelve session %s to the new session "
                       "directory %s; use that directory from now on", old, target)

    def _init_dir(self):
        with session_lock(self.path, self.lock_timeout):
            if not (self.path / MANIFEST_NAME).exists():
                (self.path / OBJECTS_DIR).mkdir(exist_ok=True)
                manifest = Manifest(created=_now(), updated=_now(),
                                    environment=_context.environment(),
                                    context=_context.notebook_context())
                write_manifest(self.path, manifest)

    def _convert_from(self, old):
        from rotifer.core.session import _legacy
        values, errors = _legacy.read_legacy(old)
        report = self._save(values, None, None, False, None, action="convert",
                            extra_failed=errors, use_filters=False)
        report.converted_from = str(old)
        self.converted_from = str(old)
        self.conversion_report = report

    def _manifest(self):
        return read_manifest(self.path)

    def _objects(self):
        return self.path / OBJECTS_DIR

    def _report(self, action):
        return Report(action=action, path=str(self.path),
                      converted_from=self.converted_from)

    def _delete_files(self, files):
        for rel in files:
            try:
                (self.path / rel).unlink(missing_ok=True)
            except OSError as exc:
                logger.warning("could not delete %s: %s", rel, exc)

    def _write_variable(self, name, value, compress, revision):
        """Write one variable under a name unique to this save.

        Returns (entry, new_files). If anything is raised, the files this
        call created are removed (never files the manifest lists).
        """
        prefix = self._objects() / f"{_file_prefix(name)}.r{revision}"
        try:
            codec, meta = codecs.encode(value, prefix, compress=compress)
            files = [f"{OBJECTS_DIR}/{f}" for f in meta["files"]]
            digest = hashlib.sha256()
            size = 0
            for rel in files:
                data = (self.path / rel).read_bytes()
                size += len(data)
                digest.update(data)
            quick = meta.pop("fingerprint", None) or fingerprint(value, codec)
            entry = VariableEntry(
                type=_type_name(value), codec=codec.name,
                codec_version=codec.version, files=files, bytes=size,
                sha256=digest.hexdigest(), fingerprint=quick, saved=_now(),
                version_bound=bool(meta.get("version_bound", False)),
                summary=_safe_describe(value), meta=meta)
        except BaseException:
            self._delete_unreferenced(
                [f"{OBJECTS_DIR}/{p.name}"
                 for p in self._objects().glob(glob_escape(prefix.name) + ".*")])
            raise
        return entry, files

    def _referenced(self):
        """Files the manifest on disk lists."""
        try:
            manifest = self._manifest()
        except Exception:
            return None
        return {f for e in manifest.variables.values() for f in e.files}

    def _delete_unreferenced(self, files):
        """Delete `files` except those the manifest on disk lists."""
        referenced = self._referenced()
        if referenced is None:      # cannot tell: delete nothing
            return
        self._delete_files([f for f in files if f not in referenced])

    def _sweep(self, manifest):
        """Delete files in objects/ that the manifest does not list."""
        keep = {Path(f).name for e in manifest.variables.values() for f in e.files}
        try:
            entries = list(self._objects().iterdir())
        except OSError:
            return
        for p in entries:
            if p.name not in keep and p.is_file():
                try:
                    p.unlink()
                except OSError as exc:
                    logger.warning("could not delete orphan %s: %s", p.name, exc)

    # ------------------------------------------------------------------ API

    def save(self, ns=None, names=None, exclude=None, prune=False, compress=None):
        """
        Save variables from a namespace into the session.

        Parameters
        ----------
        ns : dict-like, optional
            Namespace mapping names to values. Default: the user namespace,
            IPython's `user_ns` inside IPython and `__main__.__dict__`
            otherwise.
        names : str or sequence of str, optional
            Names or glob patterns of the variables to save. Default: all
            names, subject to the default filters.
        exclude : str or sequence of str, optional
            Names or glob patterns to leave out. Applied after `names`.
        prune : bool, default False
            Remove every saved variable inside the selection that this call
            does not save: gone from `ns`, or excluded by `exclude` or by the
            default filters. The selection is `names`, or every saved
            variable if `names` is None. Default: only add and update.
        compress : str or None, optional
            Compression for the pickle codecs: None (default) writes plain
            files, "gzip" writes gzipped files. Other codecs ignore it. Any
            other value puts the variables saved with pickle in `failed`.

        Returns
        -------
        Report
            The report: `saved` lists the names written, `skipped` gives the
            reason for each name in `names` that was left out (names dropped
            by the default filters or `exclude` when `names` is None are not
            listed), `failed` gives the error for each name that could not
            be written, and `removed` lists the names removed by `prune`.

        Raises
        ------
        SessionLocked
            If another process holds the session lock for longer than
            `lock_timeout`.

        See Also
        --------
        load : Put saved variables back into a namespace.

        Notes
        -----
        New files are written under names unique to this save and the
        manifest is replaced atomically; only then are the files it no longer
        lists deleted, so an interrupted save leaves the old session intact.

        The default filters skip names that start with "_", the IPython names
        (`In`, `Out`, `exit`, `quit`, `get_ipython`, `open`), modules, and
        functions or classes not defined in `__main__`.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1, "y": 2})
        >>> sorted(report.saved)
        ['x', 'y']
        >>> target = {}
        >>> report = session.load(names="x", ns=target)
        >>> target
        {'x': 1}
        """
        return self._save(ns, names, exclude, prune, compress)

    def _save(self, ns, names, exclude, prune, compress, action="save",
              extra_failed=None, use_filters=True):
        if ns is None:
            ns = _context.default_namespace()
        names = _as_list(names)
        exclude = _as_list(exclude) or []
        report = self._report(action)

        def selected(name):
            if names is not None and not _matches(name, names):
                return False
            return not _matches(name, exclude)

        todo = {}
        for name, value in ns.items():
            if not isinstance(name, str):
                continue
            if not selected(name):
                continue
            reason = _skip_reason(name, value) if use_filters else None
            if reason:
                if names is not None:
                    report.skipped[name] = reason
                continue
            todo[name] = value
        for pattern in names or []:
            if not _is_glob(pattern) and pattern not in ns:
                report.skipped[pattern] = "not in the namespace"

        with session_lock(self.path, self.lock_timeout):
            manifest = self._manifest()
            revision = manifest.revision + 1
            new_files = []        # created by this save
            stale = []            # no longer needed once the manifest is written
            try:
                for name, value in todo.items():
                    old = manifest.variables.get(name)
                    try:
                        entry, files = self._write_variable(name, value, compress,
                                                            revision)
                    except Exception as exc:
                        if old is not None:
                            stale.extend(old.files)
                            del manifest.variables[name]
                        manifest.failed[name] = _error_text(exc)
                        report.failed[name] = _error_text(exc)
                        logger.warning("could not save %s: %s", name, _error_text(exc))
                        continue
                    new_files.extend(files)
                    if old is not None:
                        stale.extend(old.files)
                    manifest.variables[name] = entry
                    manifest.failed.pop(name, None)
                    report.saved.append(name)
                for name, error in (extra_failed or {}).items():
                    manifest.failed[name] = str(error)
                    report.failed[name] = str(error)
                if prune:
                    for name in list(manifest.variables):
                        if name in todo or not (names is None or _matches(name, names)):
                            continue
                        stale.extend(manifest.variables.pop(name).files)
                        report.removed.append(name)
                    for name in list(manifest.failed):
                        if name not in todo and name not in (extra_failed or {}) \
                                and (names is None or _matches(name, names)):
                            del manifest.failed[name]
                manifest.revision = revision
                manifest.updated = _now()
                if not manifest.created:
                    manifest.created = manifest.updated
                manifest.environment = _context.environment()
                manifest.context = _context.notebook_context()
                if self.converted_from:
                    manifest.context["converted_from"] = self.converted_from
                write_manifest(self.path, manifest)
            except BaseException:
                self._delete_unreferenced(new_files)
                raise
            referenced = {f for e in manifest.variables.values() for f in e.files}
            self._delete_files([f for f in stale if f not in referenced])
            self._sweep(manifest)
        return report

    def load(self, names=None, ns=None):
        """
        Load saved variables into a namespace.

        A variable that cannot be read goes to `failed`; the other variables
        are still loaded.

        Parameters
        ----------
        names : str or sequence of str, optional
            Name or glob patterns of the variables to load. Default: all
            saved variables.
        ns : dict-like, optional
            Namespace that receives the values. Default: the user namespace,
            IPython's `user_ns` inside IPython and `__main__.__dict__`
            otherwise.

        Returns
        -------
        Report
            The report: `loaded` lists the names loaded, `skipped` gives the
            reason for each requested name that is not in the session,
            `failed` gives the error for each name that could not be read,
            and `warnings` holds any environment differences.

        See Also
        --------
        save : Write the variables of a namespace into the session.
        get : Read one variable without touching a namespace.

        Notes
        -----
        If the environment recorded in the manifest differs from the running
        one (Python, ROTIFER or the versions of the main packages), the
        difference is reported in the log and in `warnings`.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1, "y": 2})
        >>> target = {}
        >>> report = session.load(names="x", ns=target)
        >>> target
        {'x': 1}
        >>> report.loaded
        ['x']
        """
        if ns is None:
            ns = _context.default_namespace()
        names = _as_list(names)
        report = self._report("load")
        manifest = self._manifest()
        for message in _context.compare_environment(manifest.environment,
                                                    _context.environment()):
            logger.warning("%s", message)
            report.warnings.append(message)
        for pattern in names or []:
            if not _is_glob(pattern) and pattern not in manifest.variables:
                report.skipped[pattern] = "not in the session"
        for name, entry in manifest.variables.items():
            if names is not None and not _matches(name, names):
                continue
            try:
                ns[name] = self._decode(name, entry)
            except Exception as exc:
                report.failed[name] = _error_text(exc)
                logger.warning("could not load %s: %s", name, _error_text(exc))
                continue
            report.loaded.append(name)
        return report

    def _decode(self, name, entry):
        if entry.sha256:
            digest = hashlib.sha256()
            for rel in entry.files:
                digest.update((self.path / rel).read_bytes())
            if digest.hexdigest() != entry.sha256:
                raise ValueError(f"files of {name!r} are corrupted: sha256 "
                                 f"does not match the manifest")
        prefix = self._objects() / _file_prefix(name)
        meta = dict(entry.meta)
        meta.setdefault("files", [Path(f).name for f in entry.files])
        return codecs.decode(entry.codec, prefix, meta)

    def get(self, name):
        """
        Read one saved value without touching a namespace.

        Parameters
        ----------
        name : str
            Name of the variable.

        Returns
        -------
        object
            The value, decoded from the files of the session.

        Raises
        ------
        KeyError
            If there is no variable `name` in the session.
        ValueError
            If the stored files do not match the sha256 recorded in the
            manifest.

        See Also
        --------
        load : Load several variables into a namespace.
        info : The manifest entry of one variable.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": [1, 2]})
        >>> session.get("x")
        [1, 2]
        """
        entry = self.info_entry(name)
        return self._decode(name, entry)

    def info_entry(self, name):
        """
        Return the manifest entry of one variable.

        Parameters
        ----------
        name : str
            Name of the variable.

        Returns
        -------
        VariableEntry
            What the manifest records for `name`: type, codec, files, size,
            checksum, fingerprint and summary. Reading it unpickles nothing.

        Raises
        ------
        KeyError
            If there is no variable `name` in the session.

        See Also
        --------
        info : The same entry as a dict.
        get : The value of the variable.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1})
        >>> session.info_entry("x").codec
        'pickle'
        """
        manifest = self._manifest()
        try:
            return manifest.variables[name]
        except KeyError:
            raise KeyError(f"no variable {name!r} in session {self.path}") from None

    def info(self, name):
        """
        Return the manifest entry of one variable as a dict.

        Parameters
        ----------
        name : str
            Name of the variable.

        Returns
        -------
        dict
            The manifest entry with the keys "type", "codec",
            "codec_version", "files", "bytes", "sha256", "fingerprint",
            "saved", "version_bound", "summary" and "meta".

        Raises
        ------
        KeyError
            If there is no variable `name` in the session.

        See Also
        --------
        info_entry : The same entry as a `VariableEntry`.
        ls : The entries of all variables as a DataFrame.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1})
        >>> session.info("x")["codec"]
        'pickle'
        """
        return self.info_entry(name).to_dict()

    def ls(self):
        """
        List the saved variables as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            One row per variable, with the columns "name", "type", "codec",
            "bytes", "saved", "version_bound", "kind", "length" and "shape".
            Empty if the session has no variables.

        See Also
        --------
        info : The complete manifest entry of one variable.

        Notes
        -----
        The table comes from the manifest alone: no object file is read.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1, "y": 2})
        >>> sorted(session.ls()["name"])
        ['x', 'y']
        """
        import pandas as pd
        columns = ["name", "type", "codec", "bytes", "saved", "version_bound",
                   "kind", "length", "shape"]
        rows = []
        for name, entry in self._manifest().variables.items():
            summary = entry.summary or {}
            rows.append([name, entry.type, entry.codec, entry.bytes, entry.saved,
                         entry.version_bound, summary.get("kind"),
                         summary.get("length"), summary.get("shape")])
        return pd.DataFrame(rows, columns=columns)

    def remove(self, name):
        """
        Delete one variable and its files from the session.

        Parameters
        ----------
        name : str
            Name of the variable.

        Returns
        -------
        Report
            The report, with `name` in `removed`.

        Raises
        ------
        KeyError
            If there is no variable `name` in the session.
        SessionLocked
            If another process holds the session lock for longer than
            `lock_timeout`.

        See Also
        --------
        save : Add or update variables.
        ls : List the variables that are still in the session.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1})
        >>> report = session.remove("x")
        >>> report.removed
        ['x']
        >>> session.ls().empty
        True
        """
        report = self._report("remove")
        with session_lock(self.path, self.lock_timeout):
            manifest = self._manifest()
            stale = []
            if name in manifest.variables:
                stale = manifest.variables.pop(name).files
            elif name in manifest.failed:
                pass
            else:
                raise KeyError(f"no variable {name!r} in session {self.path}")
            manifest.failed.pop(name, None)
            manifest.revision += 1
            manifest.updated = _now()
            write_manifest(self.path, manifest)
            self._delete_files(stale)
        report.removed.append(name)
        return report

    def diff(self, ns=None, level="names", max_bytes=None):
        """
        Compare a namespace with the saved session.

        Parameters
        ----------
        ns : dict-like, optional
            Namespace to compare. Default: the user namespace, IPython's
            `user_ns` inside IPython and `__main__.__dict__` otherwise.
        level : str, default "names"
            Comparison level. "names" only compares which names exist:
            `added` and `removed` are filled, and a name in both goes to
            `unchanged`, meaning only "present in both". "summary" also
            compares the type and the summary (length, shape, columns,
            dtypes). "content" also compares fingerprints.
        max_bytes : int or None, optional
            Largest live object, in estimated bytes, that "content" hashes.
            A larger object goes to `unknown` and is not hashed. Default:
            None, no limit.

        Returns
        -------
        DiffResult
            The names in `added`, `removed`, `changed`, `unchanged` and
            `unknown`, as described by the level.

        Raises
        ------
        ValueError
            If `level` is not "names", "summary" or "content".

        See Also
        --------
        DiffResult : The result of the comparison.
        save : Write the variables of a namespace into the session.

        Notes
        -----
        The same default filters as `save` apply to `ns`.

        Values saved by pickle whose pickle depends on hash randomization (a
        set or frozenset, alone or directly inside a list, tuple or dict, or
        one level deeper) cannot be compared across processes: they go to
        `unknown`, never `changed`.

        Examples
        --------
        >>> import os
        >>> import tempfile
        >>> from rotifer.core.session import Session
        >>> session = Session(os.path.join(tempfile.mkdtemp(), "s"))
        >>> report = session.save({"x": 1})
        >>> result = session.diff({"x": 1, "y": 2})
        >>> result.added, result.removed
        (['y'], [])
        >>> result.unchanged
        ['x']
        >>> session.diff({"x": 2}, level="content").changed
        ['x']
        """
        if level not in DIFF_LEVELS:
            raise ValueError(f"level must be one of {DIFF_LEVELS}, not {level!r}")
        if ns is None:
            ns = _context.default_namespace()
        manifest = self._manifest()
        live = {n: v for n, v in ns.items()
                if isinstance(n, str) and _skip_reason(n, v) is None}
        result = DiffResult(level=level, path=str(self.path))
        result.added = sorted(set(live) - set(manifest.variables))
        result.removed = sorted(set(manifest.variables) - set(live))
        for name in sorted(set(live) & set(manifest.variables)):
            entry = manifest.variables[name]
            value = live[name]
            if level == "names":
                result.unchanged.append(name)
                continue
            if _type_name(value) != entry.type or \
                    _comparable_summary(_safe_describe(value)) != _comparable_summary(entry.summary or {}):
                result.changed.append(name)
                continue
            if level == "summary":
                result.unchanged.append(name)
                continue
            if max_bytes is not None and _estimated_size(value) > max_bytes:
                result.unknown.append(name)
                continue
            if (entry.fingerprint or "").startswith("pickle:") \
                    and hash_order_dependent(value):
                result.unknown.append(name)
                continue
            try:
                codec = codecs.get_codec(entry.codec)
            except KeyError:
                codec = None
            try:
                current = fingerprint(value, codec)
            except Exception:
                current = None
            if current is None or entry.fingerprint is None:
                result.unknown.append(name)
            elif current != entry.fingerprint:
                result.changed.append(name)
            else:
                result.unchanged.append(name)
        return result


def open_session(path, create=False, lock_timeout=10.0):
    """
    Open an existing session directory.

    Parameters
    ----------
    path : str or pathlib.Path
        Path of the session directory.
    create : bool, default False
        Create the session at `path` if it does not exist. Unlike `Session`,
        it is not created by default.
    lock_timeout : float, default 10.0
        Seconds a writer waits for the session lock before giving up.

    Returns
    -------
    Session
        The open session.

    Raises
    ------
    FileNotFoundError
        If `create` is False and there is no session at `path`.
    ValueError
        If `path` exists and is not a session directory.

    See Also
    --------
    Session : The class this function calls.

    Examples
    --------
    >>> import os
    >>> import tempfile
    >>> from rotifer.core.session import Session, open_session
    >>> path = os.path.join(tempfile.mkdtemp(), "s")
    >>> report = Session(path).save({"x": 1})
    >>> open_session(path).get("x")
    1
    """
    return Session(path, create=create, lock_timeout=lock_timeout)
