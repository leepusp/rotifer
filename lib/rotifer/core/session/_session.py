"""Session: save and load the variables of a working session.

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
    """What a Session call did.

    saved, loaded and removed are lists of names. skipped maps a name to the
    reason it was left out; failed maps a name to "ExcType: message".
    converted_from is the old shelve path if this session was converted from
    one.
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
    """Result of Session.diff.

    added: in the namespace, not in the session.
    removed: in the session, not in the namespace.
    changed: in both, and different at the chosen level.
    unchanged: in both, and no difference found at the chosen level. At level
        "names" it means only "present in both".
    unknown: in both, but it could not be compared (too big for max_bytes, or
        it cannot be fingerprinted).
    """

    level: str
    path: str
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    changed: list = field(default_factory=list)
    unchanged: list = field(default_factory=list)
    unknown: list = field(default_factory=list)

    def to_dict(self):
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
    """A session directory.

    Session(path) opens the session at `path`, creating it if `create` is
    true. If `path` is an old shelve session, it is converted to
    `<path>.session/` and that directory is used; the old files are never
    changed. `converted_from` then holds the old path.

    `lock_timeout` is how many seconds a writer waits for the session lock.
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
        """Convert the old shelve session `old` into the session directory `new`.

        `new` defaults to `legacy_target(old)`. The old files are not
        changed. Returns the new Session; its `conversion_report` says what
        was converted.
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
        """Save variables of `ns` (default: the user namespace) into the session.

        `names` and `exclude` are names or glob patterns. By default a call
        only adds and updates; with `prune=True` every saved variable inside
        `names` (all of them if `names` is None) that this call does not save
        is removed: gone from `ns`, or excluded by `exclude` or by the
        default filters.

        New files are written under names unique to this save and the
        manifest is replaced atomically; only then are the files it no longer
        lists deleted, so an interrupted save leaves the old session intact. Default filters
        skip names that start with "_", the IPython names (In, Out, exit,
        quit, get_ipython, open), modules, and functions or classes not
        defined in `__main__`. Returns a Report.
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
        """Load saved variables into `ns` (default: the user namespace).

        `names` is a name or glob patterns. Warns, in the log and in the
        report, if the saved environment differs from the current one. A
        variable that cannot be read goes to `report.failed`; the others
        are still loaded. Returns a Report.
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
        """Return the saved value of `name` without touching any namespace."""
        entry = self.info_entry(name)
        return self._decode(name, entry)

    def info_entry(self, name):
        manifest = self._manifest()
        try:
            return manifest.variables[name]
        except KeyError:
            raise KeyError(f"no variable {name!r} in session {self.path}") from None

    def info(self, name):
        """Return the manifest entry of `name` as a dict."""
        return self.info_entry(name).to_dict()

    def ls(self):
        """List the saved variables as a DataFrame, from the manifest only."""
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
        """Delete `name` and its files from the session. Returns a Report."""
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
        """Compare the namespace with the saved session.

        level "names": only added and removed; names in both go to
            `unchanged`, meaning just "present in both".
        level "summary": also compares the type and the summary (length,
            shape, columns, dtypes).
        level "content": also compares fingerprints. A live object larger
            than `max_bytes` (estimated) goes to `unknown` and is not hashed.
            Values saved by pickle whose pickle depends on hash randomization
            (a set or frozenset, alone or directly inside a list, tuple or
            dict, or one level deeper) cannot be compared across processes
            and go to `unknown`, never `changed`.
        The same default filters as `save` apply to the namespace.
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
    """Open the session at `path`. Unlike Session(), it does not create it by default."""
    return Session(path, create=create, lock_timeout=lock_timeout)
