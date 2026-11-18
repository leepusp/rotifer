"""Codecs: how each kind of value is written to, and read back from, disk.

A codec writes a value to one or more files that share a path prefix (no
extension) and returns a metadata dict. The dict must include "files", the
names of the files created, relative to the prefix's directory. The same
prefix and dict are later given to the codec's read function.

Built-in codecs:

pickle             stdlib pickle, file ".pkl" (or ".pkl.gz" with compress="gzip").
cloudpickle        for values that depend on code defined in __main__, file
                   ".cloudpickle.pkl" (or ".cloudpickle.pkl.gz"). Flagged
                   version_bound, since it can only be read by the same Python
                   version.
parquet-dataframe  plain pandas DataFrames, file ".parquet". Registered only
                   if pyarrow imports. Checked by reading the file back.

Write functions may accept a "compress" keyword; encode() passes it to those
that do and ignores it for the others.
"""
from __future__ import annotations

import gzip
import importlib
import inspect
import pickle
import pickletools
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from rotifer.core.logger import get_logger
from rotifer.core.session._fingerprint import fingerprint as _fingerprint
from rotifer.core.session._fingerprint import pickle_bytes_fingerprint

logger = get_logger(__name__)


@dataclass(frozen=True)
class Codec:
    """A way of saving one family of values.

    write(obj, prefix) creates files named prefix + suffix and returns a
    metadata dict with "files". read(prefix, meta) returns the object.
    Lookup is by exact type; only a codec with subclasses=True also handles
    subclasses of its types.
    """
    name: str
    version: int
    types: tuple
    write: Callable[[object, Path], dict]
    read: Callable[[Path, dict], object]
    subclasses: bool = False
    describe: Callable[[object], dict] | None = None
    fingerprint: Callable[[object], str] | None = None


_codecs: dict[str, Codec] = {}
_by_type: dict[type, str] = {}
_aliases: dict[str, str] = {}


# ---------------------------------------------------------------- registry

def register_codec(codec: Codec, *, replace: bool = False) -> None:
    """Add codec to the registry.

    Raises ValueError if the name, or one of its exact types, belongs to
    another codec, unless replace=True (then the new codec takes them over).
    """
    if codec.name in _codecs and not replace:
        raise ValueError("codec %r is already registered" % codec.name)
    for t in codec.types:
        owner = _by_type.get(t)
        if owner is not None and owner != codec.name and not replace:
            raise ValueError("type %s is already handled by codec %r" % (t.__name__, owner))
    old = _codecs.get(codec.name)
    if old is not None:
        for t in old.types:
            if _by_type.get(t) == old.name:
                del _by_type[t]
    _codecs[codec.name] = codec
    for t in codec.types:
        _by_type[t] = codec.name


_providers: dict[str, str] = {"rotifer.sequence": "rotifer.devel.beta.sequence"}


def register_codec_provider(codec_name: str, module_name: str) -> None:
    """Say which module registers codec_name when imported.

    get_codec() imports that module the first time the codec is missing, so a
    session can be read in a process that never imported the module.
    """
    _providers[codec_name] = module_name


def get_codec(name: str) -> Codec:
    """Return the codec registered under name (used when reading).

    If the codec is unknown but a provider module is known for it, the module
    is imported and the lookup retried.
    """
    if name not in _codecs and name in _providers:
        try:
            importlib.import_module(_providers[name])
        except ImportError as exc:
            logger.warning("could not import %s, which provides codec %r: %s",
                           _providers[name], name, exc)
    try:
        return _codecs[name]
    except KeyError:
        raise KeyError("no codec named %r is registered (known: %s)"
                       % (name, ", ".join(sorted(_codecs)))) from None


def register_alias(old: str, new: str) -> None:
    """Make pickles that mention class path old load class path new.

    Both are "module.QualName" strings. Used when a class moves.
    """
    _aliases[old] = new


def _resolve_alias(module: str, name: str):
    new = _aliases.get(module + "." + name)
    if new is None:
        return None
    parts = new.split(".")
    for i in range(len(parts) - 1, 0, -1):
        try:
            mod = importlib.import_module(".".join(parts[:i]))
        except ImportError:
            continue
        obj = mod
        try:
            for attr in parts[i:]:
                obj = getattr(obj, attr)
        except AttributeError:
            continue
        return obj
    raise ImportError("alias target %r for %s.%s cannot be imported" % (new, module, name))


class _AliasUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        found = _resolve_alias(module, name)
        if found is not None:
            return found
        return super().find_class(module, name)


def _load_pickle(fileobj):
    return _AliasUnpickler(fileobj).load()


# ------------------------------------------------------- __main__ detection

def _uses_main(data: bytes) -> bool:
    """True if the pickle bytes reference any name in the __main__ module."""
    memo: dict = {}
    recent: list = []     # values of the latest string pushes, newest last
    string_ops = ("SHORT_BINUNICODE", "BINUNICODE", "BINUNICODE8", "UNICODE")
    for op, arg, _pos in pickletools.genops(data):
        name = op.name
        if name in string_ops:
            recent.append(arg)
        elif name == "GLOBAL":
            if str(arg).split(" ", 1)[0] == "__main__":
                return True
        elif name == "STACK_GLOBAL":
            if len(recent) >= 2 and recent[-2] == "__main__":
                return True
        elif name == "MEMOIZE":
            if recent:
                memo[len(memo)] = recent[-1]
        elif name in ("PUT", "BINPUT", "LONG_BINPUT"):
            if recent:
                memo[arg] = recent[-1]
        elif name in ("GET", "BINGET", "LONG_BINGET"):
            recent.append(memo.get(arg))
        if len(recent) > 4:
            del recent[:-4]
    return False


def _choose_fallback_data(obj) -> tuple[str, bytes | None]:
    """Return ("pickle" or "cloudpickle", stdlib pickle bytes or None)."""
    try:
        data = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception:
        return "cloudpickle", None
    return ("cloudpickle" if _uses_main(data) else "pickle"), data


def _choose_fallback(obj) -> str:
    """Return "pickle" or "cloudpickle" for a value no typed codec handles."""
    return _choose_fallback_data(obj)[0]


# ------------------------------------------------------------- file helpers

def _atomic(path: Path, mode: str = "wb"):
    from rotifer.core.session._files import atomic_writer
    return atomic_writer(path, mode)


def _path(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def _dump_pickle(obj, prefix: Path, suffix: str, dumps, compress, data=None) -> dict:
    if compress not in (None, "gzip"):
        raise ValueError("compress must be None or 'gzip', not %r" % (compress,))
    path = _path(prefix, suffix + (".gz" if compress else ""))
    if data is None:
        data = dumps(obj)
    with _atomic(path, "wb") as f:
        if compress:
            with gzip.GzipFile(fileobj=f, mode="wb", mtime=0) as g:
                g.write(data)
        else:
            f.write(data)
    return {"files": [path.name], "compress": compress}


def _read_pickle(prefix: Path, meta: dict):
    path = Path(prefix).parent / meta["files"][0]
    opener = gzip.open if meta.get("compress") == "gzip" else open
    with opener(path, "rb") as f:
        return _load_pickle(f)


# ----------------------------------------------------------- built-in codecs

def _write_pickle(obj, prefix: Path, compress=None, _data=None) -> dict:
    meta = _dump_pickle(obj, prefix, ".pkl",
                        lambda o: pickle.dumps(o, protocol=pickle.HIGHEST_PROTOCOL),
                        compress, _data)
    meta.update(protocol=pickle.HIGHEST_PROTOCOL, version_bound=False)
    if _data is not None:       # the bytes written are the stdlib pickle
        meta["fingerprint"] = pickle_bytes_fingerprint(_data)
    return meta


def _write_cloudpickle(obj, prefix: Path, compress=None, _data=None) -> dict:
    try:
        import cloudpickle
    except ImportError as exc:
        raise ImportError("cloudpickle is needed to save values defined in __main__; "
                          "install it with 'pip install cloudpickle'") from exc
    meta = _dump_pickle(obj, prefix, ".cloudpickle.pkl",
                        lambda o: cloudpickle.dumps(o, protocol=pickle.HIGHEST_PROTOCOL), compress)
    meta.update(protocol=pickle.HIGHEST_PROTOCOL, version_bound=True,
                cloudpickle=cloudpickle.__version__)
    if _data is not None:       # stdlib pickle of the same object, already made
        meta["fingerprint"] = pickle_bytes_fingerprint(_data)
    return meta


def _write_parquet(df, prefix: Path) -> dict:
    import pandas as pd
    path = _path(prefix, ".parquet")
    try:
        with _atomic(path, "wb") as f:
            df.to_parquet(f, engine="pyarrow")
        back = pd.read_parquet(path, engine="pyarrow")
        if not (df.equals(back) and df.columns.equals(back.columns)
                and list(df.dtypes) == list(back.dtypes) and df.index.equals(back.index)
                and df.columns.names == back.columns.names
                and df.index.names == back.index.names
                and df.index.dtype == back.index.dtype
                and getattr(df.index, "freq", None) == getattr(back.index, "freq", None)
                and getattr(df.columns, "freq", None) == getattr(back.columns, "freq", None)):
            raise _RoundTripMismatch("round-trip mismatch")
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return {"files": [path.name], "engine": "pyarrow"}


class _RoundTripMismatch(Exception):
    pass


def _read_parquet(prefix: Path, meta: dict):
    import pandas as pd
    return pd.read_parquet(Path(prefix).parent / meta["files"][0], engine="pyarrow")


def _builtin_codecs():
    yield Codec("pickle", 1, (), _write_pickle, _read_pickle)
    yield Codec("cloudpickle", 1, (), _write_cloudpickle, _read_pickle)
    try:
        import pandas as pd
        import pyarrow  # noqa: F401
    except ImportError:
        logger.debug("pyarrow or pandas missing; parquet-dataframe codec not registered")
        return
    yield Codec("parquet-dataframe", 1, (pd.DataFrame,), _write_parquet, _read_parquet,
                subclasses=False,
                fingerprint=lambda df: _fingerprint(df))


for _c in _builtin_codecs():
    register_codec(_c)


# -------------------------------------------------------------- selection

def find_typed_codec(obj) -> Codec | None:
    """Return the typed codec for obj, or None. Never pickles obj.

    Exact type match first; then codecs with subclasses=True, following the
    type's MRO.
    """
    name = _by_type.get(type(obj))
    if name is not None:
        return _codecs[name]
    for cls in type(obj).__mro__[1:]:
        name = _by_type.get(cls)
        if name is not None and _codecs[name].subclasses:
            return _codecs[name]
    return None


def lookup(obj) -> Codec:
    """Return the codec that would be tried first for obj. Never fails.

    The typed codec if there is one (see find_typed_codec); otherwise the
    fallback ("cloudpickle" if obj depends on __main__, else "pickle"). The
    fallback check pickles obj, so use find_typed_codec where that is too slow.
    """
    codec = find_typed_codec(obj)
    if codec is not None:
        return codec
    return _codecs[_choose_fallback(obj)]


def _call_write(codec: Codec, obj, prefix: Path, compress, std_data=None) -> dict:
    kwargs = {}
    if std_data is not None and codec.name in ("pickle", "cloudpickle"):
        kwargs["_data"] = std_data
    if compress is not None and "compress" in inspect.signature(codec.write).parameters:
        kwargs["compress"] = compress
    meta = dict(codec.write(obj, prefix, **kwargs))
    if "files" not in meta:
        raise ValueError("codec %r did not report its files" % codec.name)
    meta.setdefault("version_bound", False)
    return meta


def encode(obj, prefix, *, compress=None) -> tuple[Codec, dict]:
    """Write obj under prefix with the best codec; return (codec, metadata).

    compress=None or "gzip" applies to the pickle codecs. If the chosen
    typed codec fails (an exception, or its read-back check), its files are
    removed and obj is saved with pickle or cloudpickle instead; the
    metadata then has "fallback_from" and "fallback_reason".
    """
    prefix = Path(prefix)
    codec = find_typed_codec(obj)
    std_data = None
    if codec is None:
        name, std_data = _choose_fallback_data(obj)
        codec = _codecs[name]
    first_error = None
    if codec.name not in ("pickle", "cloudpickle"):
        try:
            return codec, _call_write(codec, obj, prefix, compress)
        except Exception as exc:
            first_error = (codec.name, exc)
            logger.info("codec %s failed for %s (%s); falling back to pickle",
                        codec.name, type(obj).__name__, exc)
            name, std_data = _choose_fallback_data(obj)
            codec = _codecs[name]
    try:
        meta = _call_write(codec, obj, prefix, compress, std_data)
    except Exception as exc:
        if codec.name == "pickle":    # e.g. an unpicklable lambda
            codec = _codecs["cloudpickle"]
            meta = _call_write(codec, obj, prefix, compress, None)
        else:
            raise
    if first_error is not None:
        name, exc = first_error
        meta["fallback_from"] = name
        meta["fallback_reason"] = str(exc) if isinstance(exc, _RoundTripMismatch) \
            else "%s: %s" % (type(exc).__name__, exc)
    return codec, meta


def decode(codec_name: str, prefix, meta: dict):
    """Read back a value written by encode()."""
    return get_codec(codec_name).read(Path(prefix), meta)
