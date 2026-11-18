"""describe(): a cheap summary of an in-memory object.

Reads only attributes of objects already in memory (shape, columns, len).
It never reads files and never converts data.
"""

import reprlib

PREVIEW_CHARS = 200
MAX_COLUMNS = 50

_repr = reprlib.Repr()
_repr.maxstring = 80
_repr.maxother = 80
_repr.maxlist = _repr.maxtuple = _repr.maxset = _repr.maxdict = 6


def _preview(obj):
    try:
        text = _repr.repr(obj)
    except Exception:
        text = f"<{type(obj).__name__}>"
    return text[:PREVIEW_CHARS]


def _type_name(obj):
    cls = type(obj)
    return f"{cls.__module__}.{cls.__qualname__}"




def _builtin_summary(obj):
    """Summary for known in-memory types, or None. Exact builtin types only
    for containers, so subclasses with custom __len__/__repr__ are not called."""
    cls = type(obj)
    if cls in (bool, int, float, complex, str, bytes, type(None)):
        info = {"kind": "scalar"}
        if cls in (str, bytes):
            info["length"] = len(obj)
        return info
    if cls in (list, tuple, set, frozenset):
        return {"kind": "sequence", "length": len(obj)}
    if cls is dict:
        return {"kind": "mapping", "length": len(obj)}
    import sys
    pd = sys.modules.get("pandas")
    if pd is not None and isinstance(obj, pd.DataFrame):
        return {
            "kind": "dataframe",
            "length": len(obj),
            "shape": [int(n) for n in obj.shape],
            "columns": [str(c) for c in list(obj.columns)[:MAX_COLUMNS]],
            "dtypes": {str(c): str(t) for c, t in
                       list(obj.dtypes.items())[:MAX_COLUMNS]},
        }
    if pd is not None and isinstance(obj, pd.Series):
        return {"kind": "series", "length": len(obj),
                "shape": [int(n) for n in obj.shape],
                "dtypes": {str(obj.name): str(obj.dtype)}}
    np = sys.modules.get("numpy")
    if np is not None and isinstance(obj, np.ndarray):
        return {"kind": "ndarray", "shape": [int(n) for n in obj.shape],
                "length": int(obj.shape[0]) if obj.ndim else None,
                "dtypes": {"": str(obj.dtype)}}
    return None


def describe(obj):
    """
    Return a summary dict of `obj`: kind, and when they apply length, shape,
    columns (at most 50), dtypes and preview (at most 200 characters).

    For types without a built-in summary, a codec registered for the type may
    supply its own describe function; it is found without pickling `obj`.
    """
    info = _builtin_summary(obj)
    if info is not None:
        info["preview"] = _preview(obj)
        return info
    try:
        from rotifer.core.session import codecs
        codec = codecs.find_typed_codec(obj)
        custom = getattr(codec, "describe", None)
        if custom is not None:
            result = dict(custom(obj))
            result.setdefault("preview", "<%s object>" % _type_name(obj))
            return result
    except Exception:
        pass  # a failing codec must not break the summary
    name = _type_name(obj)
    return {"kind": "object", "type": name, "preview": "<%s object>" % name}
