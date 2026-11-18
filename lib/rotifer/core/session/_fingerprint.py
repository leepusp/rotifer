"""Cheap, repeatable fingerprints of values, used to tell whether a variable changed.

A fingerprint is a string "<method>:<hex digest>". Equal values give equal
fingerprints within the same library versions; it is a change detector, not a
security feature.
"""
from __future__ import annotations

import hashlib
import pickle
from typing import TYPE_CHECKING

from rotifer.core.logger import get_logger

if TYPE_CHECKING:
    from rotifer.core.session.codecs import Codec

logger = get_logger(__name__)


def _new_hash():
    return hashlib.blake2b(digest_size=20)


def _pandas_fingerprint(obj) -> str:
    import pandas as pd

    h = _new_hash()
    h.update(type(obj).__name__.encode())
    if isinstance(obj, pd.DataFrame):
        h.update(repr(list(obj.columns)).encode())
        h.update(repr([str(t) for t in obj.dtypes]).encode())
    else:
        h.update(repr(obj.name).encode())
        h.update(str(obj.dtype).encode())
    h.update(repr(obj.shape).encode())
    h.update(pd.util.hash_pandas_object(obj, index=True).values.tobytes())
    return "pandas:" + h.hexdigest()


def _numpy_fingerprint(obj) -> str:
    import numpy as np

    h = _new_hash()
    h.update(str(obj.dtype).encode())
    h.update(repr(obj.shape).encode())
    h.update(np.ascontiguousarray(obj).tobytes())
    return "numpy:" + h.hexdigest()


def pickle_bytes_fingerprint(data: bytes) -> str:
    """The "pickle:" fingerprint of already pickled bytes."""
    h = _new_hash()
    h.update(data)
    return "pickle:" + h.hexdigest()


def _pickle_fingerprint(obj) -> str:
    return pickle_bytes_fingerprint(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))


def _has_set(obj) -> bool:
    return isinstance(obj, (set, frozenset))


def _children(obj):
    if isinstance(obj, dict):
        return list(obj.keys()) + list(obj.values())
    if isinstance(obj, (list, tuple)):
        return list(obj)
    return []


def hash_order_dependent(obj) -> bool:
    """True if the pickle of obj may change between processes.

    Sets and frozensets pickle in an order that depends on hash randomization
    of their members. Checks obj, its items, and the items of those items.
    Other containers (custom classes, deeper nesting) are not inspected.
    """
    if _has_set(obj):
        return True
    for child in _children(obj):
        if _has_set(child):
            return True
        if any(_has_set(c) for c in _children(child)):
            return True
    return False


def fingerprint(obj, codec: "Codec | None" = None) -> str | None:
    """Return a fingerprint string for obj, or None if it cannot be hashed.

    Order of methods: the codec's own fingerprint function, then pandas
    hashing for DataFrame and Series, then raw bytes for numpy arrays, then
    a hash of the stdlib pickle bytes. A "pickle:" fingerprint of a value that
    contains sets is not comparable across processes (hash randomization);
    Session.diff reports such variables as unknown.
    """
    if codec is not None and codec.fingerprint is not None:
        try:
            return codec.fingerprint(obj)
        except Exception as exc:
            logger.debug("codec %s fingerprint failed: %s", codec.name, exc)
    methods = []
    try:
        import pandas as pd
        if isinstance(obj, (pd.DataFrame, pd.Series)):
            methods.append(_pandas_fingerprint)
    except ImportError:
        pass
    try:
        import numpy as np
        if isinstance(obj, np.ndarray) and obj.dtype != object:
            methods.append(_numpy_fingerprint)
    except ImportError:
        pass
    methods.append(_pickle_fingerprint)
    for method in methods:
        try:
            return method(obj)
        except Exception as exc:
            logger.debug("%s failed for %s: %s", method.__name__, type(obj).__name__, exc)
    return None
