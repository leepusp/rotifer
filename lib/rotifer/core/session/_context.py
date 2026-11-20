"""Where the session runs: namespace, environment and notebook context."""

import importlib.metadata
import os
import platform
import sys

PACKAGES = {
    "pandas": "pandas",
    "numpy": "numpy",
    "biopython": "biopython",
    "pyarrow": "pyarrow",
    "cloudpickle": "cloudpickle",
}


def _ipython():
    """The running IPython shell, or None. IPython is imported lazily."""
    try:
        from IPython import get_ipython
    except ImportError:
        return None
    return get_ipython()


def default_namespace():
    """
    Return the user namespace.

    This is IPython's `user_ns` if inside IPython, else `__main__.__dict__`.

    Returns
    -------
    dict
        The namespace dict, not a copy.

    Examples
    --------
    >>> ns = default_namespace()
    >>> type(ns).__name__
    'dict'
    """
    shell = _ipython()
    if shell is not None:
        return shell.user_ns
    return sys.modules["__main__"].__dict__


def _version(distribution):
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def environment():
    """
    Return the versions of Python, ROTIFER and the key packages.

    Returns
    -------
    dict
        Keys "python" and "rotifer" (version strings; "rotifer" is None if
        it is not installed), "packages" (maps each key package to its
        version, or None if not installed) and "tools" (always empty;
        reserved).

    See Also
    --------
    compare_environment : Compare two environment dicts.
    """
    return {
        "python": platform.python_version(),
        "rotifer": _version("rotifer"),
        "packages": {name: _version(dist) for name, dist in PACKAGES.items()},
        "tools": {},
    }


def notebook_context():
    """
    Information about the IPython session, or {} outside IPython.

    Returns
    -------
    dict
        Empty outside IPython. Otherwise "shell" (the shell class name) plus
        the keys `history_session`, `kernel_id`, `connection_file` and
        `notebook`, each present only when it can be found.
    """
    shell = _ipython()
    if shell is None:
        return {}
    info = {"shell": type(shell).__name__}
    history = getattr(shell, "history_manager", None)
    if history is not None and getattr(history, "session_number", None):
        info["history_session"] = int(history.session_number)
    try:
        from ipykernel.connect import get_connection_file
        connection_file = get_connection_file()
        info["connection_file"] = connection_file
        name = os.path.basename(connection_file)
        if name.startswith("kernel-") and name.endswith(".json"):
            info["kernel_id"] = name[len("kernel-"):-len(".json")]
    except Exception:
        pass
    notebook = (os.environ.get("JPY_SESSION_NAME")
                or getattr(shell, "user_ns", {}).get("__vsc_ipynb_file__"))
    if notebook:
        info["notebook"] = str(notebook)
    return info


def compare_environment(saved, current):
    """
    Compare two `environment()` dicts and return one message per difference.

    Parameters
    ----------
    saved : dict
        Environment recorded when the session was saved.
    current : dict
        Environment of the running process.

    Returns
    -------
    list of str
        One message per difference. The list is empty when nothing differs.

    See Also
    --------
    environment : Return the versions of Python, ROTIFER and key packages.
    """
    messages = []
    for key, label in (("python", "Python"), ("rotifer", "ROTIFER")):
        a, b = saved.get(key), current.get(key)
        if a != b:
            messages.append(f"{label} version differs: saved with {a}, running {b}")
    for section in ("packages", "tools"):
        old, new = saved.get(section) or {}, current.get(section) or {}
        for name in sorted(set(old) | set(new)):
            a, b = old.get(name), new.get(name)
            if a == b:
                continue
            if a is None:
                messages.append(f"{name} was not installed when saved, running {b}")
            elif b is None:
                messages.append(f"{name} {a} was used when saved, but it is not installed")
            else:
                messages.append(f"{name} version differs: saved with {a}, running {b}")
    return messages
