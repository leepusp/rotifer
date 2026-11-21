# Sessions

The `rotifer.core.session` package saves variables from a Python or IPython
namespace in a directory. A JSON manifest describes the saved values, and
each value has its own file or files under `objects/`. This replaces the
single shelve file used by the old session functions. You can inspect,
update or remove one variable without touching the others, and failures come
back to the caller. The full reference is in the generated API pages linked
below. The {doc}`session manager design record
<../adr/0001-session-manager>` explains the storage and compatibility
choices.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:device-floppy` Session management
:link: /api/generated/rotifer.core.session
:link-type: doc

Save, restore, inspect, compare and remove variables in a session directory.
:::

:::{grid-item-card} {iconify}`tabler:history` Compatibility wrappers
:link: /api/generated/rotifer.core.functions
:link-type: doc

The `save_session` and `restore_session` wrappers for existing callers.
:::
::::

## Quickstart

Pass a dictionary as `ns` to save or restore an explicit namespace. If `ns`
is omitted, sessions use IPython's user namespace inside IPython and the
`__main__` namespace otherwise.

::::{tab-set}

:::{tab-item} Save and load

`Session` creates its directory by default. Saving adds and updates variables
without removing others. Loading can restore the whole session or a subset.

```python
import os
import tempfile

from rotifer.core.session import Session

session = Session(os.path.join(tempfile.mkdtemp(), "s"))
report = session.save({"x": 1, "y": 2})
assert sorted(report.saved) == ["x", "y"]
assert not report.failed

target = {}
report = session.load(names="x", ns=target)
assert target == {"x": 1}
assert session.get("x") == 1
```

When saving or diffing without an explicit `names` selection, the default
filters skip names starting
with `_`, the IPython names `In`, `Out`, `exit`, `quit`, `get_ipython` and
`open`, modules, and functions or classes not defined in `__main__`.
:::

:::{tab-item} Inspect

`ls()` returns a DataFrame built from the manifest alone. `info()` returns the
complete manifest entry for one variable. Neither reads the object files.
`diff()` compares names by default; use `level="summary"` to compare types and
summaries, or `level="content"` to compare content fingerprints.

```python
import os
import tempfile

from rotifer.core.session import Session

session = Session(os.path.join(tempfile.mkdtemp(), "s"))
session.save({"x": 1, "y": 2})

listing = session.ls()
assert listing["name"].tolist() == ["x", "y"]
assert session.info("x")["codec"] == "pickle"

result = session.diff(ns={"x": 1, "z": 3})
assert result.added == ["z"]
assert result.removed == ["y"]
assert result.unchanged == ["x"]
```
:::

:::{tab-item} Convert an old session

Pass the path of an old shelve session to `Session.convert()`. The conversion
writes `<old>.session/` by default and leaves the old files untouched.

```python
import os
import shelve
import tempfile
import warnings

from rotifer.core.session import Session

old = os.path.join(tempfile.mkdtemp(), "old")
with shelve.open(old) as db:
    db["n"] = 5

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    session = Session.convert(old)

assert any(issubclass(item.category, FutureWarning) for item in caught)
assert session.path.name == "old.session"
assert session.get("n") == 5
assert not session.conversion_report.failed
```
:::
::::

## Files on disk

A session is a plain directory. Its manifest has format `rotifer-session` and
version 1. It can be read as JSON without unpickling any saved value. Run
`python -m rotifer.core.session schema` to print its JSON Schema.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Path
  - Written by
  - Content
  - Notes
* - `<dir>/manifest.json`
  - Session initialization, saves, removals and conversions
  - Format and version, revision, environment, variable entries and save failures
  - Replaced atomically. Variable entries record type, codec, files, size, checksum, fingerprint, timestamp and summary.
* - `<dir>/objects/<name>-<hash>.r<revision>.<ext>`
  - The codec selected for each saved variable
  - The serialized value, in one or more files
  - The name hash and revision make file names unique. Loading checks the sha256 recorded in the manifest.
* - `<dir>/.lock`
  - Saves, removals and conversions
  - Writer coordination
  - Only writers take the lock. A writer raises `SessionLocked` after `lock_timeout`; readers do not lock.
```

## Codecs

The manifest records the codec that actually wrote each variable, including
when a typed codec falls back to pickle.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Name
  - Handles
  - Written as
  - Notes
* - `parquet-dataframe`
  - Plain pandas DataFrames
  - `.parquet`
  - Registered when pandas and pyarrow can be imported. Every file is read back and compared with the original; a failure falls back to a pickle codec.
* - `pickle`
  - Values without a typed codec that do not depend on `__main__`
  - `.pkl`, or `.pkl.gz` with `compress="gzip"`
  - The fallback, using the standard library pickle.
* - `cloudpickle`
  - Functions, classes and instances that depend on definitions in `__main__`, and values the standard pickle cannot write
  - `.cloudpickle.pkl`, or `.cloudpickle.pkl.gz` with `compress="gzip"`
  - Marked `version_bound` because the output depends on the Python version.
```

Custom codecs are added with
`rotifer.core.session.codecs.register_codec`. Selection is by exact type unless
a codec explicitly handles subclasses.

## Common patterns

Use `prune=True` to make the saved selection mirror a namespace. Variables in
the selection that the call does not save are removed. Use `names` to limit
both saving and pruning to a subset.

```python
import os
import tempfile

from rotifer.core.session import Session

session = Session(os.path.join(tempfile.mkdtemp(), "s"))
session.save({"a": 1, "b": 2})
report = session.save({"a": 10}, prune=True)

assert report.removed == ["b"]
assert session.ls()["name"].tolist() == ["a"]
assert session.get("a") == 10
```

The `names` and `exclude` arguments accept one name, one glob pattern or a
sequence of names and patterns. The same glob syntax works when loading.

```python
import os
import tempfile

from rotifer.core.session import Session

session = Session(os.path.join(tempfile.mkdtemp(), "s"))
report = session.save(
    {"a1": 1, "a2": 2, "b1": 3},
    names=["a*", "b1"],
    exclude="a2",
)
assert sorted(report.saved) == ["a1", "b1"]

target = {}
session.load(names="a*", ns=target)
assert target == {"a1": 1}
```

The package sends messages to the `rotifer` logger instead of printing.
`save`, `load` and `remove` return a `Report`; a conversion stores one in
`session.conversion_report`. The report's `action` and `path` identify the
operation. `saved`, `loaded` and `removed` list successful names, while
`skipped` maps names to reasons and `failed` maps names to errors. It also
holds environment `warnings` and the `converted_from` path when relevant.
Check `report.failed` after saving or loading. `report.to_dict()` returns a
JSON-serializable dictionary.

```python
import os
import tempfile

from rotifer.core.session import Session


class Unsaveable:
    def __reduce__(self):
        raise RuntimeError("cannot save this value")


session = Session(os.path.join(tempfile.mkdtemp(), "s"))
report = session.save(
    {"good": 1, "bad": Unsaveable(), "_private": 2},
    names=["good", "bad", "missing", "_private"],
)

assert report.saved == ["good"]
assert "bad" in report.failed
assert report.skipped["missing"] == "not in the namespace"
assert report.skipped["_private"] == "name starts with underscore"
assert session.get("good") == 1
```

Loading a pickle or cloudpickle file can execute code. Only load sessions from
sources you trust. Reading the manifest runs nothing.

## Migrating

Opening an old shelve path with `Session(old_path)` converts it to
`<old>.session/` and opens that directory. `Session.convert(old, new=None)`
does the same and lets you choose the destination. When the conversion runs,
both emit a `FutureWarning`, leave the old files untouched and record the
outcome in `conversion_report`. If `<old>.session/` is already complete,
`Session(old_path)` just opens it. Use the new session directory from then on. The legacy
reader is scheduled for removal on 2026-11-30.

`rotifer.core.functions.save_session` and `restore_session` remain as thin
wrappers over `Session`, and they now write and read session directories.
`save_session` uses `prune=True` to preserve its old
whole-namespace behavior. New code should use `Session` directly.

```{toctree}
:hidden:

../adr/0001-session-manager
```
