# ADR 0001: Session manager (rotifer.core.session)

Status: Proposed

Date: 2026-11-18

Issue: #139

## Context

`save_session()` and `restore_session()` in `rotifer.core.functions` store the whole `__main__` namespace in one `shelve` file. That design has these problems:

- To see, copy or delete one variable, the whole store has to be loaded.
- The file grows and nothing says which variable is responsible.
- Saving fails on systems that use the gdbm backend. Both functions open the file twice, and the second open fails with `EAGAIN`.
- If a variable fails to save, its previous value stays in the file and comes back on restore without a warning.
- Each dbm backend writes different files, so a session does not move reliably between machines.
- Failures are printed to stderr, and the caller never sees them.

More interfaces are planned on top of sessions: IPython magics, a command-line interface, a live variable viewer in ROTIFER's own interface, and extensions for code editors. If each one reads and writes sessions its own way, they will drift apart. They need one library to build on, and that library must not depend on any of them.

## Decision

The shelve code is replaced by a new package, `rotifer.core.session`. The old functions stay as wrappers over it.

### Storage layout

A session is a directory with a `manifest.json` and one or more files per variable under `objects/`, named `<name>-<hash>.r<revision>.<ext>`. The hash suffix prevents collisions on filesystems that ignore case, and the revision suffix gives every save its own file names. A save writes the new object files, then replaces the manifest, and only then deletes the files the manifest no longer lists. A crash or an interrupt at any point leaves the previous manifest pointing at files that still exist; leftover files are removed by the next successful save. Each file is written under a temporary name and moved into place with `os.replace`, and `load` and `get` check the sha256 recorded in the manifest. A lock file (`filelock`) keeps two writers out of the same session. Readers do not lock.

### Entry point

`Session` is the public way to work with a session. It has `save`, `load`, `get`, `ls`, `info`, `remove`, `diff` and a `convert` classmethod, and every method returns a report with `to_dict()`. `save` adds and updates by default; `prune=True` makes the session mirror the namespace. The `save_session` wrapper passes `prune=True`, which keeps its old behavior.

By default `save` skips private names, modules, IPython's own names, and functions or classes imported from elsewhere. Functions and classes defined in `__main__` are saved.

### Namespaces

Any dict-like object works as a namespace, so there are no namespace adapters. `default_namespace()` returns IPython's user namespace when IPython is running and `__main__.__dict__` otherwise. IPython is imported lazily, only for that check.

### Codecs

A codec turns one type of value into files and back. The registry picks a codec by exact type, and a codec covers subclasses only if it opts in, so an `MSA` stays on pickle.

Pickle is the fallback. Values that depend on something defined in `__main__` (functions, classes and their instances) go through cloudpickle, because plain pickle stores them by reference and they would not load in a new process. The manifest flags those entries `version_bound`, since cloudpickle output depends on the Python version.

Plain DataFrames are written as Parquet. The Parquet codec reads the file back and compares it with the original. If writing fails or the copy differs, the variable is saved with pickle, and the manifest records the codec that was used. The Parquet codec is registered only when pyarrow imports.

The codec for the sequence class is defined in that class's module, `rotifer/devel/beta/sequence.py`, so it moves with the class. It writes `df` (and `ndf` when present) as Parquet, falling back to pickle for a table Parquet cannot hold, plus a JSON sidecar with the plain attributes and the checksum, which is verified on load. It stores no module path, and cursors are rebuilt from configuration on load. A session loaded in a new process finds this codec through a small table of codec providers in `codecs.py`, which imports the module the first time the codec is needed.

### Storage backend

There is no `Store` abstraction. One private module, `_files.py`, handles the local disk. If remote or zip storage is ever needed, fsspec is the likely route.

### Rules for the core

The library never writes to stdout and never prompts. Results come back in reports, and messages go through `rotifer.core.logger.get_logger`, which writes to stderr. The core never imports interface code.

### Environment

The manifest records the Python version, the ROTIFER version and the versions of pandas, numpy, biopython, pyarrow and cloudpickle. Loading warns when they differ from the running environment. A `tools` field is reserved for a future tool manager and stays empty for now.

### Diff

`diff(level="names" | "summary" | "content", max_bytes=None)` compares a namespace with a saved session. `"names"` is the default. `"summary"` compares types and summaries from the manifest. `"content"` compares fingerprints: each one is stored at save time, so only the live value is hashed, and values larger than `max_bytes` are reported as unknown. Fingerprints use `pd.util.hash_pandas_object` for DataFrames, raw bytes for numpy arrays, a codec's own function when it has one, and a hash of the pickled value otherwise.

### Legacy sessions

Nothing writes shelve files anymore. A private `_legacy.py` reads old shelve files without modifying them. Opening an old shelve path with `Session` converts it into `<path>.session/` and leaves the old files untouched. The conversion is written to a temporary directory and renamed into place only when it is complete, and reading one emits a `FutureWarning` that asks users to convert. `_legacy.py` will be removed on 2026-11-30, tracked in a separate issue.

### Plugins

Codec registration through entry points (group `rotifer.session.codecs`) waits for the pyproject packaging in #135. Until then, codecs register in code.

## Consequences

What this gives:

- Each variable has its own files, so users can see sizes and copy or delete one variable without loading the rest.
- The manifest can be read without unpickling anything. A viewer, a CLI or an editor extension can list a session with summaries, shapes and sizes from the manifest alone.
- Saving works on gdbm hosts, and a variable that fails to save is reported in the manifest and in the report, with no old value left behind.
- Sessions are plain directories that can be copied between machines.
- Magics, a CLI and a viewer can be added as separate layers, because the core has no interface dependencies.
- DataFrames, the most common large value, are stored in a format other tools can read.

What it costs:

- Loading a pickle executes code, so a session should only be loaded from a trusted source. Editor extensions should gate loading on workspace trust. Reading the manifest executes nothing.
- Values saved with cloudpickle are tied to the Python version that wrote them and may fail to load in another one.
- Parquet cannot hold every DataFrame (mixed-type object columns, for example). Those fall back to pickle, so two DataFrames in one session can use different formats.
- The manifest schema and the codec contract are now public contracts. Both are versioned from v1, and changes need an upgrade path.
- The sequence class is expected to move out of `rotifer.devel` (#150). Pickles written before the move name the old module path, so the old-path shim has to stay until aliases for the new path are registered in the pickle codec.
- `filelock` and `cloudpickle` become dependencies, and #135 has to carry them into the new packaging.
- The legacy reader is temporary code with a fixed removal date.

## Appendix: Manifest v1

The file is `<dir>/manifest.json`. `python -m rotifer.core.session schema` prints its JSON Schema, which is also available as a Python dict from `manifest_schema()`.

```
Session fields
  format        identifier of the manifest format ("rotifer-session")
  version       manifest version (1)
  revision      integer, goes up on every save
  created       timestamp
  updated       timestamp
  environment
    python      Python version
    rotifer     ROTIFER version
    packages    {name: version} for pandas, numpy, biopython, pyarrow, cloudpickle
    tools       {} reserved for the tool manager
  context       {} notebook or kernel context

variables       {name: entry}
  entry fields
    type            qualified type name
    codec           codec actually used
    codec_version   integer
    files           list of files under objects/
    bytes           size on disk
    sha256          checksum of the stored files
    fingerprint     content fingerprint used by diff
    saved           timestamp
    version_bound   true for cloudpickle values
    summary
      kind, length, shape, columns, dtypes, preview
    meta            codec-specific metadata

failed          {name: error} for variables that could not be saved
```

On disk:

```
<dir>/manifest.json
<dir>/objects/<name>-<hash>.r<revision>.<ext>
```

## Appendix: Codec contract

```python
@dataclass(frozen=True)
class Codec:
    name: str
    version: int
    types: tuple[type, ...]
    write: Callable[[object, Path], dict]   # owns every file under the prefix; returns its metadata
    read: Callable[[Path, dict], object]
    subclasses: bool = False                # exact-type lookup unless the codec opts in
    describe: Callable[[object], dict] | None = None
    fingerprint: Callable[[object], str] | None = None
```

Selection is by exact type first. A codec covers subclasses only when `subclasses=True`. Values that depend on `__main__` go to cloudpickle, and anything else without a matching codec goes to pickle.
