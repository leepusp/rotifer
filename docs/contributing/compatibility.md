# Compatibility and deprecation

For contributors who change existing code, and for users who want to know what they can rely on. MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

## What counts as public API

The following are public, and changing them is a compatibility matter:

- Functions and classes importable from `rotifer.*` whose names do not start with an underscore.
- Command line tools in `bin/`, and their options.
- Output file formats and their column names.

Anything with a leading underscore, and anything not importable from a documented module path, is internal and MAY change in any release.

`lib/rotifer/devel/` was never public API. Code there could change or disappear at any time, and the transition away from it is described in {doc}`personal-codebases`. Personal-module functions that are removed and still had users get one release of deprecated alias, as described below.

## Deprecation cycle

Removing or incompatibly changing public API takes three steps.

1. Keep the old behaviour working and emit a warning that names the replacement.
2. Add a `deprecated` changelog fragment that names the replacement (see [Changelog](https://leepusp.github.io/rotifer/contributing/changelog.html)).
3. Remove it in the next release, with a `removed` fragment.

Use `DeprecationWarning` for code that other code calls. Use `FutureWarning` when the change affects what a user of a command line tool or a notebook sees. Set `stacklevel=2` so the warning points at the caller and not at the library.

A renamed function keeps its old name as an alias:

```python
import warnings


def merge_hits(hits):
    """Merge overlapping hits."""
    ...


def mergeHits(*args, **kwargs):
    """Deprecated alias of :func:`merge_hits`."""
    warnings.warn(
        "mergeHits is deprecated and will be removed in a future release, "
        "use merge_hits instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return merge_hits(*args, **kwargs)
```

The same applies to a module that moves: leave a small module at the old path that imports from the new one and warns when imported.

## How long a deprecation lasts

A deprecated item spends exactly one release deprecated and is removed in the next release. If it is deprecated in 0.3.0, it is still present, with its warning, in 0.3.0 and gone in 0.4.0. PATCH releases do not count, since they only fix bugs.

Only in extraordinary circumstances MAY the maintainers keep a deprecated item for longer, for example when a removal would break a published pipeline that has no replacement yet. The reason MUST be recorded in the Deprecated section of the release that keeps it, together with the release that will remove it.

Do not remove a deprecated name in the same pull request that deprecates it. The release manager checks before each release that everything deprecated in the previous release is removed (see [Releasing](https://leepusp.github.io/rotifer/contributing/releasing.html)).

## Versioning

ROTIFER follows [semantic versioning](https://semver.org) (MAJOR.MINOR.PATCH). The NumPy project's [backwards compatibility and deprecation policy](https://numpy.org/neps/nep-0023-backwards-compatibility.html) is the model for the cycle above.

- MAJOR releases make incompatible changes.
- MINOR releases add features, introduce deprecations and remove items deprecated in the previous release.
- PATCH releases fix bugs and nothing else.

ROTIFER is below 1.0, where Semantic Versioning allows incompatible changes in MINOR releases. How versions work from 1.0 on will be decided in an ADR (see [Architecture decision records](https://leepusp.github.io/rotifer/contributing/adr.html)).

A bug fix that changes the numbers a function returns can still break someone's analysis. State it in the pull request and in the release notes, even though it ships in a PATCH release.
