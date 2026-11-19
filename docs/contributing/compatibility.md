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
2. List the deprecation in the changelog or release notes.
3. Remove it in a later release, after the period below.

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

A deprecated name stays for at least one release after the release that introduced the warning. If releases are irregular, it stays for at least 6 months. A maintainer MAY extend that period when known users need more time.

Do not remove a deprecated name in the same pull request that deprecates it.

## Versioning

ROTIFER follows [semantic versioning](https://semver.org) (MAJOR.MINOR.PATCH). The NumPy project's [backwards compatibility and deprecation policy](https://numpy.org/neps/nep-0023-backwards-compatibility.html) is the model for the cycle above.

- MAJOR releases remove deprecated API or make incompatible changes.
- MINOR releases add features and introduce deprecations.
- PATCH releases fix bugs and nothing else.

A bug fix that changes the numbers a function returns can still break someone's analysis. State it in the pull request and in the release notes, even though it ships in a PATCH release.
