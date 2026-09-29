# Building the docs

TL;DR, I want to build the docs:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
pip install -e .
cd docs
make livehtml
```

## Environment setup

The docs build needs Python 3.12 or newer, the pinned documentation
dependencies, and the rotifer package itself:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
pip install -e .
```

The build is fully offline. Intersphinx inventories are vendored in
`docs/_intersphinx/`, and the DataTables and iconify assets are
vendored in `docs/_static/vendor/`, so no network access, account or
API credential is needed.

## Makefile targets

All targets run from the `docs/` directory.

```{list-table}
:header-rows: 1

* - Target
  - What it does
* - `make html`
  - Builds the HTML site into `_build/html`.
* - `make livehtml`
  - Serves the site with sphinx-autobuild and rebuilds on every
    file save. Open the printed URL and edit; the browser reloads
    itself.
* - `make strict`
  - The CI gate: `-W` turns warnings into errors, `--keep-going`
    reports all of them instead of stopping at the first, `-n`
    (nitpicky) fails on broken cross-references.
* - `make doctest`
  - Runs every docstring and page example that is not marked
    `# doctest: +SKIP`.
* - `make linkcheck`
  - Checks external links. The only target that uses the network;
    expect occasional false failures from rate limiting servers.
* - `make clean`
  - Removes `_build`.
```

Before opening a pull request, run `make strict` and
`make doctest`. Both must exit cleanly.

## Live preview

```bash
cd docs
make livehtml
```

sphinx-autobuild watches the source tree, rebuilds changed pages
and reloads the browser. Docstring edits in `lib/rotifer` are
picked up on the next rebuild of the page that renders them; when
in doubt, touch the corresponding file under `docs/api` or restart
the server.

## How the hosted build works

The site is published to GitHub Pages at
<https://leepusp.github.io/rotifer/> by
`.github/workflows/docs.yml`. On every push to `docs/v1` that
touches `docs/`, `lib/rotifer/` or `setup.py`, the workflow runs
`sphinx-build` on Ubuntu with Python 3.12 and uploads
`docs/_build/html` as the Pages artifact.

It installs `docs/requirements.txt` and then `pip install -e .` —
and only those. Never add the repository-root `requirements.txt` to
the workflow: it lists conda and system tools (`blast`, `datamash`,
`famsa`) that do not exist on PyPI, plus `pygraphviz`, which needs
Graphviz headers. The editable install is what makes numpy and
pandas available to autodoc, which imports them for real, and what
lets `conf.py` read the version from the package metadata; without
it the published title would read `0.0.0`.

The hosted build does **not** pass `-W`, so unlike `make strict` a
warning will not fail the deploy. Run `make strict` locally before
opening a pull request; the hosted build is not a substitute for it.

## Troubleshooting

```{list-table}
:header-rows: 1

* - Symptom
  - Cause and fix
* - `TabError` or `SyntaxError` while autodoc imports a module
  - The source file cannot be compiled (mixed tabs and spaces, or
    a genuine syntax error). Do not edit the code as part of a
    docs change: remove the module from `docs/api/index.rst`,
    open an issue on the
    [tracker](https://github.com/leepusp/rotifer/issues), and let
    the code owner fix the file.
* - `ModuleNotFoundError` for a scientific dependency during the
    build
  - autodoc imported a module whose dependency is not installed.
    Add the missing top-level package to `MOCKED_IMPORTS` in the
    tunables block of `conf.py`. Never pip-install heavy
    scientific packages just to build docs.
* - Duplicate object description warnings
  - The same object is documented twice, usually by an explicit
    `automodule`/`autoclass` directive on a hand-written page that
    autosummary already generates. Keep exactly one documenting
    location per object: hand-written pages link to the generated
    pages, they do not re-document objects.
* - Broken cross-reference under nitpicky mode
    (`reference target not found`)
  - The role points at something that is not documented. Fix the
    reference, or, if the target genuinely cannot exist (private
    stdlib classes, mocked externals), add it to
    `NITPICK_IGNORES` in `conf.py` with a comment saying why.
* - A page renders a docstring as one unformatted blob
  - The docstring is Google style or free form. Rewrite it to the
    NumPy standard; do not re-add napoleon.
```
