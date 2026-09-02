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

## How the Read the Docs build works

`.readthedocs.yaml` at the repository root drives the hosted
build: Ubuntu 24.04, Python 3.12, `pip install -r
docs/requirements.txt` plus `pip install .`, then Sphinx with
`docs/conf.py`. `fail_on_warning: true` makes the hosted build as
strict as `make strict`, and htmlzip and PDF downloads are built
alongside the HTML. Nothing else is special: if `make strict`
passes locally, Read the Docs will pass.

## Troubleshooting

```{list-table}
:header-rows: 1

* - Symptom
  - Cause and fix
* - `TabError` or `SyntaxError` while autodoc imports a module
  - The source file cannot be compiled (mixed tabs and spaces, or
    a genuine syntax error). Do not edit the code as part of a
    docs change: remove the module from `docs/api/index.rst`,
    record it in `docs/OPEN_QUESTIONS.md`, and let the code owner
    fix the file.
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
