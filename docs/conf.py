# Configuration file for the Sphinx documentation builder.
#
# Everything a maintainer may want to change lives in the TUNABLES block
# below. Nothing tunable appears further down this file.

import os
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as pkg_version

# ---------------------------------------------------------------------------
# --- TUNABLES ---
# ---------------------------------------------------------------------------

# Project metadata. The version is read from the installed package metadata
# and is never hardcoded here.
PROJECT_NAME = "ROTIFER"
PROJECT_AUTHOR = "Robson F. de Souza and contributors"
GITHUB_URL = "https://github.com/leepbioinfo/rotifer"

# Shibuya accent family (radix color name). This only controls surfaces the
# theme owns; the exact brand tokens live in _static/theme.css.
ACCENT_COLOR = "teal"

# Third-party packages that autodoc must NOT import for real. These are the
# heavy scientific dependencies of the rotifer package. Extend this list if
# a new dependency breaks the documentation build with ModuleNotFoundError.
MOCKED_IMPORTS = [
    "Bio",
    "BCBio",
    "ete3",
    "ete4",
    "tqdm",
    "requests",
    "psycopg2",
    "matplotlib",
    "seaborn",
    "networkx",
    "sqlalchemy",
    "yaml",
    "scipy",
    "IPython",
    "igraph",
    "leidenalg",
    "pyhmmer",
    "pynvml",
    "psutil",
    "dash",
    "dash_bio",
    "dash_cytoscape",
    "reportlab",
    "weasyprint",
    "pdfplumber",
    "PyPDF2",
    "gspread",
    "oauth2client",
    "argcomplete",
    "ascii_graph",
    "clickhouse_driver",
    "community",
    "joblib",
    "pygraphviz",
    "termcolor",
    "textdistance",
    "unipressed",
]

# Intersphinx targets. Inventories are vendored in docs/_intersphinx so the
# build never touches the network. To refresh one, download the objects.inv
# from the URL below and overwrite the local file.
INTERSPHINX_TARGETS = {
    "python": ("https://docs.python.org/3/", "_intersphinx/python.inv"),
    "numpy": ("https://numpy.org/doc/stable/", "_intersphinx/numpy.inv"),
    "pandas": ("https://pandas.pydata.org/docs/", "_intersphinx/pandas.inv"),
    "biopython": ("https://biopython.org/docs/latest/", "_intersphinx/biopython.inv"),
}

# References that genuinely cannot be resolved, one per line, each with the
# reason it is ignored.
NITPICK_IGNORES = [
    ("py:class", "Bio.SeqRecord.SeqRecord"),  # Biopython inventory does not index SeqRecord as py:class
    ("py:class", "tempfile._TemporaryFileWrapper"),  # private stdlib class, not in the Python inventory
    ("py:class", "argparse._AppendAction"),  # private stdlib class exposed by rotifer.core.cli inheritance
]

# Module members that autosummary must skip because autodoc cannot process
# them (documenting them crashes or warns without a code fix).
AUTOSUMMARY_SKIP_MEMBERS = [
    "from_genbank",  # module-level @classmethod in rotifer.genome, not a callable object
]

# Legacy docstrings contain malformed reST that cannot be fixed without
# editing code, which is outside the scope of a docs change. Suppressing the
# docutils warning category keeps the strict build clean; the affected
# docstrings are listed in docs/OPEN_QUESTIONS.md. Set to False while
# writing new docstrings to lint them.
SUPPRESS_LEGACY_DOCSTRING_WARNINGS = True

# Legacy docstrings also contain bare ">>>" example blocks that are not
# runnable (some are not even valid Python). Only examples inside explicit
# ".. doctest::" directives are executed by "make doctest" until the
# docstrings are rewritten to the NumPy standard.
DOCTEST_ONLY_EXPLICIT_BLOCKS = True

# Feature toggles.
ENABLE_DATATABLES = True  # sortable/searchable tables via sphinx-datatables
ENABLE_ICONIFY = True  # icons via sphinx-iconify (assets self-hosted in _static/vendor)
COPYBUTTON_STRIP_PROMPTS = True  # strip ">>> " and "$ " prompts when copying code

# ---------------------------------------------------------------------------
# --- END OF TUNABLES ---
# ---------------------------------------------------------------------------

# Make the package importable when it is not installed (local builds run
# against the checkout; Read the Docs installs the package with pip).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))

project = PROJECT_NAME
author = PROJECT_AUTHOR
try:
    release = pkg_version("rotifer")
except PackageNotFoundError:
    release = "0.0.0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "numpydoc",
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
]
if ENABLE_ICONIFY:
    extensions.append("sphinx_iconify")
if ENABLE_DATATABLES:
    # sphinxcontrib.jquery ships with sphinx-datatables but must be loaded
    # explicitly: sphinx-datatables 1.0.0 sets it up during page rendering,
    # which is too late for jQuery to be installed into the page.
    extensions.append("sphinxcontrib.jquery")
    extensions.append("sphinx_datatables")

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "OPEN_QUESTIONS.md",  # working notes, not part of the site
    "STEP-BY-STEP.md",  # superseded scratch notes kept for history
    "_intersphinx",
]

# -- Markdown (MyST) ---------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

# -- Autodoc and autosummary -------------------------------------------------

autodoc_mock_imports = MOCKED_IMPORTS
autodoc_member_order = "groupwise"
autosummary_generate = True
autosummary_context = {
    "skip_members": AUTOSUMMARY_SKIP_MEMBERS,
}

if SUPPRESS_LEGACY_DOCSTRING_WARNINGS:
    suppress_warnings = ["docutils"]

# -- numpydoc ----------------------------------------------------------------
# numpydoc is the only docstring parser. Do not add sphinx.ext.napoleon:
# running both produces duplicate Parameters sections.

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

# -- Cross-references --------------------------------------------------------

intersphinx_mapping = INTERSPHINX_TARGETS
nitpicky = True
nitpick_ignore = NITPICK_IGNORES

extlinks = {
    "ghsrc": (GITHUB_URL + "/blob/master/%s", "%s"),
}

# -- Doctests ----------------------------------------------------------------
# Examples run against the real numpy/pandas but the heavy scientific
# dependencies are replaced with lightweight stubs, mirroring
# autodoc_mock_imports. Examples that reach a live database endpoint must be
# marked with ``# doctest: +SKIP``.

if DOCTEST_ONLY_EXPLICIT_BLOCKS:
    doctest_test_doctest_blocks = ""

doctest_global_setup = f"""
import sys
from unittest.mock import MagicMock
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec

_MOCKS = {MOCKED_IMPORTS!r}

class _StubLoader(Loader):
    def create_module(self, spec):
        m = MagicMock()
        m.__name__ = spec.name
        m.__path__ = []
        m.__version__ = "0.0"
        return m
    def exec_module(self, module):
        pass

class _StubFinder(MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in _MOCKS:
            return ModuleSpec(name, _StubLoader(), is_package=True)
        return None

if not any(isinstance(f, _StubFinder) for f in sys.meta_path):
    sys.meta_path.insert(0, _StubFinder())
"""

# -- Copy button -------------------------------------------------------------

if COPYBUTTON_STRIP_PROMPTS:
    copybutton_prompt_text = r">>> |\.\.\. |\$ "
    copybutton_prompt_is_regexp = True

# -- DataTables (self-hosted, no CDN) ----------------------------------------

datatables_js = "vendor/datatables.min.js"
datatables_css = "vendor/datatables.min.css"
datatables_options = {
    "paging": False,
    "info": False,
}

# -- Iconify (self-hosted, no CDN) -------------------------------------------

iconify_script_url = "vendor/iconify-icon.min.js"

# -- HTML output -------------------------------------------------------------

html_theme = "shibuya"
html_static_path = ["_static"]
html_css_files = ["theme.css"]
html_title = f"{PROJECT_NAME} {release}"

html_theme_options = {
    "accent_color": ACCENT_COLOR,
    "github_url": GITHUB_URL,
    "show_ai_links": False,  # keep the site free of external services
    "nav_links": [
        {"title": "Data access", "url": "db/index"},
        {"title": "API", "url": "api/index"},
        {"title": "Contributing", "url": "contributing/index"},
    ],
}

# Drop the theme sidebars that call external services (repo stats, ads).
html_sidebars = {
    "**": ["sidebars/localtoc.html"],
}


def setup(app):
    # The icon preload must be registered before the iconify web component
    # script so that icons render without contacting the iconify API.
    if ENABLE_ICONIFY:
        app.add_js_file("vendor/iconify-preload.js", priority=400)
