"""Give the landing page its copy and its background sequences.

The landing page (docs/index.md, layout "rotifer-landing") is rendered by
templates that contain markup only. This extension loads the two data files
those templates read and hands them over as template variables, for the
landing page only:

- ``landing``: the copy, from docs/_data/landing.<language>.yml.
- ``drift_rows``: rows of protein sequence for the hero background, built
  from docs/_static/landing/sequences.json.

Configuration (in conf.py):

- ``landing_page``: the page name that gets the context (default "index").
- ``landing_rows``: how many background rows to build (default 14).
- ``landing_row_length``: minimum characters in one row before it is
  doubled for the seamless loop (default 480).
"""

import json
from pathlib import Path

import yaml

DATA_DIR = Path(__file__).resolve().parents[1] / "_data"
SEQUENCES = Path(__file__).resolve().parents[1] / "_static" / "landing" / "sequences.json"


def load_copy(language):
    """The landing copy for ``language``, falling back to English."""
    for lang in (language, "en"):
        path = DATA_DIR / f"landing.{lang}.yml"
        if lang and path.exists():
            return yaml.safe_load(path.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"no landing copy in {DATA_DIR}")


def build_rows(records, count, min_length):
    """Rows of sequence text, cycling through the records.

    Each row repeats one sequence (separated by spaces) until it reaches
    ``min_length``. The template prints the row twice so a CSS translation
    of -50% loops without a seam. Speed and direction alternate so
    neighbouring rows never move in lockstep.
    """
    rows = []
    for i in range(count):
        rec = records[i % len(records)]
        seq = rec["sequence"]
        text = seq
        while len(text) < min_length:
            text += "   " + seq
        rows.append(
            {
                "text": text,
                "entry": rec["entry"],
                "reverse": i % 2 == 1,
                # 140s to 230s per loop: slow enough to read as texture.
                "duration": 140 + (i * 37) % 91,
            }
        )
    return rows


def on_page_context(app, pagename, templatename, context, doctree):
    if pagename != app.config.landing_page:
        return
    context["landing"] = load_copy(app.config.language)
    records = json.loads(SEQUENCES.read_text(encoding="utf-8"))
    context["drift_rows"] = build_rows(
        records, app.config.landing_rows, app.config.landing_row_length
    )


def setup(app):
    app.add_config_value("landing_page", "index", "html")
    app.add_config_value("landing_rows", 14, "html")
    app.add_config_value("landing_row_length", 480, "html")
    app.connect("html-page-context", on_page_context)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
