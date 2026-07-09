"""
operon_fig.py
=============

Draw gene-neighborhood ("operon") figures from a long-format table of
genes/proteins that have already been grouped into genomic blocks (e.g.
one block per hit of interest plus its surrounding genomic context).

Each block is drawn as one row of the figure:

    [ text label ]  [gene] [gene] [QUERY] [gene] [gene] ...

The text label is a 3-line box with the reference query protein id, the
block id and the organism name. Every neighbor gene is drawn pointing
right, regardless of its real strand, so a row never has arrows pointing
in different directions; only the query gene's own arrow still reflects
its real strand. Genes are colored by a domain/annotation label (Pfam,
Aravind, or free-text product, depending on `label_col`).

If a block contains more than one query gene, the one closest to the
middle of the block (by gene order) is used as that block's "reference"
query -- it is what orientation-normalization, query-centering and the
row label are anchored to. Every actual query gene is still outlined in
red; only the anchor choice is affected.

Required input columns
-----------------------
pid     : protein id (used in the row label and to find the query gene)
strand  : +1 / -1 (controls arrow direction and orientation normalization)
query   : 1/'1'/True marks the query gene of a block; everything else is
          treated as a non-query gene. Optional -- if absent, no gene is
          treated as a query.
... plus whatever `group_col`, `org_col` and `label_col` point at.

`genome_overview_fig` additionally uses `nucleotide`, `start`, `end`
(per-gene genomic coordinates) and, optionally, `nlen` (total contig
length, for scaling); see its docstring.

Public entry points
--------------------
neighborhood_figure(df, ...) -> pandas.DataFrame
    Builds the per-block neighborhood figure (one row per block, gene
    detail) and writes it to `output_file`.
genome_overview_fig(df, ...) -> pandas.DataFrame
    Builds a genome-wide companion figure: one horizontal track per
    contig, with every block marked at its actual genomic position --
    "where are these neighborhoods", as opposed to neighborhood_figure's
    "what's in each neighborhood".
build_html_report(df, ...) -> str
    Runs both of the above and assembles a single, self-contained HTML
    page with the genome overview, the neighborhood figure, and the
    original input table (search box included), for sharing/viewing
    everything at once.

Everything else below is a small, independently testable helper
function. None of them are defined *inside* another function (the
original had a `pad_and_escape` closure) -- they all live at module
level so they can be imported and unit tested on their own, and so a
future change to one of them does not require re-reading the whole
pipeline.

    resolve_domain_labels        fill the 'domain' column (labels/colors)
    flag_query_rows              add the boolean 'is_query' column
    rename_label_values          apply a user rename dict to the raw
                                  label column before domain resolution
    prepare_dataframe            run the helpers above + housekeeping
    compute_label_width          shared padding width for row labels
    pad_and_escape               pad + HTML-escape one label string
    build_row_label_html         the 3-line HTML label for one block
    build_color_map              domain -> fill color, incl. user overrides
    select_reference_query_index which query anchors a multi-query block
    normalize_block_strand       optionally mirror a block to a common
                                  query orientation
    gene_node_style               Graphviz node attributes for one gene
    add_block_to_graph            add one full row (label + genes) to the
                                  graph, with optional left/right padding
    chain_align_nodes             pull one node per row into the same
                                  visual column via high-weight invisible
                                  edges
    neighborhood_figure            main neighborhood-figure orchestrator
    compute_block_extents         one row per block: contig, span, ref query
    assign_label_lanes            stagger overlapping labels into lanes
    build_genome_overview_svg     static genome-wide SVG from extents
    build_genome_overview_interactive_html  zoomable overview (used in report)
    genome_overview_fig           genome-wide-figure orchestrator
    render_dataframe_html         a dataframe as a plain <table>
    render_neighborhood_svgs_by_block  one SVG per block (per-result views)
    build_gene_tooltip_html       per-protein hover "info window" body
    annotate_neighborhood_svg     inject those tooltips into a graphviz SVG
    build_tabs_and_panels_html    (removed; replaced by build_neighborhood_panels)
    build_neighborhood_panels     pop-up selector + per-block figure/table panels
    render_table_card             sortable/filterable/downloadable table widget
    compute_domain_stats          reference-query + full-architecture domain counts
    build_bar_chart_svg           horizontal bar chart SVG from a counts table
    build_stats_section_html      toggleable statistics section HTML
    build_html_report             combine everything into one HTML page

How the drawing actually works (no GUI toolkit involved)
--------------------------------------------------------
There is no Qt, matplotlib, or other plotting/GUI library here. Two
very different rendering paths are used:

  * The neighborhood figures (the gene-arrow rows) are laid out by
    **Graphviz** -- the same C graph-layout engine behind `dot` -- which
    we drive from Python through the **pygraphviz** binding. Each gene
    is a Graphviz node (`shape=rarrow`/`larrow`/`triangle`), each row is
    a same-rank subgraph, and invisible weighted edges nudge things into
    alignment; Graphviz's `dot` engine does the placement and writes
    **SVG**. We then do light text surgery on that SVG to attach
    per-protein hover windows (`annotate_neighborhood_svg`).
  * The genome overview and the whole report page are just **hand-written
    SVG/HTML/CSS/JavaScript** strings -- no library. Zoom/pan and the
    hover tooltips are a small vanilla-JS block embedded in the page.

Supporting libraries are **pandas**/**numpy** (data wrangling) and
**seaborn** (only to pick color palettes). The output is a static,
self-contained `.html` file that runs in any browser; nothing runs
server-side or in a desktop toolkit.
"""

import html
import os
import re
import shutil
import tempfile
from string import Template

import numpy as np
import pandas as pd
import pygraphviz as pgv
import seaborn as sns

# Domain/annotation values that never get an automatic color, because
# they are generic/structural rather than informative (signal peptides,
# transmembrane regions, etc.) or simply mean "no annotation".
DEFAULT_IGNORE_DOMAINS = ['TM', 'SP', 'LP', 'LIPO', 'SIG']

# Default branding logos, pre-prepared for inline embedding.
# Pass to `build_html_report` as `header_logo` / `footer_logo`, or use
# `read_svg_logo(path)` to swap in your own SVG file.
SHARP_HEADER_LOGO = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" zoomAndPan="magnify" viewBox="0 0 810 1012.49997" preserveAspectRatio="xMidYMid meet" version="1.0"><defs><filter x="0%" y="0%" width="100%" height="100%" id="050de061e5"><feColorMatrix values="0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0" color-interpolation-filters="sRGB"/></filter><filter x="0%" y="0%" width="100%" height="100%" id="4fcb26cbf3"><feColorMatrix values="0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0.2126 0.7152 0.0722 0 0" color-interpolation-filters="sRGB"/></filter><g/><clipPath id="2b33e3c2d3"><path d="M 48.925781 145.703125 L 723.054688 145.703125 L 723.054688 819.835938 L 48.925781 819.835938 Z M 48.925781 145.703125 " clip-rule="nonzero"/></clipPath><clipPath id="0475d7ad1b"><path d="M 385.992188 145.703125 C 199.835938 145.703125 48.925781 296.613281 48.925781 482.769531 C 48.925781 668.925781 199.835938 819.835938 385.992188 819.835938 C 572.148438 819.835938 723.054688 668.925781 723.054688 482.769531 C 723.054688 296.613281 572.148438 145.703125 385.992188 145.703125 Z M 385.992188 145.703125 " clip-rule="nonzero"/></clipPath><clipPath id="d20bebb775"><path d="M 0.925781 0.703125 L 675.054688 0.703125 L 675.054688 674.835938 L 0.925781 674.835938 Z M 0.925781 0.703125 " clip-rule="nonzero"/></clipPath><clipPath id="4f3de7ccd1"><path d="M 337.992188 0.703125 C 151.835938 0.703125 0.925781 151.613281 0.925781 337.769531 C 0.925781 523.925781 151.835938 674.835938 337.992188 674.835938 C 524.148438 674.835938 675.054688 523.925781 675.054688 337.769531 C 675.054688 151.613281 524.148438 0.703125 337.992188 0.703125 Z M 337.992188 0.703125 " clip-rule="nonzero"/></clipPath><clipPath id="9c0b549aea"><rect x="0" width="676" y="0" height="675"/></clipPath><clipPath id="dbcc122e99"><path d="M 0.199219 305 L 809.800781 305 L 809.800781 679.808594 L 0.199219 679.808594 Z M 0.199219 305 " clip-rule="nonzero"/></clipPath><mask id="c9d8498719"><g filter="url(#050de061e5)"><g filter="url(#4fcb26cbf3)" transform="matrix(1.189189, 0, 0, 1.189889, -153.303221, 304.994196)"><image x="0" y="0" width="851" xlink:href="data:image/png;base64,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" height="315" preserveAspectRatio="xMidYMid meet"/></g></g></mask><clipPath id="a5ce00304d"><rect x="0" width="256" y="0" height="51"/></clipPath></defs><g clip-path="url(#2b33e3c2d3)"><g clip-path="url(#0475d7ad1b)"><g transform="matrix(1, 0, 0, 1, 48, 145)"><g clip-path="url(#9c0b549aea)"><g clip-path="url(#d20bebb775)"><g clip-path="url(#4f3de7ccd1)"><path fill="#8f003c" d="M 0.925781 0.703125 L 675.054688 0.703125 L 675.054688 674.835938 L 0.925781 674.835938 Z M 0.925781 0.703125 " fill-opacity="1" fill-rule="nonzero"/></g></g></g></g></g></g><g clip-path="url(#dbcc122e99)"><g mask="url(#c9d8498719)"><g transform="matrix(1.189189, 0, 0, 1.189889, -153.303221, 304.994196)"><image x="0" y="0" width="851" xlink:href="data:image/png;base64,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" height="315" preserveAspectRatio="xMidYMid meet"/></g></g></g><g transform="matrix(1, 0, 0, 1, 257, 94)"><g clip-path="url(#a5ce00304d)"><g fill="#ffffff" fill-opacity="1"><g transform="translate(7.721789, 44.304424)"><g><path d="M 10.453125 -21.0625 L 15.921875 -5.671875 L 13.015625 -4.625 L 2.015625 -11.25 L 5.34375 -1.90625 L 1.828125 -0.65625 L -3.640625 -16.046875 L -0.703125 -17.09375 L 10.265625 -10.46875 L 6.953125 -19.8125 Z M 10.453125 -21.0625 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(30.352066, 36.494205)"><g><path d="M 13.203125 -6.875 L 14.03125 -3.953125 L 1.859375 -0.53125 L -2.5625 -16.25 L 9.328125 -19.59375 L 10.15625 -16.671875 L 1.875 -14.34375 L 2.828125 -10.9375 L 10.140625 -12.984375 L 10.9375 -10.15625 L 3.625 -8.109375 L 4.640625 -4.46875 Z M 13.203125 -6.875 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(50.238438, 30.868056)"><g><path d="M 22.96875 -21.15625 L 20.875 -4.109375 L 16.890625 -3.328125 L 11.234375 -13.5 L 9.734375 -1.921875 L 5.765625 -1.140625 L -2.65625 -16.109375 L 1.203125 -16.875 L 7.03125 -6.328125 L 8.59375 -18.328125 L 12.015625 -19 L 17.921875 -8.390625 L 19.40625 -20.453125 Z M 22.96875 -21.15625 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(81.945524, 25.143312)"><g/></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(93.666381, 23.632309)"><g><path d="M 7.71875 -16.96875 C 9.164062 -17.082031 10.4375 -16.9375 11.53125 -16.53125 C 12.632812 -16.132812 13.503906 -15.515625 14.140625 -14.671875 C 14.785156 -13.835938 15.15625 -12.816406 15.25 -11.609375 C 15.34375 -10.410156 15.132812 -9.34375 14.625 -8.40625 C 14.125 -7.46875 13.359375 -6.722656 12.328125 -6.171875 C 11.304688 -5.617188 10.070312 -5.285156 8.625 -5.171875 L 5.359375 -4.921875 L 5.703125 -0.4375 L 1.9375 -0.140625 L 0.6875 -16.421875 Z M 8.1875 -8.234375 C 9.320312 -8.316406 10.164062 -8.625 10.71875 -9.15625 C 11.269531 -9.6875 11.507812 -10.40625 11.4375 -11.3125 C 11.363281 -12.226562 11.015625 -12.90625 10.390625 -13.34375 C 9.765625 -13.789062 8.882812 -13.972656 7.75 -13.890625 L 4.6875 -13.65625 L 5.125 -8 Z M 8.1875 -8.234375 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(115.785637, 22.106735)"><g><path d="M 12.34375 -0.109375 L 9.140625 -4.640625 L 8.953125 -4.625 L 5.671875 -4.59375 L 5.71875 -0.046875 L 1.9375 -0.015625 L 1.78125 -16.34375 L 8.84375 -16.40625 C 10.289062 -16.425781 11.550781 -16.203125 12.625 -15.734375 C 13.695312 -15.265625 14.523438 -14.585938 15.109375 -13.703125 C 15.691406 -12.816406 15.988281 -11.769531 16 -10.5625 C 16.007812 -9.34375 15.722656 -8.289062 15.140625 -7.40625 C 14.566406 -6.53125 13.75 -5.851562 12.6875 -5.375 L 16.40625 -0.15625 Z M 12.171875 -10.515625 C 12.160156 -11.429688 11.859375 -12.128906 11.265625 -12.609375 C 10.671875 -13.097656 9.804688 -13.335938 8.671875 -13.328125 L 5.59375 -13.296875 L 5.640625 -7.609375 L 8.71875 -7.640625 C 9.851562 -7.648438 10.710938 -7.90625 11.296875 -8.40625 C 11.890625 -8.90625 12.179688 -9.609375 12.171875 -10.515625 Z M 12.171875 -10.515625 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(138.000846, 22.03285)"><g><path d="M 9.8125 0.890625 C 8.113281 0.785156 6.609375 0.324219 5.296875 -0.484375 C 3.992188 -1.296875 2.992188 -2.363281 2.296875 -3.6875 C 1.609375 -5.019531 1.3125 -6.488281 1.40625 -8.09375 C 1.507812 -9.6875 1.984375 -11.097656 2.828125 -12.328125 C 3.679688 -13.566406 4.804688 -14.507812 6.203125 -15.15625 C 7.597656 -15.800781 9.144531 -16.070312 10.84375 -15.96875 C 12.53125 -15.863281 14.03125 -15.40625 15.34375 -14.59375 C 16.65625 -13.78125 17.660156 -12.707031 18.359375 -11.375 C 19.054688 -10.050781 19.351562 -8.59375 19.25 -7 C 19.15625 -5.394531 18.679688 -3.972656 17.828125 -2.734375 C 16.972656 -1.503906 15.84375 -0.566406 14.4375 0.078125 C 13.039062 0.722656 11.5 0.992188 9.8125 0.890625 Z M 10.015625 -2.328125 C 10.984375 -2.273438 11.863281 -2.441406 12.65625 -2.828125 C 13.457031 -3.210938 14.101562 -3.789062 14.59375 -4.5625 C 15.082031 -5.332031 15.359375 -6.222656 15.421875 -7.234375 C 15.484375 -8.242188 15.316406 -9.15625 14.921875 -9.96875 C 14.523438 -10.789062 13.957031 -11.445312 13.21875 -11.9375 C 12.476562 -12.425781 11.625 -12.695312 10.65625 -12.75 C 9.6875 -12.8125 8.796875 -12.644531 7.984375 -12.25 C 7.179688 -11.863281 6.535156 -11.285156 6.046875 -10.515625 C 5.554688 -9.753906 5.28125 -8.867188 5.21875 -7.859375 C 5.15625 -6.847656 5.320312 -5.929688 5.71875 -5.109375 C 6.113281 -4.285156 6.6875 -3.628906 7.4375 -3.140625 C 8.1875 -2.660156 9.046875 -2.390625 10.015625 -2.328125 Z M 10.015625 -2.328125 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(162.736992, 23.835322)"><g><path d="M 4.875 0.90625 C 3.78125 0.757812 2.800781 0.425781 1.9375 -0.09375 C 1.070312 -0.625 0.390625 -1.300781 -0.109375 -2.125 L 2.28125 -4.375 C 3.019531 -3.070312 3.957031 -2.351562 5.09375 -2.21875 C 6.601562 -2.019531 7.46875 -2.8125 7.6875 -4.59375 L 8.671875 -12.296875 L 2.984375 -13.015625 L 3.359375 -16.03125 L 12.796875 -14.84375 L 11.46875 -4.328125 C 11.21875 -2.347656 10.535156 -0.921875 9.421875 -0.046875 C 8.304688 0.828125 6.789062 1.144531 4.875 0.90625 Z M 4.875 0.90625 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(180.365157, 26.17341)"><g><path d="M 14.890625 -0.3125 L 14.328125 2.671875 L 1.90625 0.359375 L 4.890625 -15.6875 L 17.03125 -13.4375 L 16.484375 -10.46875 L 8.03125 -12.03125 L 7.375 -8.53125 L 14.84375 -7.140625 L 14.3125 -4.265625 L 6.84375 -5.65625 L 6.15625 -1.9375 Z M 14.890625 -0.3125 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(200.75288, 30.117753)"><g><path d="M 9.390625 2.671875 C 7.773438 2.253906 6.40625 1.53125 5.28125 0.5 C 4.15625 -0.519531 3.375 -1.75 2.9375 -3.1875 C 2.5 -4.625 2.484375 -6.128906 2.890625 -7.703125 C 3.285156 -9.265625 4.015625 -10.566406 5.078125 -11.609375 C 6.148438 -12.660156 7.425781 -13.375 8.90625 -13.75 C 10.382812 -14.125 11.941406 -14.101562 13.578125 -13.6875 C 14.929688 -13.332031 16.101562 -12.769531 17.09375 -12 C 18.082031 -11.238281 18.835938 -10.316406 19.359375 -9.234375 L 16.453125 -7.65625 C 15.691406 -9.164062 14.523438 -10.125 12.953125 -10.53125 C 11.960938 -10.78125 11.03125 -10.785156 10.15625 -10.546875 C 9.28125 -10.304688 8.53125 -9.859375 7.90625 -9.203125 C 7.289062 -8.546875 6.859375 -7.726562 6.609375 -6.75 C 6.359375 -5.769531 6.347656 -4.835938 6.578125 -3.953125 C 6.804688 -3.078125 7.25 -2.328125 7.90625 -1.703125 C 8.5625 -1.078125 9.382812 -0.640625 10.375 -0.390625 C 11.945312 0.0117188 13.429688 -0.269531 14.828125 -1.25 L 16.640625 1.5 C 15.648438 2.226562 14.535156 2.691406 13.296875 2.890625 C 12.066406 3.085938 10.765625 3.015625 9.390625 2.671875 Z M 9.390625 2.671875 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(222.238595, 35.794573)"><g><path d="M 9.140625 -10.984375 L 4.171875 -12.578125 L 5.109375 -15.5 L 18.65625 -11.125 L 17.71875 -8.203125 L 12.734375 -9.8125 L 8.65625 2.796875 L 5.0625 1.625 Z M 9.140625 -10.984375 "/></g></g></g><g fill="#ffffff" fill-opacity="1"><g transform="translate(240.874326, 42.039073)"><g/></g></g></g></g></svg>'

SHARP_FOOTER_LOGO = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" zoomAndPan="magnify" viewBox="0 0 240 239.999995" preserveAspectRatio="xMidYMid meet" version="1.0"><defs><filter x="0%" y="0%" width="100%" height="100%" id="1a64d36e95"><feColorMatrix values="0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0" color-interpolation-filters="sRGB"/></filter><filter x="0%" y="0%" width="100%" height="100%" id="999fee5746"><feColorMatrix values="0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0.2126 0.7152 0.0722 0 0" color-interpolation-filters="sRGB"/></filter><clipPath id="d3a2c2e922"><path d="M 22 2 L 213.882812 2 L 213.882812 240 L 22 240 Z M 22 2 " clip-rule="nonzero"/></clipPath><mask id="972b4a3af0"><g filter="url(#1a64d36e95)"><g filter="url(#999fee5746)" transform="matrix(0.197531, 0, 0, 0.197531, 21.883769, 1.683369)"><image x="0" y="0" width="972" xlink:href="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA8wAAAS/CAAAAAD6Ui5TAAAAAmJLR0QA/4ePzL8AACAASURBVHic7N13gFxV+T7w95x7Z3fTew8hhBRDiITemxQFBAEVEZCiggj4VQEpirSfKEVAQKUKiESRDkpVqrRAaCaUEJIQEkJ6ssm2ufec9/n9sW3u7MzuzO7sZmOezx+wU+6dM5P73tPPESEiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiIi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ERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERH9//bgQAAAAABAkL/1CAtUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAC9Xc8YgNqRfTAAAAAElFTkSuQmCC" height="1215" preserveAspectRatio="xMidYMid meet"/></g></g></mask></defs><g clip-path="url(#d3a2c2e922)"><g mask="url(#972b4a3af0)"><g transform="matrix(0.197531, 0, 0, 0.197531, 21.883769, 1.683369)"><image x="0" y="0" width="972" xlink:href="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA8wAAAS/CAIAAABQW+bYAAAABmJLR0QA/wD/AP+gvaeTAAAgAElEQVR4nOzde3Bc133g+e8593Y3XgRfAiRKfECkZFmAYslWnJh2EkNljxWX7TiZKWhqZuOJppKqeJJMZr3OpCpJbbH1z+7OzMapeayr9pEtu2az2RUzU1OJJ4kzdgZMYlvexGPZHsCS9TAkUiRFSHyAD6Af957943aTVOlhioIIAvp+CiWBYPe9597bRP/ur3/nd0CSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmSJEmS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AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACcTL8DEV+iHAermvsAAAAASUVORK5CYII=" height="1215" preserveAspectRatio="xMidYMid meet"/></g></g></g></svg>'


# ---------------------------------------------------------------------------
# Column preparation
# ---------------------------------------------------------------------------

def resolve_domain_labels(df, label_col='pfam'):
    """
    Build the 'domain' column that gene node labels/colors are based on.

    When `label_col == 'pfam'` (the default), a gene with no Pfam hit
    falls back to its 'aravind' annotation, then to its free-text
    'product' description, and finally to the literal string 'unk' if
    none of those are available either. This mirrors the common
    situation where only some genes in a neighborhood have a Pfam domain.

    For any other `label_col`, that column is used as-is (missing values
    become 'unk'); if the column does not exist at all, every gene gets
    'unk'.

    Parameters
    ----------
    df : pandas.DataFrame
    label_col : str

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with a 'domain' column added.
    """
    out = df.copy()

    if label_col == 'pfam':
        domain = out['pfam'] if 'pfam' in out.columns else pd.Series(np.nan, index=out.index)
        if 'aravind' in out.columns:
            domain = domain.fillna(out['aravind'])
        if 'product' in out.columns:
            domain = domain.fillna(out['product'])
        out['domain'] = domain.fillna('unk')
    elif label_col in out.columns:
        out['domain'] = out[label_col].fillna('unk')
    else:
        out['domain'] = 'unk'

    return out


def flag_query_rows(df):
    """
    Add a boolean 'is_query' column derived from the raw 'query' column.

    1, '1' and True are all treated as "this is the query gene"; a
    missing 'query' column means no gene in `df` is a query.
    """
    out = df.copy()
    if 'query' in out.columns:
        out['is_query'] = out['query'].isin([1, '1', True])
    else:
        out['is_query'] = False
    return out


def rename_label_values(df, label_col='pfam', rename_map=None):
    """
    Rename values in `label_col` through a user-supplied dictionary,
    before domain resolution/coloring happen downstream.

    `label_col` is "the architecture column" -- whatever column holds
    each gene's domain/annotation label. It defaults to 'pfam', but can
    be pointed at any column name (e.g. 'aravind' or 'product'); this
    function always renames whichever column that is, never a column
    hardcoded to literally be called 'pfam'.

    Values in this column are often a multi-domain "architecture"
    string with parts joined by '+' (e.g. 'HTH_1+LysR_substrate'), so
    renaming is done piece-by-piece on each '+'-separated component, not
    only on an exact whole-string match -- this lets a single dict entry
    like `{'GntR': 'MyFavoriteRegulator'}` rename that domain everywhere
    it shows up, whether alone or combined with other domains.

    Parameters
    ----------
    df : pandas.DataFrame
    label_col : str
    rename_map : dict[str, str] or None
        {old_name: new_name}. Pieces not present in the dict are left
        untouched. No-op if `rename_map` is empty/None or `label_col`
        is not a column in `df`.

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with the renamed column (or `df` itself, unchanged,
        if there was nothing to rename).

    Notes
    -----
    Downstream color/label lookups (`custom_colors`, `ignore_domains`)
    should reference the *renamed* values, since renaming happens before
    those steps run.
    """
    if not rename_map or label_col not in df.columns:
        return df

    def _rename_one(value):
        if pd.isna(value):
            return value
        parts = str(value).split('+')
        return '+'.join(rename_map.get(part, part) for part in parts)

    out = df.copy()
    out[label_col] = out[label_col].apply(_rename_one)
    return out


def prepare_dataframe(df, group_col='block_id', org_col='organism', label_col='pfam',
                       rename_map=None):
    """
    Normalize a raw input table into the columns the rest of this module
    relies on: 'ID' (block id), 'org_name', 'domain', 'is_query' and
    'pid_order' (an integer 0..n_blocks-1, in first-seen order, used to
    group rows into figure rows).

    Parameters mirror `neighborhood_figure`.

    Returns
    -------
    pandas.DataFrame
    """
    out = df.copy()
    out['ID'] = out[group_col] if group_col in out.columns else 'Unknown_Block'
    out['org_name'] = out[org_col] if org_col in out.columns else 'Unknown Organism'
    out = rename_label_values(out, label_col=label_col, rename_map=rename_map)
    out = resolve_domain_labels(out, label_col)
    out = flag_query_rows(out)
    out['pid_order'] = pd.factorize(out['ID'])[0]
    return out.reset_index(drop=True)


def compute_label_width(df):
    """
    Character width every row's 3-line label should be padded to, so the
    query-id / block-id / organism-name boxes line up across rows.
    """
    max_id_len = df['ID'].astype(str).str.len().max()
    max_org_len = df['org_name'].astype(str).str.len().max()
    query_pids = df.loc[df['is_query'], 'pid'].astype(str)
    max_q_len = query_pids.str.len().max() if not query_pids.empty else 8
    return max(max_id_len, max_org_len, max_q_len)


# ---------------------------------------------------------------------------
# Row label
# ---------------------------------------------------------------------------

def pad_and_escape(text, width):
    """
    Right-pad `text` to `width` characters, HTML-escape it, then turn the
    padding spaces into non-breaking spaces (&nbsp;).

    Graphviz HTML-like labels collapse normal spaces, so padding only
    works if it survives as &nbsp;. Doing this with a monospace font
    (see `build_row_label_html`) is what makes the label boxes line up
    across rows.
    """
    return html.escape(str(text).ljust(width)).replace(" ", "&nbsp;")


def build_row_label_html(query_pid, block_id, org_name, width, font_size):
    """
    Build the 3-line Graphviz HTML-like label for one block's left-hand
    label box: bold query protein id, block id, italic organism name.
    """
    query_str = pad_and_escape(query_pid, width)
    block_str = pad_and_escape(block_id, width)
    org_str = pad_and_escape(org_name, width)
    return (
        '<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}">'
        f'<B>{query_str}</B></FONT></TD></TR>'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}">'
        f'{block_str}</FONT></TD></TR>'
        f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas italic" POINT-SIZE="{font_size}">'
        f'{org_str}</FONT></TD></TR>'
        '</TABLE>>'
    )


# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

def build_color_map(df, max_colors=5, ignore_domains=None, custom_colors=None):
    """
    Decide which fill color each domain/annotation value gets.

    Priority order:

    1. Anything explicitly listed in `custom_colors` keeps that exact
       color. These do not count against `max_colors`.
    2. Domains seen on a query gene get an automatic color next (these
       are usually the reason the figure exists).
    3. The remaining `max_colors` slots go to the most frequent of the
       remaining domains, by gene count.
    4. Everything else -- including anything in `ignore_domains`, the
       literal value 'unk'/blank/'-'/'?', and anything containing the
       word "hypothetical" -- is left white (no color).

    Parameters
    ----------
    df : pandas.DataFrame
        Must already have 'domain' and 'is_query' columns (see
        `prepare_dataframe`).
    max_colors : int
        Number of *automatically chosen* colors, on top of any fixed via
        `custom_colors`.
    ignore_domains : list[str] or None
        Domain values that never get an automatic color (case
        insensitive). Defaults to `DEFAULT_IGNORE_DOMAINS`.
    custom_colors : dict[str, str] or None
        Explicit {domain_value: color} overrides. Keys should match the
        values that appear in the 'domain' column -- i.e. typically the
        raw Pfam identifiers when `label_col='pfam'` (the default).
        Values can be any color Graphviz/seaborn understands, e.g.
        '#ff8800' or 'tomato'.
        Example: `custom_colors={'LysR_substrate': '#ff8800', 'PrpF': 'tomato'}`

    Returns
    -------
    dict[str, str]
        Mapping from domain value to color string.
    """
    if ignore_domains is None:
        ignore_domains = DEFAULT_IGNORE_DOMAINS
    custom_colors = dict(custom_colors) if custom_colors else {}

    ignore_lower = [d.lower() for d in ignore_domains] + ['unk', ' ', '-', '?']
    domain_str = df['domain'].astype(str)
    is_ignorable = (
        domain_str.str.lower().isin(ignore_lower)
        | domain_str.str.lower().str.contains('hypothetical', na=False)
    )
    already_colored = domain_str.isin(custom_colors)

    query_domains = (
        df.loc[df['is_query'] & ~is_ignorable & ~already_colored, 'domain']
        .unique()
        .tolist()
    )

    remaining_slots = max(0, max_colors - len(query_domains))
    freq_domains = (
        df.loc[~is_ignorable & ~already_colored & ~domain_str.isin(query_domains), 'domain']
        .value_counts()
        .head(remaining_slots)
        .index
        .tolist()
    )

    auto_domains = query_domains + freq_domains
    palette = sns.color_palette('pastel', len(auto_domains)).as_hex() if auto_domains else []

    color_map = dict(custom_colors)
    color_map.update(dict(zip(auto_domains, palette)))
    return color_map


# ---------------------------------------------------------------------------
# Orientation
# ---------------------------------------------------------------------------

def select_reference_query_index(block_df):
    """
    Pick which query gene to use as a block's reference point for
    orientation, centering and the row label, when the block contains
    more than one.

    The query closest to the middle of the block, by gene order, is
    used -- it best represents "the middle of the region". Ties prefer
    the more upstream (lower position) one.

    Returns the row's label in `block_df.index` (not a bare position),
    so the same gene can still be found correctly after `block_df` is
    reversed (e.g. by `normalize_block_strand`) -- `.iloc[::-1]` keeps
    each row's original index label attached to it even as the row
    order changes.

    Parameters
    ----------
    block_df : pandas.DataFrame
        Rows for a single block; must have 'is_query'.

    Returns
    -------
    Any or None
        An index label from `block_df.index`, or None if the block has
        no query gene at all.
    """
    query_index_labels = block_df.index[block_df['is_query']]
    if len(query_index_labels) == 0:
        return None
    if len(query_index_labels) == 1:
        return query_index_labels[0]

    positions = np.flatnonzero(block_df['is_query'].to_numpy())
    middle = (len(block_df) - 1) / 2
    ranked = sorted(zip(positions, query_index_labels), key=lambda p: (abs(p[0] - middle), p[0]))
    return ranked[0][1]


def normalize_block_strand(block_df, normalize_orientation=True):
    """
    Optionally mirror a block so its reference query gene always points
    the same way (strand +1, drawn as a right-pointing arrow).

    Neighborhoods are usually pulled out with no regard for which strand
    the query happens to land on, which makes "upstream"/"downstream"
    mean different things from row to row. When `normalize_orientation`
    is True and this block's reference query (see
    `select_reference_query_index`) is on the minus strand, the whole
    block is reversed -- gene order *and* every strand sign flip -- which
    is equivalent to flipping the picture so it reads in the same
    direction as every other block.

    Parameters
    ----------
    block_df : pandas.DataFrame
        Rows for a single block, already in genomic (left-to-right) order.
    normalize_orientation : bool

    Returns
    -------
    pandas.DataFrame
        `block_df` unchanged, or a reversed/strand-flipped copy.
    """
    if not normalize_orientation:
        return block_df

    ref_idx = select_reference_query_index(block_df)
    if ref_idx is None or block_df.loc[ref_idx, 'strand'] != -1:
        return block_df

    flipped = block_df.iloc[::-1].copy()
    flipped['strand'] = -flipped['strand']
    return flipped


# Graphviz 'orientation' (degrees) that makes a `shape=triangle` node
# point left. Since every neighbor gene is now always drawn pointing
# right (see `gene_node_style`), a collapsed opposite-strand neighbor's
# triangle always points the other way -- left -- regardless of which
# strand the query itself happens to be on. (Verified empirically by
# rendering test shapes with Graphviz 2.43: a plain triangle points up
# by default, and `orientation=90` rotates it to point left.)
COLLAPSED_TRIANGLE_ORIENTATION = 90


# ---------------------------------------------------------------------------
# Per-gene node styling
# ---------------------------------------------------------------------------

def gene_node_style(row, query_canonical_strand, color_map, highlight_query=True,
                     collapse_opposite_strand=False, font_size=10):
    """
    Decide the Graphviz node attributes (shape/color/label/size) for one
    gene.

    Neighbor genes (anything that is not the query) are always drawn as
    a right-pointing arrow, regardless of their real strand -- this
    keeps every row reading in one consistent direction instead of a mix
    of arrows pointing every which way. The query gene is the one
    exception: its shape still reflects its real strand (right for +1,
    left for -1, a plain box for anything else), since the query's own
    orientation is usually exactly the thing worth seeing at a glance.
    Every gene is filled per `color_map` and labeled with its
    domain/annotation; the query additionally gets a red, thicker
    outline when `highlight_query` is True.

    When `collapse_opposite_strand` is True, neighbors on the strand
    *opposite* the query's (`query_canonical_strand`) are drawn instead
    as small, unlabeled, grey triangles pointing left (the opposite of
    the direction every other neighbor points) -- a lightweight
    "something is here, transcribed the other way" cue instead of giving
    them the same visual weight as same-strand neighbors. The query gene
    itself is never collapsed.

    Parameters
    ----------
    row : pandas.Series
        One gene row; needs 'strand', 'domain' and 'is_query'.
    query_canonical_strand : int
        The strand (1 or -1) this block's reference query gene has.
    color_map : dict[str, str]
    highlight_query : bool
    collapse_opposite_strand : bool
    font_size : int

    Returns
    -------
    dict
        Keyword arguments for `AGraph.add_node`.
    """
    strand_val = row.get('strand', 1)
    is_target = bool(row['is_query'])

    opposite_strand = (
        collapse_opposite_strand
        and not is_target
        and strand_val != query_canonical_strand
    )

    if opposite_strand:
        return dict(
            label='',
            shape='triangle',
            orientation=COLLAPSED_TRIANGLE_ORIENTATION,
            style='filled',
            fixedsize='true',
            width='0.18',
            height='0.18',
            fillcolor='#cccccc',
            color='#888888',
            penwidth='1',
        )

    if is_target:
        node_shape = 'rarrow' if strand_val == 1 else 'larrow' if strand_val == -1 else 'box'
    else:
        node_shape = 'rarrow'  # neighbors always point right, regardless of real strand

    return dict(
        label=str(row['domain']),
        shape=node_shape,
        style='filled',
        fixedsize='false',  # text length dictates the box size naturally
        margin='0.1,0.05',
        height='0.4',
        color='red' if (highlight_query and is_target) else 'black',
        penwidth='3' if (highlight_query and is_target) else '1',
        fillcolor=color_map.get(row['domain'], '#ffffff'),
        fontsize=font_size,
        fontname='Consolas',
    )


# ---------------------------------------------------------------------------
# Graph assembly
# ---------------------------------------------------------------------------

def add_block_to_graph(graph, block_df, block_index, label_width, color_map,
                        highlight_query=True, collapse_opposite_strand=False,
                        font_size=10, left_pad=0, right_pad=0, spacer_width=0.6):
    """
    Add one full row (label box + every gene) to `graph`, wiring
    everything together with invisible same-rank edges so it is drawn as
    a single left-to-right row.

    `left_pad`/`right_pad` invisible spacer nodes are inserted before the
    first gene / after the last gene respectively, so that rows with
    fewer genes than the widest row still take up the same amount of
    horizontal space on each side of the query. This is what lets
    `neighborhood_figure(..., align_query_center=True)` line the query
    gene up in (approximately) the same column on every row -- "approximately"
    because real gene boxes have label-dependent widths, so this is a
    layout heuristic, not a pixel-exact guarantee. See `chain_align_nodes`
    for the second part of that trick.

    Parameters
    ----------
    graph : pygraphviz.AGraph
    block_df : pandas.DataFrame
        Rows for one block, in the left-to-right order to draw them in
        (already normalized/reversed by `normalize_block_strand` if
        that was requested).
    block_index : int
        Used to build unique node ids for this row.
    label_width : int
        From `compute_label_width`.
    color_map : dict[str, str]
    highlight_query, collapse_opposite_strand, font_size :
        Forwarded to `gene_node_style`.
    left_pad, right_pad : int
        Number of invisible spacer nodes to add on each side.
    spacer_width : float
        Width (inches) of each spacer node.

    Returns
    -------
    dict
        {'label_node': node id, 'query_node': node id of the reference
         query (see `select_reference_query_index`), or None if this
         block has no query gene, 'gene_nodes': [node ids, left-to-right],
         'gene_meta': {node id -> per-gene info dict (pid, start, end,
         strand, domain, product, plen, is_query)}}.
    """
    ref_idx = select_reference_query_index(block_df)
    if ref_idx is not None:
        query_pid = block_df.loc[ref_idx, 'pid']
        query_canonical_strand = block_df.loc[ref_idx, 'strand']
    else:
        query_pid = 'No Query'
        query_canonical_strand = 1

    label_node_id = f'label_{block_index}'
    label_html = build_row_label_html(
        query_pid=query_pid,
        block_id=block_df['ID'].iloc[0],
        org_name=block_df['org_name'].iloc[0],
        width=label_width,
        font_size=font_size,
    )
    graph.add_node(label_node_id, label=label_html, shape='none', margin=0.1)

    gene_node_ids = []
    gene_meta = {}
    query_node_id = None
    for row_position, (row_idx, row) in enumerate(block_df.iterrows()):
        node_id = f'gene_{block_index}_{row_position}'
        style = gene_node_style(
            row,
            query_canonical_strand=query_canonical_strand,
            color_map=color_map,
            highlight_query=highlight_query,
            collapse_opposite_strand=collapse_opposite_strand,
            font_size=font_size,
        )
        graph.add_node(node_id, **style)
        gene_node_ids.append(node_id)
        # Keep everything a tooltip might want; the report enriches each
        # gene node in the SVG with this (see annotate_neighborhood_svg).
        gene_meta[node_id] = dict(
            pid=row.get('pid'),
            start=row.get('start'),
            end=row.get('end'),
            strand=row.get('strand'),
            domain=row.get('domain'),
            product=row.get('product'),
            plen=row.get('plen'),
            is_query=bool(row['is_query']),
        )
        if row_idx == ref_idx:
            query_node_id = node_id

    spacer_ids_left = [f'spacer_{block_index}_L{i}' for i in range(left_pad)]
    spacer_ids_right = [f'spacer_{block_index}_R{i}' for i in range(right_pad)]
    for spacer_id in spacer_ids_left + spacer_ids_right:
        graph.add_node(spacer_id, label='', shape='box', style='invis',
                        width=spacer_width, height=0.01)

    row_node_ids = [label_node_id] + spacer_ids_left + gene_node_ids + spacer_ids_right
    graph.add_subgraph(row_node_ids, rank='same')
    for a, b in zip(row_node_ids[:-1], row_node_ids[1:]):
        graph.add_edge(a, b, style='invis', penwidth=0)

    return {'label_node': label_node_id, 'query_node': query_node_id,
            'gene_nodes': gene_node_ids, 'gene_meta': gene_meta}


def chain_align_nodes(graph, node_ids, weight=10000):
    """
    Pull a list of nodes -- one per row, in row order -- into the same
    visual column, by connecting consecutive nodes with a very high
    weight invisible edge.

    `dot` lays a graph out by minimizing total (edge weight x edge
    length); putting a large weight on an edge that is never actually
    drawn (`style='invis'`) is the standard trick to bias two nodes on
    different ranks towards the same horizontal position. The original
    code already used this for the row-label boxes; this function pulls
    that out so it can be reused for the query column too.

    Parameters
    ----------
    graph : pygraphviz.AGraph
    node_ids : list
        Skips silently over any `None` entries (e.g. a block with no
        query gene).
    weight : int
    """
    clean_ids = [n for n in node_ids if n is not None]
    for a, b in zip(clean_ids[:-1], clean_ids[1:]):
        graph.add_edge(a, b, style='invis', penwidth=0, weight=weight)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def neighborhood_figure(df, group_col='block_id', label_col='pfam', org_col='organism',
                        output_file='operon_fig_out.svg', max_colors=5,
                        highlight_query=True, font_size=10, ignore_domains=None,
                        custom_colors=None, rename_map=None, normalize_orientation=True,
                        align_query_center=False, collapse_opposite_strand=False,
                        spacer_width=0.6, color_map=None, collect_node_meta=False):
    """
    Draw a gene-neighborhood ("operon") figure, one row per block, and
    write it to `output_file`.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format table, one row per gene. See the module docstring
        for required columns.
    group_col : str
        Column that identifies which block/row a gene belongs to.
    label_col : str
        The "architecture" column -- whatever column is used for node
        labels/colors (and for `rename_map`, below). Defaults to
        'pfam', but accepts any column name (e.g. 'aravind', 'product',
        or a custom column of your own). 'pfam' (the default) uses the
        pfam -> aravind -> product fallback chain; see
        `resolve_domain_labels`. Any other value is used as-is.
    org_col : str
        Column with the organism name shown in the row label.
    output_file : str
        Path Graphviz writes the rendered figure to (extension controls
        format, e.g. '.svg', '.png', '.pdf').
    max_colors : int
        Max number of *automatically* assigned colors; see `build_color_map`.
    highlight_query : bool
        Outline every query gene in red.
    font_size : int
    ignore_domains : list[str] or None
        Domain values to never auto-color (see `build_color_map`).
    custom_colors : dict[str, str] or None
        Explicit {domain_value: color} overrides -- keys should match
        values in the 'domain' column (typically raw pfam ids when
        `label_col='pfam'`, *after* `rename_map` has been applied). See
        `build_color_map`.
    rename_map : dict[str, str] or None
        {old_name: new_name} overrides applied to the `label_col`
        column before anything else happens, so renamed values are what
        get shown, colored, and matched against `ignore_domains`/
        `custom_colors` everywhere downstream. Matches are done
        component-by-component on '+'-joined architecture strings (e.g.
        'GntR+FCD'), not only on an exact whole-string match. See
        `rename_label_values`.
    normalize_orientation : bool, default True
        If True, blocks whose reference query (see
        `select_reference_query_index`) is on the minus strand are
        mirrored so every reference query is drawn pointing the same
        way (strand +1). Set to False to draw every block in its
        original orientation.
    align_query_center : bool, default False
        If True, pad each row with invisible spacer nodes (and an extra
        alignment edge) so the reference query gene falls in roughly
        the same column on every row. The default is False: each row
        simply starts flush left (rows are not centered on the query).
        See `add_block_to_graph` and `chain_align_nodes` for the
        mechanics and its limits.
    collapse_opposite_strand : bool, default False
        If True, neighbors on the strand opposite the reference query
        are drawn as small unlabeled grey triangles pointing left,
        instead of full domain-labeled arrows. Useful for decluttering
        when those genes are not of interest. Query genes are never
        collapsed.
    spacer_width : float, default 0.6
        Width (inches) of the invisible spacer nodes used for
        `align_query_center`. Tune this if real gene boxes in your
        figure are consistently much narrower/wider than this.
    color_map : dict[str, str] or None
        Precomputed domain -> color mapping. If given, it is used
        verbatim and `build_color_map` is skipped (so `max_colors`,
        `ignore_domains` and `custom_colors` are ignored). This lets a
        caller fix one global color scheme and reuse it across several
        figures -- e.g. `build_html_report` renders one figure per
        block but wants the same domain to be the same color in every
        one, and the same as in the genome overview.

    Notes
    -----
    Every neighbor gene (anything that is not a query) is always drawn
    pointing right, regardless of its real strand -- see
    `gene_node_style`. A block with more than one query gene uses the
    one closest to the middle of the block as its reference for
    orientation/centering/labeling, but every query gene is still
    outlined in red -- see `select_reference_query_index`.

    Returns
    -------
    pandas.DataFrame, or (pandas.DataFrame, dict)
        Normally the working copy of `df` actually used to build the
        figure (with 'ID', 'org_name', 'domain', 'is_query', 'pid_order'
        added), mainly useful for debugging. If `collect_node_meta=True`,
        returns a `(working, node_meta)` tuple instead, where `node_meta`
        maps each gene node's id (e.g. 'gene_0_3') to its per-gene info
        dict -- this is what `build_html_report` uses to attach a
        hover "info window" to every protein in a neighborhood. The node
        ids match the `<title>` elements graphviz writes into the SVG.

    Other parameters
    ----------------
    collect_node_meta : bool, default False
        If True, also return the per-gene metadata dict (see Returns).
    """
    working = prepare_dataframe(
        df, group_col=group_col, org_col=org_col, label_col=label_col, rename_map=rename_map
    )
    label_width = compute_label_width(working)
    if color_map is None:
        color_map = build_color_map(
            working, max_colors=max_colors, ignore_domains=ignore_domains, custom_colors=custom_colors
        )

    blocks = [
        normalize_block_strand(block_df, normalize_orientation=normalize_orientation)
        for _, block_df in working.groupby('pid_order', sort=True)
    ]

    # First pass: how far left/right of the reference query does each
    # block extend? Needed up front so every row can be padded to the
    # same width before any nodes are added.
    left_counts, right_counts = [], []
    for block_df in blocks:
        ref_idx = select_reference_query_index(block_df)
        q_pos = block_df.index.get_loc(ref_idx) if ref_idx is not None else 0
        left_counts.append(q_pos)
        right_counts.append(len(block_df) - 1 - q_pos)
    max_left = max(left_counts) if (align_query_center and left_counts) else 0
    max_right = max(right_counts) if (align_query_center and right_counts) else 0

    graph = pgv.AGraph(directed=True)
    graph.graph_attr.update(nodesep=0.05, ranksep=0.15)

    blocks_info = []
    for block_index, (block_df, n_left, n_right) in enumerate(zip(blocks, left_counts, right_counts)):
        left_pad = (max_left - n_left) if align_query_center else 0
        right_pad = (max_right - n_right) if align_query_center else 0
        info = add_block_to_graph(
            graph, block_df, block_index, label_width, color_map,
            highlight_query=highlight_query,
            collapse_opposite_strand=collapse_opposite_strand,
            font_size=font_size,
            left_pad=left_pad,
            right_pad=right_pad,
            spacer_width=spacer_width,
        )
        blocks_info.append(info)

    # Align the label boxes into one column (original behavior)...
    chain_align_nodes(graph, [info['label_node'] for info in blocks_info])
    # ...and, optionally, the reference query genes into their own column.
    if align_query_center:
        chain_align_nodes(graph, [info['query_node'] for info in blocks_info])

    graph.draw(output_file, prog='dot')

    if collect_node_meta:
        node_meta = {}
        for info in blocks_info:
            node_meta.update(info['gene_meta'])
        return working, node_meta
    return working


# ---------------------------------------------------------------------------
# Genome-wide overview
# ---------------------------------------------------------------------------
#
# neighborhood_figure answers "what's in each neighborhood"; the
# functions below answer the complementary question, "where are these
# neighborhoods in the genome" -- one horizontal track per contig, with
# every block marked at its real genomic position.

def compute_block_extents(working, nucleotide_col='nucleotide', start_col='start',
                           end_col='end', length_col='nlen'):
    """
    Collapse a prepared gene table (see `prepare_dataframe`) down to one
    row per block: its nucleotide/contig, genomic span, contig length,
    and reference query (see `select_reference_query_index`).

    A block is assumed to sit on a single contig; its span is the
    min(start)/max(end) across all of its genes.

    Parameters
    ----------
    working : pandas.DataFrame
        Output of `prepare_dataframe` (needs 'ID', 'pid_order',
        'is_query', 'domain', 'org_name', plus `nucleotide_col`,
        `start_col`, `end_col`, and ideally `length_col`).
    nucleotide_col, start_col, end_col : str
        Per-gene columns giving its contig and genomic span.
    length_col : str
        Column with the contig's total length. If missing, or blank
        for some contig, that contig's length is approximated as the
        furthest gene/block end seen on it.

    Returns
    -------
    pandas.DataFrame
        One row per block, columns: ID, nucleotide, block_start,
        block_end, contig_length, query_pid, query_domain, org_name,
        n_genes. `query_pid`/`query_domain` are None for a block with
        no query.
    """
    records = []
    for _, block_df in working.groupby('pid_order', sort=True):
        nucleotide = block_df[nucleotide_col].iloc[0] if nucleotide_col in block_df.columns else 'Unknown'
        block_start = block_df[start_col].min() if start_col in block_df.columns else np.nan
        block_end = block_df[end_col].max() if end_col in block_df.columns else np.nan
        contig_length = block_df[length_col].iloc[0] if length_col in block_df.columns else np.nan

        ref_idx = select_reference_query_index(block_df)
        if ref_idx is not None:
            query_pid = block_df.loc[ref_idx, 'pid']
            query_domain = block_df.loc[ref_idx, 'domain']
        else:
            query_pid = None
            query_domain = None

        records.append(dict(
            ID=block_df['ID'].iloc[0],
            nucleotide=nucleotide,
            block_start=block_start,
            block_end=block_end,
            contig_length=contig_length,
            query_pid=query_pid,
            query_domain=query_domain,
            org_name=block_df['org_name'].iloc[0],
            n_genes=len(block_df),
        ))

    extents = pd.DataFrame.from_records(records)
    if not extents.empty:
        # Fall back to "furthest block end seen on this contig" for any
        # contig whose real length wasn't supplied.
        fallback_length = extents.groupby('nucleotide')['block_end'].transform('max')
        extents['contig_length'] = extents['contig_length'].fillna(fallback_length)
    return extents


def assign_label_lanes(x_positions, min_gap=280, n_lanes=3):
    """
    Stagger a set of x positions into `n_lanes` rows so that labels
    placed near each other don't overlap.

    Positions are processed left to right; each one goes into the
    first lane whose most-recently-placed position is at least
    `min_gap` away, or -- if every lane is still "busy" -- the lane
    whose last position is furthest behind (least likely to still be in
    the way). This is a simple greedy heuristic, not an optimal packing,
    but is more than enough to keep a typical handful of neighboring
    blocks legible.

    Parameters
    ----------
    x_positions : sequence of float
    min_gap : float
        Minimum spacing (in the same units as `x_positions`) before two
        labels in the same lane are considered to be clear of each other.
    n_lanes : int

    Returns
    -------
    list[int]
        Lane index (0 .. n_lanes-1) for each input position, in the
        same order as `x_positions`.
    """
    lane_last_x = [-float('inf')] * n_lanes
    lane_of = [0] * len(x_positions)
    order = sorted(range(len(x_positions)), key=lambda i: x_positions[i])

    for i in order:
        x = x_positions[i]
        free_lane = next((lane for lane in range(n_lanes) if x - lane_last_x[lane] >= min_gap), None)
        chosen = free_lane if free_lane is not None else min(range(n_lanes), key=lambda l: lane_last_x[l])
        lane_last_x[chosen] = x
        lane_of[i] = chosen

    return lane_of


def build_genome_overview_svg(extents, color_map=None, highlight_color='#c0392b',
                               marker_color='#2a6f77', track_width=760, left_margin=190,
                               top_margin=30, row_height=92, track_height=10,
                               font_size=11, label_lanes=3, label_lane_gap=18,
                               max_labels_per_track=25):
    """
    Render a *static* SVG showing where every block/neighborhood sits
    along its nucleotide (contig), one horizontal track per distinct
    nucleotide.

    This is the standalone/fallback renderer (what `genome_overview_fig`
    writes to a file). For a crowded genome the interactive version in
    `build_html_report` is far more usable -- this static one can only
    fit so many text labels before they collide, which is exactly why
    `max_labels_per_track` exists.

    Each track is scaled independently to its own `contig_length` -- a
    50 kb plasmid and a 9 Mb chromosome both draw at the same pixel
    width -- since the point of this figure is "where is this block
    relative to its own contig", not a comparison of absolute distance
    across contigs of very different sizes.

    Parameters
    ----------
    extents : pandas.DataFrame
        Output of `compute_block_extents`.
    color_map : dict[str, str] or None
        Domain -> color (e.g. from `build_color_map`), used to color
        each marker by its reference query's domain, for visual
        consistency with `neighborhood_figure`. A block whose domain
        has no entry (or no `color_map` given) uses `marker_color`.
    highlight_color : str
        Border color for every block marker.
    marker_color : str
        Fallback marker fill.
    track_width, left_margin, top_margin, row_height, track_height : float
        Layout, in SVG user units (effectively pixels).
    font_size : int
    label_lanes, label_lane_gap : int, float
        Up to this many staggered rows are used above each track for
        block labels; see `assign_label_lanes`. Increase `row_height`
        if labels still collide with the row above.
    max_labels_per_track : int
        If a track has more than this many blocks, its text labels are
        omitted entirely (markers are still drawn) -- a crowded track's
        labels just turn into noise, as in a whole-chromosome view with
        hundreds of hits. Set very high to always label.

    Returns
    -------
    str
        A full, self-contained `<svg>...</svg>` document.
    """
    color_map = color_map or {}
    nucleotides = list(dict.fromkeys(extents['nucleotide'])) if not extents.empty else []

    fig_width = left_margin + track_width + 40
    fig_height = top_margin + len(nucleotides) * row_height + 20

    parts = [
        f'<svg viewBox="0 0 {fig_width:.0f} {fig_height:.0f}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="Consolas, \'SF Mono\', Menlo, monospace" font-size="{font_size}">',
        f'<rect x="0" y="0" width="{fig_width:.0f}" height="{fig_height:.0f}" fill="white"/>',
    ]

    for row_i, nucleotide in enumerate(nucleotides):
        row_blocks = extents[extents['nucleotide'] == nucleotide]
        contig_length = row_blocks['contig_length'].iloc[0]
        if not contig_length or pd.isna(contig_length) or contig_length <= 0:
            contig_length = max(row_blocks['block_end'].max(), 1)

        track_y = top_margin + row_i * row_height + label_lanes * label_lane_gap
        track_x0 = left_margin

        def to_x(pos, _x0=track_x0, _len=contig_length):
            return _x0 + (pos / _len) * track_width

        parts.append(
            f'<text x="{track_x0 - 10:.0f}" y="{track_y + track_height / 2 + 4:.0f}" '
            f'text-anchor="end" fill="#222">{html.escape(str(nucleotide))}</text>'
        )
        parts.append(
            f'<text x="{track_x0 - 10:.0f}" y="{track_y + track_height / 2 + 4 + font_size + 2:.0f}" '
            f'text-anchor="end" fill="#888" font-size="{font_size - 2}">{contig_length:,.0f} bp</text>'
        )
        parts.append(
            f'<rect x="{track_x0:.1f}" y="{track_y:.1f}" width="{track_width:.1f}" height="{track_height:.1f}" '
            f'fill="#e3e3e3" stroke="#999" stroke-width="0.5" rx="2"/>'
        )

        show_labels = len(row_blocks) <= max_labels_per_track
        mid_x = [to_x((b['block_start'] + b['block_end']) / 2) for _, b in row_blocks.iterrows()]
        lanes = assign_label_lanes(mid_x, min_gap=label_lane_gap * 12, n_lanes=label_lanes)

        for (_, block), x_mid, lane in zip(row_blocks.iterrows(), mid_x, lanes):
            x0 = to_x(block['block_start'])
            x1 = to_x(block['block_end'])
            marker_w = max(4.0, x1 - x0)
            fill = color_map.get(block['query_domain'], marker_color)

            parts.append(
                f'<rect x="{x0:.1f}" y="{track_y - 3:.1f}" width="{marker_w:.1f}" '
                f'height="{track_height + 6:.1f}" fill="{fill}" stroke="{highlight_color}" '
                f'stroke-width="1.2" rx="1.5"/>'
            )

            if not show_labels:
                continue
            label_y = track_y - 10 - lane * label_lane_gap
            parts.append(
                f'<line x1="{x_mid:.1f}" y1="{label_y + 4:.1f}" x2="{x_mid:.1f}" y2="{track_y - 3:.1f}" '
                f'stroke="#bbb" stroke-width="1"/>'
            )
            label_text = block['query_pid'] if block['query_pid'] is not None else block['ID']
            parts.append(
                f'<text x="{x_mid:.1f}" y="{label_y:.1f}" text-anchor="middle" fill="#222">'
                f'{html.escape(str(label_text))}</text>'
            )

    parts.append('</svg>')
    return '\n'.join(parts)


def _slug(text):
    """
    Turn an arbitrary block id into a string safe to use as an HTML
    `id`/`data-` value (letters, digits, dash, underscore only). Used
    to link a genome-overview marker to its per-block panel in the
    report. Not meant to be reversible -- just stable and collision-
    resistant enough for the block ids seen here.
    """
    out = ''.join(c if (c.isalnum() or c in '-_') else '-' for c in str(text))
    return out or 'block'


def build_genome_overview_interactive_html(extents, color_map=None,
                                            marker_color='#2a6f77',
                                            highlight_color='#c0392b',
                                            base_track_height=14):
    """
    Build the *interactive* genome overview used by `build_html_report`:
    one horizontal track per nucleotide (contig), where each block is a
    marker positioned at its genomic span, and where the surrounding
    report provides zoom/pan and hover tooltips.

    Unlike `build_genome_overview_svg`, this does NOT bake any text
    labels into the figure -- that is exactly what turns a
    whole-chromosome view with hundreds of hits into unreadable noise.
    Instead every marker carries its info in a `data-tip` attribute that
    the report's JavaScript shows on hover, and a `data-block` slug that
    links it to that block's own panel in the neighborhoods section
    (click a marker to jump to it).

    Each marker stores its position as fractions of its contig length
    (`data-fs`/`data-fe`), so the report's JS can re-place every marker
    at any zoom level without this function knowing the final pixel
    width.

    Parameters
    ----------
    extents : pandas.DataFrame
        Output of `compute_block_extents` (needs the `n_genes` column).
    color_map : dict[str, str] or None
        Domain -> color; a block whose reference-query domain is absent
        falls back to `marker_color`.
    marker_color, highlight_color : str
        Fallback fill and the marker border color.
    base_track_height : int
        Marker/track height in pixels.

    Returns
    -------
    str
        An HTML fragment (a `<div class="go-wrap">...`). It depends on
        the CSS/JS that `build_html_report` injects, so it is not
        standalone on its own.
    """
    color_map = color_map or {}
    nucleotides = list(dict.fromkeys(extents['nucleotide'])) if not extents.empty else []

    rows = []
    for nucleotide in nucleotides:
        row_blocks = extents[extents['nucleotide'] == nucleotide]
        contig_length = row_blocks['contig_length'].iloc[0]
        if not contig_length or pd.isna(contig_length) or contig_length <= 0:
            contig_length = max(row_blocks['block_end'].max(), 1)

        markers = []
        for _, block in row_blocks.iterrows():
            fs = max(0.0, min(1.0, block['block_start'] / contig_length))
            fe = max(0.0, min(1.0, block['block_end'] / contig_length))
            if fe < fs:
                fs, fe = fe, fs
            fill = color_map.get(block['query_domain'], marker_color)

            label = block['query_pid'] if block['query_pid'] is not None else block['ID']
            domain = block['query_domain'] if block['query_domain'] is not None else '-'
            tip_html = (
                f"<b>{html.escape(str(label))}</b>"
                f"<span class='t-row'>block&nbsp;&middot;&nbsp;{html.escape(str(block['ID']))}</span>"
                f"<span class='t-row'>domain&nbsp;&middot;&nbsp;{html.escape(str(domain))}</span>"
                f"<span class='t-row'>organism&nbsp;&middot;&nbsp;{html.escape(str(block['org_name']))}</span>"
                f"<span class='t-row'>position&nbsp;&middot;&nbsp;{block['block_start']:,.0f}&ndash;{block['block_end']:,.0f} bp</span>"
                f"<span class='t-row'>genes&nbsp;&middot;&nbsp;{int(block['n_genes'])}</span>"
            )
            tip_attr = html.escape(tip_html, quote=True)

            markers.append(
                f'<div class="go-marker" data-block="{_slug(block["ID"])}" '
                f'data-fs="{fs:.6f}" data-fe="{fe:.6f}" data-tip="{tip_attr}" '
                f'style="background:{fill};border-color:{highlight_color};"></div>'
            )

        rows.append(
            '<div class="go-track" style="--track-h:%dpx;">'
            '<div class="go-track-label"><div class="go-contig">%s</div>'
            '<div class="go-len">%s bp</div></div>'
            '<div class="go-viewport"><div class="go-inner"><div class="go-axis"></div>%s</div></div>'
            '</div>' % (
                base_track_height,
                html.escape(str(nucleotide)),
                f'{contig_length:,.0f}',
                ''.join(markers),
            )
        )

    controls = (
        '<div class="go-controls">'
        '<button type="button" id="go-zoom-out" title="Zoom out">&minus;</button>'
        '<span id="go-zoom-val">1x</span>'
        '<button type="button" id="go-zoom-in" title="Zoom in">+</button>'
        '<button type="button" id="go-zoom-reset" title="Reset zoom and pan">reset</button>'
        '<span class="go-hint">scroll to pan &middot; ctrl/&#8984;+scroll or buttons to zoom &middot; click a marker to open it</span>'
        '</div>'
    )

    return f'<div class="go-wrap">{controls}{"".join(rows)}</div><div id="go-tooltip" class="go-tooltip"></div>'


def genome_overview_fig(df, group_col='block_id', org_col='organism', label_col='pfam',
                         rename_map=None, nucleotide_col='nucleotide', start_col='start',
                         end_col='end', length_col='nlen', output_file='genome_overview.svg',
                         custom_colors=None, max_colors=5, ignore_domains=None,
                         color_map=None, **svg_kwargs):
    """
    Draw a genome-wide overview: one horizontal track per nucleotide
    (contig), with every block/neighborhood marked at its position
    along that contig.

    `df`, `group_col`, `org_col`, `label_col`, `rename_map`,
    `custom_colors`, `max_colors` and `ignore_domains` mean exactly what
    they mean in `neighborhood_figure` -- pass the same values to both if
    you want the two figures to agree on domain colors/renamed names.

    This writes the *static* SVG version (see `build_genome_overview_svg`).
    The interactive, zoomable version with hover tooltips lives in
    `build_html_report`.

    Parameters
    ----------
    nucleotide_col, start_col, end_col : str
        Columns giving each gene's contig name and genomic span.
    length_col : str
        Column with the contig's total length, used to scale each
        track; see `compute_block_extents` for the fallback when it's
        missing.
    output_file : str
        Path to write the SVG to.
    color_map : dict[str, str] or None
        Precomputed domain -> color mapping; if given, `build_color_map`
        is skipped. See the same parameter on `neighborhood_figure`.
    **svg_kwargs :
        Forwarded to `build_genome_overview_svg` (layout/color tuning,
        e.g. `track_width`, `row_height`, `marker_color`,
        `max_labels_per_track`).

    Returns
    -------
    pandas.DataFrame
        One row per block; see `compute_block_extents`.
    """
    working = prepare_dataframe(
        df, group_col=group_col, org_col=org_col, label_col=label_col, rename_map=rename_map
    )
    if color_map is None:
        color_map = build_color_map(
            working, max_colors=max_colors, ignore_domains=ignore_domains, custom_colors=custom_colors
        )
    extents = compute_block_extents(
        working, nucleotide_col=nucleotide_col, start_col=start_col,
        end_col=end_col, length_col=length_col
    )
    svg = build_genome_overview_svg(extents, color_map=color_map, **svg_kwargs)

    with open(output_file, 'w') as f:
        f.write(svg)

    return extents


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def _strip_svg_prolog(svg):
    """
    Return `svg` starting at its first `<svg` tag, dropping any XML
    declaration / DOCTYPE that precedes it.

    Graphviz writes a full XML document (`<?xml ...?>` + `<!DOCTYPE ...>`
    before `<svg>`); those are invalid inline in an HTML body and some
    browsers render them as stray text, so they are stripped before an
    SVG is embedded into the report.
    """
    idx = svg.find('<svg')
    return svg[idx:] if idx != -1 else svg


def _fmt_int(value):
    """Format a number with thousands separators, or '?' if missing."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return '?'
    try:
        return f'{int(round(float(value))):,}'
    except (TypeError, ValueError):
        return html.escape(str(value))


def build_gene_tooltip_html(meta):
    """
    Build the inner HTML of the hover "info window" for a single protein
    in a neighborhood figure.

    The title is the protein id. Then a role line:

      * the query gene shows a "query" marker (this is the
        "change domain for query" behavior -- where a neighbor would
        list its domain, the query instead announces that it *is* the
        query);
      * a neighbor lists its domain/architecture.

    Followed by genomic coordinates, strand, length and product when
    those fields are available.

    Parameters
    ----------
    meta : dict
        One value from the `node_meta` dict returned by
        `neighborhood_figure(..., collect_node_meta=True)`.

    Returns
    -------
    str
        HTML for the tooltip body (same shape as the genome-overview
        tooltips: a bold title plus `.t-row` spans).
    """
    pid = meta.get('pid')
    title = pid if pid is not None and not (isinstance(pid, float) and pd.isna(pid)) else 'protein'

    rows = [f"<b>{html.escape(str(title))}</b>"]

    if meta.get('is_query'):
        rows.append("<span class='t-row t-query'>&#9733;&nbsp;query</span>")
    else:
        domain = meta.get('domain')
        domain = '-' if domain is None or (isinstance(domain, float) and pd.isna(domain)) else domain
        rows.append(f"<span class='t-row'>domain&nbsp;&middot;&nbsp;{html.escape(str(domain))}</span>")

    start, end = meta.get('start'), meta.get('end')
    if start is not None or end is not None:
        rows.append(
            f"<span class='t-row'>coords&nbsp;&middot;&nbsp;{_fmt_int(start)}&ndash;{_fmt_int(end)} bp</span>"
        )

    strand = meta.get('strand')
    strand_str = '+' if strand == 1 else '&minus;' if strand == -1 else '?'
    rows.append(f"<span class='t-row'>strand&nbsp;&middot;&nbsp;{strand_str}</span>")

    plen = meta.get('plen')
    if plen is not None and not (isinstance(plen, float) and pd.isna(plen)):
        rows.append(f"<span class='t-row'>length&nbsp;&middot;&nbsp;{_fmt_int(plen)} aa</span>")

    product = meta.get('product')
    if product is not None and not (isinstance(product, float) and pd.isna(product)):
        rows.append(f"<span class='t-row'>product&nbsp;&middot;&nbsp;{html.escape(str(product))}</span>")

    return ''.join(rows)


# Matches one graphviz node group header: `<g id="nodeN" class="node">`
# immediately followed by `<title>NODE_NAME</title>`. The node name is
# what `add_block_to_graph` assigned (e.g. 'gene_0_3'); graphviz writes
# it into the <title>, which is how we map an SVG element back to its
# gene metadata.
_NODE_GROUP_RE = re.compile(
    r'<g id="(?P<gid>[^"]*)" class="node">\s*<title>(?P<name>[^<]+)</title>'
)


def annotate_neighborhood_svg(svg, node_meta):
    """
    Enrich each gene node in a rendered neighborhood SVG so the report
    can show a per-protein info window on hover.

    For every `<g class="node">` whose `<title>` names a gene that has
    metadata, this adds `class="node nb-gene"` and a `data-tip`
    attribute (the HTML from `build_gene_tooltip_html`). Label boxes,
    spacers and any node without metadata are left untouched.

    This is plain text surgery on graphviz's SVG output rather than a
    graphviz feature: graphviz can attach a `tooltip` (which becomes a
    slow native browser tooltip) but not arbitrary `data-*` attributes
    or a styled popup, so the report injects its own.

    Parameters
    ----------
    svg : str
        A rendered neighborhood SVG (ideally prolog-stripped already).
    node_meta : dict
        node id -> per-gene info, from
        `neighborhood_figure(..., collect_node_meta=True)`.

    Returns
    -------
    str
        The SVG with gene nodes annotated.
    """
    def repl(match):
        name = match.group('name')
        meta = node_meta.get(name)
        if not meta:
            return match.group(0)
        tip = html.escape(build_gene_tooltip_html(meta), quote=True)
        gid = match.group('gid')
        new_header = f'<g id="{gid}" class="node nb-gene" data-tip="{tip}">'
        return match.group(0).replace(f'<g id="{gid}" class="node">', new_header, 1)

    return _NODE_GROUP_RE.sub(repl, svg)


def render_neighborhood_svgs_by_block(df, group_col, color_map, operon_kwargs, tmp_dir):
    """
    Render one neighborhood figure per block and return their SVGs.

    Each block is rendered on its own (the df filtered to just that
    block), so it can be viewed independently in the report's tabbed
    neighborhoods section. A shared `color_map` is passed to
    every render so the same domain is the same color in every block's
    figure (and in the genome overview).

    Parameters
    ----------
    df : pandas.DataFrame
        The full raw input table.
    group_col : str
        Column whose distinct values define blocks.
    color_map : dict[str, str]
        Shared domain -> color mapping (see `build_color_map`).
    operon_kwargs : dict
        Extra keyword arguments for `neighborhood_figure`. `output_file`,
        `color_map` and `group_col` are managed here and ignored if
        present.
    tmp_dir : str
        Directory the intermediate SVGs are written to.

    Returns
    -------
    dict[str, str]
        Maps each block's slug (see `_slug`) to its SVG markup, with XML
        prolog stripped and every protein annotated for hover info
        windows (see `annotate_neighborhood_svg`). Insertion order
        follows the first appearance of each block in `df`.
    """
    kw = {k: v for k, v in (operon_kwargs or {}).items()
          if k not in ('output_file', 'color_map', 'group_col', 'collect_node_meta')}

    svgs = {}
    for block_id, block_df in df.groupby(group_col, sort=False):
        slug = _slug(block_id)
        path = os.path.join(tmp_dir, f'nb_{slug}.svg')
        _working, node_meta = neighborhood_figure(
            block_df, group_col=group_col, output_file=path,
            color_map=color_map, collect_node_meta=True, **kw)
        with open(path) as f:
            svg = _strip_svg_prolog(f.read())
        svgs[slug] = annotate_neighborhood_svg(svg, node_meta)
    return svgs


def render_dataframe_html(df, table_id='data-table', max_rows=None, css_class=None):
    """
    Render `df` as a plain HTML `<table>` (every value HTML-escaped),
    suitable for dropping into a larger page.

    Parameters
    ----------
    df : pandas.DataFrame
    table_id : str or None
        `id` attribute on the `<table>`, so the rest of a page (CSS, a
        search box's JS) can target it. `build_html_report`'s main
        search box expects the default, 'data-table'. Pass None to omit
        the id (used for the per-neighborhood sub-tables, so they don't
        collide with the main table's id).
    max_rows : int or None
        If given, only the first `max_rows` rows are rendered, with a
        note below the table saying how many were left out. `None`
        (the default) renders every row.
    css_class : str or None
        Optional `class` attribute on the `<table>`.

    Returns
    -------
    str
        `<table>...</table>` markup (plus a trailing `<p>` note if
        `max_rows` truncated anything).
    """
    shown = df if max_rows is None else df.head(max_rows)

    header_cells = ''.join(f'<th>{html.escape(str(c))}</th>' for c in shown.columns)
    body_rows = (
        '<tr>' + ''.join('<td>' + ('' if pd.isna(v) else html.escape(str(v))) + '</td>' for v in row) + '</tr>'
        for row in shown.itertuples(index=False, name=None)
    )

    attrs = ''
    if table_id:
        attrs += f' id="{table_id}"'
    if css_class:
        attrs += f' class="{css_class}"'
    table_html = (
        f'<table{attrs}>'
        f'<thead><tr>{header_cells}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody>'
        f'</table>'
    )
    if max_rows is not None and len(df) > max_rows:
        table_html += f'<p class="table-note">Showing the first {max_rows:,} of {len(df):,} rows.</p>'
    return table_html


def render_table_card(df, filename='table.csv', max_rows=None):
    """
    Wrap a dataframe in an interactive "table card": a small toolbar
    (text filter, live row count, "download CSV" button) above the
    table itself. The report's JavaScript makes every such card's
    columns click-to-sort and its button export the (filtered) rows.

    Parameters
    ----------
    df : pandas.DataFrame
    filename : str
        Suggested name for the CSV download.
    max_rows : int or None
        Passed through to `render_dataframe_html`.

    Returns
    -------
    str
        HTML for one `.tbl-card`.
    """
    table = render_dataframe_html(df, table_id=None, max_rows=max_rows, css_class='data-tbl')
    return (
        f'<div class="tbl-card" data-filename="{html.escape(filename)}">'
        '<div class="tbl-controls">'
        '<input type="text" class="tbl-filter" placeholder="Filter rows...">'
        '<span class="tbl-count"></span>'
        '<div class="tbl-dl-wrap">''<button type="button" class="tbl-download tbl-dl-btn" data-fmt="csv">&#8681; CSV</button>''<button type="button" class="tbl-download tbl-dl-btn" data-fmt="tsv">TSV</button>''<button type="button" class="tbl-download tbl-dl-btn" data-fmt="json">JSON</button>''</div>'
        '</div>'
        f'<div class="tbl-scroll">{table}</div>'
        '</div>'
    )


def compute_domain_stats(working, extents, ignore_domains=None):
    """
    Count how often each domain appears, two ways:

      * "reference" -- the domain of each block's reference query gene,
        i.e. one count per neighborhood (how many neighborhoods are
        built around each domain);
      * "architecture" -- every domain across every protein, splitting
        '+'-joined architecture strings into their components (how often
        each domain shows up anywhere in the data).

    Generic/uninformative values (`ignore_domains`, plus 'unk'/'-'/'?'
    and anything containing "hypothetical") are dropped from both.

    Parameters
    ----------
    working : pandas.DataFrame
        Prepared table (needs the 'domain' column); see
        `prepare_dataframe`.
    extents : pandas.DataFrame
        Per-block table (needs 'query_domain'); see
        `compute_block_extents`.
    ignore_domains : list[str] or None
        Extra values to drop. Defaults to `DEFAULT_IGNORE_DOMAINS`.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        (reference_counts, architecture_counts), each with columns
        'domain' and 'count', sorted by count descending.
    """
    ignore = {d.lower() for d in (ignore_domains or DEFAULT_IGNORE_DOMAINS)} | {'unk', '-', '?', ''}

    def keep(value):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return False
        s = str(value).lower()
        return s not in ignore and 'hypothetical' not in s

    ref_values = [d for d in extents['query_domain'].tolist() if keep(d)]
    ref_counts = (
        pd.Series(ref_values, dtype=object).value_counts()
        .rename_axis('domain').reset_index(name='count')
    )

    components = []
    for value in working['domain'].tolist():
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        for part in str(value).split('+'):
            if keep(part):
                components.append(part)
    arch_counts = (
        pd.Series(components, dtype=object).value_counts()
        .rename_axis('domain').reset_index(name='count')
    )

    return ref_counts, arch_counts


def build_bar_chart_svg(counts_df, color_map=None, top=15, marker_color='#2a6f77',
                         width=680, label_w=210, row_h=26, pad=14, font_size=12):
    """
    Render a horizontal bar chart (SVG) of the `top` most frequent
    domains in a `{domain, count}` table, colored to match the figures.

    Returns a self-contained `<svg>` string (or a small "nothing to
    show" message if the table is empty).
    """
    color_map = color_map or {}
    if counts_df.empty:
        return '<p class="nb-empty">No domains to summarize.</p>'

    shown = counts_df.head(top)
    max_count = int(shown['count'].max())
    bar_area = width - label_w - 64
    height = pad * 2 + len(shown) * row_h

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="Consolas, \'SF Mono\', Menlo, monospace" font-size="{font_size}">'
    ]
    for i, (_, row) in enumerate(shown.iterrows()):
        y = pad + i * row_h
        cy = y + row_h / 2
        domain = str(row['domain'])
        count = int(row['count'])
        w = (count / max_count) * bar_area if max_count else 0
        fill = color_map.get(domain, marker_color)

        shown_label = domain if len(domain) <= 26 else domain[:24] + '\u2026'
        parts.append(
            f'<text x="{label_w - 8:.0f}" y="{cy + 4:.0f}" text-anchor="end" fill="#222">'
            f'{html.escape(shown_label)}</text>'
        )
        parts.append(
            f'<rect x="{label_w:.0f}" y="{y + 4:.0f}" width="{max(w, 1):.1f}" height="{row_h - 10:.0f}" '
            f'fill="{fill}" stroke="#c0392b" stroke-width="0.6" rx="2"/>'
        )
        parts.append(
            f'<text x="{label_w + w + 8:.0f}" y="{cy + 4:.0f}" fill="#555">{count}</text>'
        )
    parts.append('</svg>')
    return '\n'.join(parts)


def build_stats_section_html(ref_counts, arch_counts, color_map=None):
    """
    Build the statistics section's inner HTML: a toggle between
    "Reference domains" and "Architectures", each showing a bar chart of
    the most frequent domains plus a full sortable/downloadable table.

    Parameters
    ----------
    ref_counts, arch_counts : pandas.DataFrame
        From `compute_domain_stats`.
    color_map : dict[str, str] or None
        Domain -> color, to match the figures.

    Returns
    -------
    str
        HTML fragment for the statistics section body.
    """
    ref_chart = build_bar_chart_svg(ref_counts, color_map=color_map)
    arch_chart = build_bar_chart_svg(arch_counts, color_map=color_map)
    ref_card = render_table_card(ref_counts, filename='reference_domain_counts.csv')
    arch_card = render_table_card(arch_counts, filename='architecture_domain_counts.csv')

    n_ref = len(ref_counts)
    n_arch = len(arch_counts)

    return (
        '<div class="stats-toggle">'
        '<button type="button" class="stats-btn active" data-stats="ref">Reference domains</button>'
        '<button type="button" class="stats-btn" data-stats="arch">Architectures</button>'
        '</div>'
        '<div class="stats-block active" data-stats="ref">'
        f'<p class="stats-note">Counts the domain of each neighborhood\'s reference query &mdash; {n_ref} distinct domains across the queries.</p>'
        f'<div class="stats-chart">{ref_chart}</div>{ref_card}</div>'
        '<div class="stats-block" data-stats="arch">'
        f'<p class="stats-note">Counts every domain across every protein, splitting \'+\'-joined architectures &mdash; {n_arch} distinct domains.</p>'
        f'<div class="stats-chart">{arch_chart}</div>{arch_card}</div>'
    )


# Kept as a module-level constant (rather than building the string
# inline inside `build_html_report`) so the template can be read,
# tweaked, or unit-tested on its own. Uses `string.Template`'s
# `$placeholder` syntax rather than `str.format`'s `{placeholder}`,
# since the CSS/JS below are full of literal `{`/`}` braces that would
# otherwise have to be escaped throughout. (Substituted values are
# inserted verbatim and are NOT re-scanned for `$`, so embedded SVG or
# table data containing `$` is safe -- only this literal template text
# must avoid stray `$`.)
HTML_REPORT_TEMPLATE = Template(r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>$title</title>
<style>
  :root {
    --bg: #f5f4f0;
    --panel: #ffffff;
    --ink: #1f2430;
    --muted: #6b7280;
    --line: #e4e1d8;
    --accent: #2a6f77;
    --accent-soft: #e4f0f1;
    --selected: #fdecdc;
    --selected-line: #e07b39;
  }
  *{box-sizing:border-box;}
  body{margin:0;background:var(--bg);color:var(--ink);
    font-family:-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;line-height:1.45;}

  /* ---- top nav ---- */
  .top-nav{
    display:flex;align-items:stretch;border-bottom:2px solid var(--accent);
    background:var(--panel);position:sticky;top:0;z-index:400;
  }
  .top-brand{
    display:flex;align-items:center;gap:10px;padding:0 20px;
    font-size:15px;font-weight:700;letter-spacing:-.01em;color:var(--accent);
    border-right:1px solid var(--line);min-width:140px;
  }
  .top-brand .brand-sub{font-size:11px;font-weight:400;color:var(--muted);margin-left:4px;}
  /* header logo container: fixed height, width scales via viewBox */
  .top-logo-wrap{height:40px;flex-shrink:0;display:flex;align-items:center;}
  .top-logo-wrap svg{height:40px;width:auto;display:block;}
  .brand-name{font-size:15px;font-weight:700;color:var(--accent);letter-spacing:-.01em;}
  /* footer logo container — deliberately generous so the emblem reads clearly */
  .footer-logo-wrap{height:110px;display:flex;align-items:center;justify-content:center;margin-bottom:12px;}
  .footer-logo-wrap svg{height:110px;width:auto;display:block;}
  .top-tabs{display:flex;flex:1;gap:0;overflow-x:auto;}
  .top-tab{
    display:flex;align-items:center;gap:7px;padding:13px 20px;
    font-size:13px;cursor:pointer;border:none;background:none;color:var(--muted);
    border-bottom:3px solid transparent;margin-bottom:-2px;white-space:nowrap;
    transition:color .15s,border-color .15s;
  }
  .top-tab:hover{color:var(--ink);}
  .top-tab.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:600;}
  .top-tab .tab-icon{font-size:15px;opacity:.7;}
  .top-tab .tab-badge{
    font-size:10px;background:var(--accent-soft);color:var(--accent);
    border-radius:10px;padding:1px 6px;font-weight:600;
  }
  .top-meta{
    display:flex;align-items:center;padding:0 22px;font-size:12px;
    color:var(--muted);border-left:1px solid var(--line);gap:14px;white-space:nowrap;
  }
  .top-meta b{color:var(--accent);}

  /* ---- page sections (one active at a time) ---- */
  .page-section{display:none;}
  .page-section.active{display:block;}

  /* ---- shared layout ---- */
  .sec-inner{max-width:1280px;margin:0 auto;padding:28px 36px 72px;}
  .sec-title{font-size:22px;font-weight:700;margin:0 0 4px;}
  .sec-desc{color:var(--muted);font-size:14px;margin:0 0 22px;max-width:780px;}
  .panel{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:20px;}

  /* ---- tooltips (genome markers + proteins) ---- */
  .go-tooltip{
    position:fixed;display:none;z-index:1000;pointer-events:none;
    background:#1f2430;color:#fff;border-radius:8px;padding:10px 13px;
    font-size:12px;max-width:330px;box-shadow:0 6px 22px rgba(0,0,0,.28);
    font-family:Consolas,"SF Mono",Menlo,monospace;line-height:1.5;
  }
  .go-tooltip b{display:block;margin-bottom:3px;font-size:12.5px;}
  .go-tooltip .t-row{display:block;color:#cdd3da;}
  .go-tooltip .t-query{color:#ffd28a;font-weight:700;}
  .nb-gene{cursor:pointer;}
  .nb-gene:hover polygon,.nb-gene:hover ellipse{stroke-width:2.4px;}

  /* ---- genome overview ---- */
  .go-controls{display:flex;align-items:center;gap:8px;margin-bottom:16px;flex-wrap:wrap;}
  .go-controls button{
    border:1px solid var(--line);background:#fff;color:var(--ink);
    border-radius:6px;padding:5px 12px;font-size:13px;cursor:pointer;
  }
  .go-controls button:hover{background:var(--accent-soft);border-color:var(--accent);}
  #go-zoom-val{font-family:Consolas,"SF Mono",Menlo,monospace;font-size:13px;min-width:38px;text-align:center;}
  .go-hint{color:var(--muted);font-size:12px;margin-left:6px;}
  .go-track{display:grid;grid-template-columns:150px 1fr;align-items:center;margin:14px 0;gap:12px;}
  .go-track-label{text-align:right;font-family:Consolas,"SF Mono",Menlo,monospace;overflow:hidden;}
  .go-contig{font-size:13px;color:#222;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
  .go-len{font-size:11px;color:var(--muted);}
  .go-viewport{overflow-x:auto;overflow-y:hidden;cursor:grab;padding:18px 0;}
  .go-viewport.grabbing{cursor:grabbing;}
  .go-inner{position:relative;height:var(--track-h,14px);}
  .go-axis{
    position:absolute;left:0;right:0;top:50%;transform:translateY(-50%);
    height:6px;background:#e3e3e3;border:.5px solid #bbb;border-radius:3px;
  }
  .go-marker{
    position:absolute;top:0;height:100%;border:1.2px solid;
    border-radius:2px;cursor:pointer;transition:transform .06s ease;
  }
  .go-marker:hover{transform:scaleY(1.45);z-index:3;}
  .go-marker.selected{box-shadow:0 0 0 2px var(--selected-line);z-index:2;}

  /* ---- neighborhoods: toolbar ---- */
  .nb-toolbar{
    display:flex;align-items:center;gap:8px;flex-wrap:wrap;
    margin-bottom:16px;
  }
  .nb-icon-btn{
    display:flex;align-items:center;gap:6px;
    border:1px solid var(--line);background:#fff;color:var(--ink);
    border-radius:8px;padding:7px 14px;font-size:13px;cursor:pointer;
    transition:background .12s,border-color .12s;
  }
  .nb-icon-btn:hover{background:var(--accent-soft);border-color:var(--accent);}
  .nb-icon-btn svg{width:15px;height:15px;stroke:currentColor;fill:none;stroke-width:2;stroke-linecap:round;}
  .nb-zoom-group{display:flex;align-items:center;gap:0;margin-left:auto;
    border:1px solid var(--line);border-radius:8px;overflow:hidden;}
  .nb-zoom-group button{
    border:none;background:#fff;color:var(--ink);
    padding:7px 11px;font-size:15px;cursor:pointer;border-right:1px solid var(--line);
  }
  .nb-zoom-group button:last-child{border-right:none;}
  .nb-zoom-group button:hover{background:var(--accent-soft);}
  #nb-zoom-val{
    font-family:Consolas,"SF Mono",Menlo,monospace;font-size:12px;
    min-width:44px;text-align:center;padding:0 4px;background:#fff;border-right:1px solid var(--line);
  }

  /* current-block breadcrumb */
  .nb-crumb{
    font-family:Consolas,"SF Mono",Menlo,monospace;font-size:12.5px;
    color:var(--muted);margin-bottom:12px;display:flex;align-items:center;gap:8px;flex-wrap:wrap;
  }
  .nb-crumb b{color:var(--ink);}
  .nb-crumb .nb-nav-arrows{display:flex;gap:4px;margin-left:auto;}
  .nb-crumb .nb-nav-arrows button{
    border:1px solid var(--line);background:#fff;border-radius:6px;
    padding:3px 9px;font-size:13px;cursor:pointer;
  }
  .nb-crumb .nb-nav-arrows button:hover{background:var(--accent-soft);}

  /* sub-tabs inside each panel (Figure / Table) */
  /* per-panel Figure/Table sub-tabs */
  .nb-panel-subtabs{display:flex;gap:0;border-bottom:2px solid var(--line);margin-bottom:14px;}
  .nb-subtab{
    padding:6px 16px;font-size:13px;cursor:pointer;border:none;background:none;
    color:var(--muted);border-bottom:3px solid transparent;margin-bottom:-2px;
  }
  .nb-subtab.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:600;}
  .nb-subtab:hover{color:var(--ink);}
  /* multiple active panels stack with a divider */
  .nb-panel.active+.nb-panel.active{margin-top:20px;padding-top:20px;border-top:2px solid var(--line);}

  /* panels */
  .nb-panel{display:none;}
  .nb-panel.active{display:block;}
  .nb-view{display:none;}
  .nb-view.active{display:block;}
  .nb-panel-head{
    font-family:Consolas,"SF Mono",Menlo,monospace;font-size:12.5px;
    color:var(--muted);margin-bottom:10px;
  }
  .nb-panel-head b{color:var(--ink);}

  /* figure wrapper: full width, zoom stretches it */
  .nb-fig-scroll{overflow-x:auto;border:1px solid var(--line);border-radius:8px;padding:16px;background:#fff;}
  .nb-fig{width:100%;}
  .nb-fig svg{display:block;width:100%;height:auto;}
  .nb-empty{color:var(--muted);font-size:13px;padding:20px;}

  /* ---- selector pop-up ---- */
  .nb-sel-backdrop{
    display:none;position:fixed;inset:0;background:rgba(20,24,32,.38);
    z-index:900;align-items:flex-start;justify-content:center;padding-top:60px;
  }
  .nb-sel-backdrop.open{display:flex;}
  .nb-sel-modal{
    background:#fff;border-radius:14px;width:min(480px,94vw);max-height:72vh;
    display:flex;flex-direction:column;box-shadow:0 20px 56px rgba(0,0,0,.28);overflow:hidden;
  }
  .nb-sel-header{
    display:flex;align-items:center;padding:16px 18px;
    border-bottom:1px solid var(--line);gap:10px;
  }
  .nb-sel-header h3{margin:0;font-size:15px;flex:1;}
  .nb-sel-header .nb-sel-close{
    border:none;background:none;font-size:22px;line-height:1;
    cursor:pointer;color:var(--muted);padding:0;
  }
  .nb-sel-search{
    margin:12px 16px 0;padding:8px 12px;border:1px solid var(--line);
    border-radius:7px;font-size:13px;font-family:Consolas,"SF Mono",Menlo,monospace;
  }
  .nb-sel-actions{display:flex;gap:8px;padding:8px 16px;align-items:center;}
  .nb-sel-actions button{
    font-size:12px;padding:5px 12px;border:1px solid var(--line);
    background:#fff;border-radius:6px;cursor:pointer;
  }
  .nb-sel-actions button:hover{background:var(--accent-soft);}
  .nb-sel-show{
    background:var(--accent)!important;color:#fff!important;
    border-color:var(--accent)!important;margin-left:auto!important;font-weight:600!important;
  }
  .nb-sel-show:hover{opacity:.88!important;}
  .nb-sel-item input[type=checkbox]{
    width:16px;height:16px;cursor:pointer;margin-left:auto;flex-shrink:0;
    accent-color:var(--accent);
  }
  .nb-sel-list{overflow-y:auto;padding:4px 8px 16px;}
  .nb-sel-group-title{
    font-family:Consolas,"SF Mono",Menlo,monospace;font-size:11px;
    color:var(--accent);margin:12px 8px 4px;letter-spacing:.05em;
  }
  .nb-sel-count{color:var(--muted);}
  .nb-sel-item{
    display:flex;align-items:center;gap:10px;padding:8px 10px;
    border-radius:7px;cursor:pointer;font-size:13px;
    font-family:Consolas,"SF Mono",Menlo,monospace;
  }
  .nb-sel-item:hover{background:var(--accent-soft);}
  .nb-sel-item.active{background:var(--selected);color:#7a3b12;}
  .nb-sel-item.hidden-row{display:none;}
  .nb-sel-icon{font-size:16px;width:20px;text-align:center;flex-shrink:0;}
  .nb-sel-name{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
  .nb-sel-sub{font-size:11px;color:var(--muted);}
  .nb-sel-all{font-weight:600;color:var(--accent);}

  /* ---- table cards (shared for all sortable/filterable tables) ---- */
  .tbl-card{background:var(--panel);border:1px solid var(--line);border-radius:10px;overflow:hidden;}
  .tbl-controls{
    display:flex;align-items:center;gap:8px;padding:10px 14px;
    background:#f8f7f4;border-bottom:1px solid var(--line);flex-wrap:wrap;
  }
  .tbl-filter{
    flex:1 1 200px;padding:6px 10px;border:1px solid var(--line);
    border-radius:6px;font-size:12.5px;font-family:Consolas,"SF Mono",Menlo,monospace;
  }
  .tbl-filter:focus{outline:2px solid var(--accent);outline-offset:1px;}
  .tbl-count{font-size:12px;color:var(--muted);white-space:nowrap;}
  .tbl-dl-wrap{display:flex;gap:2px;}
  .tbl-dl-btn{
    border:1px solid var(--accent);background:var(--accent);color:#fff;
    border-radius:0;padding:6px 10px;font-size:12px;cursor:pointer;white-space:nowrap;
  }
  .tbl-dl-btn:first-child{border-radius:6px 0 0 6px;}
  .tbl-dl-btn:last-child{border-radius:0 6px 6px 0;}
  .tbl-dl-btn:hover{opacity:.88;}
  /* "download" alias kept so Python helper still works */
  .tbl-download{display:none;}
  .tbl-scroll{max-height:420px;overflow:auto;}
  table{border-collapse:collapse;width:100%;
    font-size:12.5px;font-family:Consolas,"SF Mono",Menlo,monospace;}
  /* header row 1: column names (sticky row 0) */
  thead tr.tbl-head-labels th{
    position:sticky;top:0;z-index:3;background:var(--accent-soft);color:var(--ink);
    text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);
    white-space:nowrap;cursor:pointer;user-select:none;
  }
  thead tr.tbl-head-labels th:hover{background:#d5eaec;}
  thead tr.tbl-head-labels th.sort-asc::after{content:" ▲";font-size:10px;}
  thead tr.tbl-head-labels th.sort-desc::after{content:" ▼";font-size:10px;}
  /* header row 2: per-column filter inputs (sticky row 1) */
  thead tr.tbl-head-filters th{
    position:sticky;top:31px;z-index:3;background:#eef6f7;
    padding:3px 4px;border-bottom:2px solid var(--line);
  }
  .tbl-col-filter{
    width:100%;padding:3px 6px;font-size:11.5px;
    font-family:Consolas,"SF Mono",Menlo,monospace;
    border:1px solid var(--line);border-radius:4px;background:#fff;
  }
  .tbl-col-filter:focus{outline:1px solid var(--accent);}
  tbody td{padding:6px 10px;border-bottom:1px solid var(--line);white-space:nowrap;}
  tbody tr:nth-child(even){background:#fbfbf9;}
  tbody tr.hidden-row{display:none;}
  .table-note{color:var(--muted);font-size:12px;margin-top:8px;}
  /* nb-crumb shows count of visible panels */
  #nb-crumb-count{font-size:11px;background:var(--accent-soft);color:var(--accent);
    border-radius:10px;padding:1px 8px;font-weight:600;}

  /* ---- stats section ---- */
  .stats-toggle{display:flex;gap:0;border-bottom:2px solid var(--line);margin-bottom:18px;}
  .stats-btn{
    padding:9px 22px;font-size:13px;cursor:pointer;
    border:none;background:none;color:var(--muted);
    border-bottom:3px solid transparent;margin-bottom:-2px;
  }
  .stats-btn.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:600;}
  .stats-btn:hover{color:var(--ink);}
  .stats-block{display:none;}
  .stats-block.active{display:block;}
  .stats-note{color:var(--muted);font-size:13px;margin:0 0 16px;}
  .stats-chart{margin-bottom:20px;overflow-x:auto;}
  .stats-chart svg{display:block;max-width:100%;height:auto;}

  /* ---- column filter popup ---- */
  .cfp-backdrop{
    display:none;position:fixed;inset:0;z-index:800;
  }
  .cfp-backdrop.open{display:block;}
  .cfp-popup{
    position:fixed;z-index:801;background:#fff;border:1px solid var(--line);
    border-radius:10px;box-shadow:0 8px 32px rgba(0,0,0,.18);
    min-width:240px;max-width:300px;overflow:hidden;
    font-size:13px;font-family:Consolas,"SF Mono",Menlo,monospace;
  }
  .cfp-head{
    display:flex;align-items:center;padding:10px 14px 8px;
    border-bottom:1px solid var(--line);gap:8px;
  }
  .cfp-head span{font-weight:600;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
  .cfp-head button{border:none;background:none;cursor:pointer;color:var(--muted);font-size:18px;line-height:1;padding:0;}
  /* numeric mode */
  .cfp-num-ops{display:flex;flex-wrap:wrap;gap:5px;padding:10px 14px;}
  .cfp-op-btn{
    padding:4px 9px;border:1px solid var(--line);background:#fff;
    border-radius:20px;font-size:12px;cursor:pointer;font-family:inherit;
  }
  .cfp-op-btn.active{background:var(--accent);color:#fff;border-color:var(--accent);}
  .cfp-op-btn:hover:not(.active){background:var(--accent-soft);}
  .cfp-num-inputs{padding:0 14px 10px;display:flex;gap:6px;align-items:center;}
  .cfp-num-inputs input{
    flex:1;padding:6px 9px;border:1px solid var(--line);border-radius:6px;
    font-size:13px;font-family:inherit;
  }
  .cfp-num-inputs input:focus{outline:1px solid var(--accent);}
  .cfp-num-inputs .cfp-between-sep{color:var(--muted);font-size:11px;}
  /* text/categorical mode */
  .cfp-text-search{
    margin:10px 14px 6px;padding:6px 10px;border:1px solid var(--line);
    border-radius:6px;font-size:12.5px;font-family:inherit;display:block;width:calc(100% - 28px);
  }
  .cfp-text-search:focus{outline:1px solid var(--accent);}
  .cfp-val-actions{display:flex;gap:6px;padding:0 14px 6px;}
  .cfp-val-actions button{font-size:11px;padding:2px 8px;border:1px solid var(--line);background:#fff;border-radius:4px;cursor:pointer;}
  .cfp-val-actions button:hover{background:var(--accent-soft);}
  .cfp-val-list{max-height:180px;overflow-y:auto;padding:2px 6px 8px;}
  .cfp-val-item{
    display:flex;align-items:center;gap:8px;padding:4px 8px;
    border-radius:5px;cursor:pointer;
  }
  .cfp-val-item:hover{background:var(--accent-soft);}
  .cfp-val-item.hidden-row{display:none;}
  .cfp-val-item input{width:14px;height:14px;accent-color:var(--accent);cursor:pointer;}
  .cfp-val-item .cfp-val-text{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
  .cfp-val-item .cfp-val-n{font-size:11px;color:var(--muted);}
  /* footer */
  .cfp-footer{display:flex;gap:8px;padding:10px 14px;border-top:1px solid var(--line);background:#f8f7f4;}
  .cfp-apply{
    flex:1;padding:6px;background:var(--accent);color:#fff;border:none;
    border-radius:6px;cursor:pointer;font-size:12.5px;font-weight:600;
  }
  .cfp-apply:hover{opacity:.88;}
  .cfp-clear{
    padding:6px 12px;background:#fff;color:var(--muted);border:1px solid var(--line);
    border-radius:6px;cursor:pointer;font-size:12.5px;
  }
  .cfp-clear:hover{background:var(--accent-soft);}
  /* active filter indicator on th */
  thead tr.tbl-head-labels th.col-filtered::after{
    content:"●";color:var(--accent);font-size:9px;margin-left:4px;vertical-align:super;
  }
  /* filter icon button in th */
  .tbl-col-filter-btn{
    border:none;background:none;cursor:pointer;padding:0 0 0 4px;
    color:var(--muted);font-size:12px;opacity:.6;vertical-align:middle;
  }
  .tbl-col-filter-btn:hover,.col-filtered .tbl-col-filter-btn{opacity:1;color:var(--accent);}

  /* ---- footer ---- */
  footer{
    text-align:center;color:var(--muted);font-size:12px;
    padding:24px 0 48px;border-top:1px solid var(--line);margin-top:40px;
  }
  footer b{color:var(--accent);}

  @media(max-width:820px){
    .go-track{grid-template-columns:100px 1fr;}
    .nb-zoom-group{margin-left:0;}
    .top-meta{display:none;}
  }
</style>
</head>
<body>

<nav class="top-nav">
  <div class="top-brand"><div class="top-logo-wrap">$header_logo_html</div><span class="brand-name">S(H)ARP</span></div>
  <div class="top-tabs">
    <button type="button" class="top-tab active" data-page="overview">
      <span class="tab-icon">&#127757;</span> Overview
    </button>
    <button type="button" class="top-tab" data-page="neighborhoods">
      <span class="tab-icon">&#9654;</span> Neighborhoods
      <span class="tab-badge">$n_blocks</span>
    </button>
    <button type="button" class="top-tab" data-page="statistics">
      <span class="tab-icon">&#9638;</span> Statistics
    </button>
    <button type="button" class="top-tab" data-page="data">
      <span class="tab-icon">&#9776;</span> Data
      <span class="tab-badge">$n_genes</span>
    </button>
  </div>
  <div class="top-meta"><b>$n_genes</b> genes &middot; <b>$n_blocks</b> neighborhoods</div>
</nav>

<!-- ===================== 01 OVERVIEW ===================== -->
<div class="page-section active" data-page="overview">
<div class="sec-inner">
  <h1 class="sec-title">$title</h1>
  <p class="sec-desc">Genome-wide position of each neighborhood. One track per contig, scaled to its own length. Hover a marker for details; click to open that neighborhood.</p>
  <div class="panel">
$genome_overview
  </div>
</div>
</div>

<!-- ===================== 02 NEIGHBORHOODS ===================== -->
<div class="page-section" data-page="neighborhoods">
<div class="sec-inner">
  <h1 class="sec-title">Neighborhoods</h1>
  <p class="sec-desc">Use the <b>&#9776; Select</b> button to choose which neighborhoods are in view. The <b>Figure</b> tab shows the gene map; the <b>Table</b> tab shows the input rows for that block. Hover any gene arrow for its info window.</p>

  <div class="nb-toolbar">
    <button type="button" class="nb-icon-btn" id="nb-sel-open">
      <svg viewBox="0 0 24 24"><circle cx="12" cy="5" r="1.5"/><circle cx="12" cy="12" r="1.5"/><circle cx="12" cy="19" r="1.5"/><line x1="5" y1="5" x2="19" y2="5"/><line x1="5" y1="12" x2="19" y2="12"/><line x1="5" y1="19" x2="19" y2="19"/></svg>
      Select
    </button>

    <div class="nb-zoom-group">
      <button type="button" id="nb-zoom-out" title="Zoom out">&minus;</button>
      <span id="nb-zoom-val">100%</span>
      <button type="button" id="nb-zoom-in" title="Zoom in">+</button>
      <button type="button" id="nb-zoom-reset" title="Reset">&#8635;</button>
      <button type="button" id="nb-zoom-25" title="25%">25%</button>
      <button type="button" id="nb-zoom-50" title="50%">50%</button>
      <button type="button" id="nb-zoom-200" title="200%">2×</button>
      <button type="button" id="nb-zoom-fit" title="Fit width">&#8596;</button>
    </div>
  </div>

  <div class="nb-crumb" id="nb-crumb">
    <span id="nb-crumb-count"></span>
    <div class="nb-nav-arrows">
      <button type="button" id="nb-prev" title="Previous neighborhood">&#8249;</button>
      <button type="button" id="nb-next" title="Next neighborhood">&#8250;</button>
    </div>
  </div>

  <div id="nb-panels">
$nb_panels
  </div>

</div>
</div>

<!-- ===================== 03 STATISTICS ===================== -->
<div class="page-section" data-page="statistics">
<div class="sec-inner">
  <h1 class="sec-title">Domain statistics</h1>
  <p class="sec-desc">How often each domain appears. <b>Reference domains</b> counts one per neighborhood (the reference query gene's domain). <b>Architectures</b> counts every occurrence across all proteins, splitting &lsquo;+&rsquo;-joined strings.</p>
  <div class="panel">
$stats_html
  </div>
</div>
</div>

<!-- ===================== 04 DATA ===================== -->
<div class="page-section" data-page="data">
<div class="sec-inner">
  <h1 class="sec-title">Input data</h1>
  <p class="sec-desc">All rows used to build this report. Filter, sort, and download as CSV.</p>
$data_card
</div>
</div>

<!-- selector pop-up -->
<div class="nb-sel-backdrop" id="nb-sel-modal">
  <div class="nb-sel-modal">
    <div class="nb-sel-header">
      <h3>Select neighborhoods</h3>
      <button type="button" class="nb-sel-close" id="nb-sel-close">&times;</button>
    </div>
    <input type="text" class="nb-sel-search" id="nb-sel-search" placeholder="Search by id or domain&hellip;">
    <div class="nb-sel-actions">
      <button type="button" id="nb-sel-all">&#9745; All</button>
      <button type="button" id="nb-sel-none">&#9744; None</button>
      <button type="button" class="nb-sel-show" id="nb-sel-apply">Show selected</button>
    </div>
    <div class="nb-sel-list" id="nb-sel-list">
$nb_selector
    </div>
  </div>
</div>

<div class="go-tooltip" id="go-tooltip"></div>

<!-- column-filter popup (singleton, repositioned by JS) -->
<div class="cfp-backdrop" id="cfp-backdrop"></div>
<div class="cfp-popup" id="cfp-popup" style="display:none">
  <div class="cfp-head">
    <span id="cfp-col-name"></span>
    <button type="button" id="cfp-close-btn">&times;</button>
  </div>
  <!-- numeric mode -->
  <div id="cfp-num-mode" style="display:none">
    <div class="cfp-num-ops" id="cfp-num-ops"></div>
    <div class="cfp-num-inputs">
      <input type="number" id="cfp-num-val" placeholder="value">
      <span class="cfp-between-sep" id="cfp-between-sep" style="display:none">and</span>
      <input type="number" id="cfp-num-val2" placeholder="value" style="display:none">
    </div>
  </div>
  <!-- text/categorical mode -->
  <div id="cfp-text-mode" style="display:none">
    <input type="text" class="cfp-text-search" id="cfp-text-search" placeholder="Search values&hellip;">
    <div class="cfp-val-actions">
      <button type="button" id="cfp-val-all">All</button>
      <button type="button" id="cfp-val-none">None</button>
    </div>
    <div class="cfp-val-list" id="cfp-val-list"></div>
  </div>
  <div class="cfp-footer">
    <button type="button" class="cfp-apply" id="cfp-apply-btn">Apply</button>
    <button type="button" class="cfp-clear" id="cfp-clear-btn">Clear filter</button>
  </div>
</div>

<footer><div class="footer-logo-wrap">$footer_logo_html</div>Made by <b>S(H)ARP</b> &mdash; Biosynthetic Gene Cluster Analysis</footer>

<script>
(function () {
  // ── helpers ──────────────────────────────────────────────────────────
  function qs(sel, ctx) { return (ctx || document).querySelector(sel); }
  function qsa(sel, ctx) { return Array.prototype.slice.call((ctx || document).querySelectorAll(sel)); }
  function on(el, ev, fn) { if (el) el.addEventListener(ev, fn); }

  // ── top-level page tabs ───────────────────────────────────────────────
  var pageTabs = qsa('.top-tab');
  var pageSecs = qsa('.page-section');
  function showPage(name) {
    pageTabs.forEach(function(t){ t.classList.toggle('active', t.dataset.page === name); });
    pageSecs.forEach(function(s){ s.classList.toggle('active', s.dataset.page === name); });
  }
  pageTabs.forEach(function(t){ on(t, 'click', function(){ showPage(t.dataset.page); }); });

  // ── genome overview ───────────────────────────────────────────────────
  var goZoom = 1;
  var goTracks = qsa('.go-track');
  var tip = qs('#go-tooltip');
  var goZoomVal = qs('#go-zoom-val');

  function goLayout() {
    goTracks.forEach(function (tr) {
      var vp = tr.querySelector('.go-viewport');
      var inner = tr.querySelector('.go-inner');
      var base = tr._base || vp.clientWidth || 600;
      var w = base * goZoom;
      inner.style.width = w + 'px';
      qsa('.go-marker', inner).forEach(function (m) {
        var left = parseFloat(m.dataset.fs) * w;
        var ww = Math.max(5, (parseFloat(m.dataset.fe) - parseFloat(m.dataset.fs)) * w);
        m.style.left = left + 'px'; m.style.width = ww + 'px';
      });
    });
    if (goZoomVal) goZoomVal.textContent = (goZoom < 10 ? goZoom.toFixed(1) : Math.round(goZoom)) + 'x';
  }
  function goSetZoom(z) { goZoom = Math.min(500, Math.max(1, z)); goLayout(); }
  function initBases() { goTracks.forEach(function (tr) { tr._base = tr.querySelector('.go-viewport').clientWidth || 600; }); }

  on(qs('#go-zoom-in'),    'click', function(){ goSetZoom(goZoom * 1.6); });
  on(qs('#go-zoom-out'),   'click', function(){ goSetZoom(goZoom / 1.6); });
  on(qs('#go-zoom-reset'), 'click', function(){ goSetZoom(1); goTracks.forEach(function(tr){ tr.querySelector('.go-viewport').scrollLeft=0; }); });

  goTracks.forEach(function (tr) {
    var vp = tr.querySelector('.go-viewport');
    on(vp, 'wheel', function (e) {
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        var rect = vp.getBoundingClientRect();
        var base = tr._base || vp.clientWidth;
        var frac = (vp.scrollLeft + e.clientX - rect.left) / (base * goZoom);
        goSetZoom(goZoom * (e.deltaY < 0 ? 1.2 : 1/1.2));
        vp.scrollLeft = frac * base * goZoom - (e.clientX - rect.left);
      } else { e.preventDefault(); vp.scrollLeft += e.deltaY + e.deltaX; }
    }, { passive: false });
    var drag=false, sx=0, ss=0;
    on(vp, 'mousedown', function(e){ if(e.target.classList.contains('go-marker')) return; drag=true; sx=e.clientX; ss=vp.scrollLeft; vp.classList.add('grabbing'); });
    on(window, 'mousemove', function(e){ if(drag) vp.scrollLeft=ss-(e.clientX-sx); });
    on(window, 'mouseup', function(){ drag=false; vp.classList.remove('grabbing'); });
  });

  // ── shared tooltip ────────────────────────────────────────────────────
  function moveTip(e) {
    var x=e.clientX+14, y=e.clientY+14, r=tip.getBoundingClientRect();
    if (x+r.width>window.innerWidth)  x=e.clientX-r.width-14;
    if (y+r.height>window.innerHeight) y=e.clientY-r.height-14;
    tip.style.left=x+'px'; tip.style.top=y+'px';
  }
  function attachTip(el) {
    on(el, 'mouseenter', function(e){ tip.innerHTML=el.dataset.tip; tip.style.display='block'; moveTip(e); });
    on(el, 'mousemove',  moveTip);
    on(el, 'mouseleave', function(){ tip.style.display='none'; });
  }
  var goMarkers = qsa('.go-marker');
  goMarkers.forEach(function(m){ attachTip(m); on(m,'click',function(){ openNeighborhood(m.dataset.block); }); });
  qsa('.nb-gene').forEach(attachTip);

  // ── neighborhood panels (multi-select) ───────────────────────────────
  var nbPanels = qsa('.nb-panel');
  var selItems = qsa('.nb-sel-item');

  // Which slugs are currently shown (set by applySelection)
  var activeSet = {};

  // Per-panel sub-tab wiring (each panel has its own Figure/Table buttons)
  document.addEventListener('click', function(e) {
    var btn = e.target.closest && e.target.closest('.nb-panel-subtabs .nb-subtab');
    if (!btn) return;
    var panel = btn.closest('.nb-panel');
    if (!panel) return;
    qsa('.nb-subtab', panel).forEach(function(t){ t.classList.toggle('active', t===btn); });
    qsa('.nb-view',   panel).forEach(function(v){ v.classList.toggle('active', v.classList.contains(btn.dataset.view)); });
  });

  function updateCrumb() {
    var n = nbPanels.filter(function(p){ return p.classList.contains('active'); }).length;
    var el = qs('#nb-crumb-count');
    if (el) {
      el.innerHTML = n > 0 ?
        (n + '&nbsp;neighborhood' + (n>1?'s':'') + '&nbsp;shown') : '';
    }
  }

  function applySelection(slugs) {
    activeSet = {};
    slugs.forEach(function(s){ activeSet[s]=true; });
    nbPanels.forEach(function(p){ p.classList.toggle('active', !!activeSet[p.dataset.block]); });
    selItems.forEach(function(i){
      i.classList.toggle('active', !!activeSet[i.dataset.block]);
      var cb = i.querySelector('input[type=checkbox]');
      if (cb) cb.checked = !!activeSet[i.dataset.block];
    });
    goMarkers.forEach(function(m){ m.classList.toggle('selected', !!activeSet[m.dataset.block]); });
    updateCrumb();
  }

  function selectBlock(slug) {
    applySelection(slug ? [slug] : []);
  }
  window._selectBlock = selectBlock;

  function openNeighborhood(slug) {
    showPage('neighborhoods');
    selectBlock(slug);
  }

  // inject checkboxes into selector items
  selItems.forEach(function(i){
    var cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.addEventListener('click', function(e){ e.stopPropagation(); });
    on(cb, 'change', function(e){ e.stopPropagation(); });
    i.appendChild(cb);
  });

  // prev / next (cycle through *all* panels, show only that one)
  var allSlugs = nbPanels.map(function(p){ return p.dataset.block; });
  function navigate(dir) {
    var activeSlugs = Object.keys(activeSet);
    var ref = activeSlugs[0] || allSlugs[0];
    var idx = allSlugs.indexOf(ref);
    var next = allSlugs[(idx + dir + allSlugs.length) % allSlugs.length];
    if (next) selectBlock(next);
  }
  on(qs('#nb-prev'), 'click', function(){ navigate(-1); });
  on(qs('#nb-next'), 'click', function(){ navigate(1); });

  // ── neighborhood zoom ─────────────────────────────────────────────────
  var nbZoom = 1;
  function nbLayout() {
    qsa('.nb-fig').forEach(function(f){ f.style.width=(nbZoom*100)+'%'; });
    var el = qs('#nb-zoom-val'); if(el) el.textContent=Math.round(nbZoom*100)+'%';
  }
  function nbSetZoom(z){ nbZoom=Math.min(10,Math.max(.25,z)); nbLayout(); }
  on(qs('#nb-zoom-in'),    'click', function(){ nbSetZoom(nbZoom*1.3); });
  on(qs('#nb-zoom-out'),   'click', function(){ nbSetZoom(nbZoom/1.3); });
  on(qs('#nb-zoom-reset'), 'click', function(){ nbSetZoom(1); });
  on(qs('#nb-zoom-25'),    'click', function(){ nbSetZoom(.25); });
  on(qs('#nb-zoom-50'),    'click', function(){ nbSetZoom(.5); });
  on(qs('#nb-zoom-200'),   'click', function(){ nbSetZoom(2); });
  on(qs('#nb-zoom-fit'),   'click', function(){ nbSetZoom(1); });

  // ── selector modal ────────────────────────────────────────────────────
  var selModal = qs('#nb-sel-modal');
  function openSel() {
    // sync checkboxes with current active state
    selItems.forEach(function(i){
      var cb = i.querySelector('input[type=checkbox]');
      if (cb) cb.checked = !!activeSet[i.dataset.block];
    });
    if(selModal) selModal.classList.add('open');
    var s=qs('#nb-sel-search'); if(s){s.value=''; s.dispatchEvent(new Event('input')); s.focus();}
  }
  function closeSel() { if(selModal) selModal.classList.remove('open'); }
  on(qs('#nb-sel-open'),  'click', openSel);
  on(qs('#nb-sel-close'), 'click', closeSel);
  on(selModal, 'click', function(e){ if(e.target===selModal) closeSel(); });
  on(qs('#nb-sel-search'), 'input', function(){
    var term = this.value.toLowerCase();
    selItems.forEach(function(i){
      i.classList.toggle('hidden-row', i.textContent.toLowerCase().indexOf(term)===-1);
    });
  });
  on(qs('#nb-sel-all'), 'click', function(){
    selItems.forEach(function(i){
      if(!i.classList.contains('hidden-row')){
        var cb=i.querySelector('input[type=checkbox]'); if(cb) cb.checked=true;
      }
    });
  });
  on(qs('#nb-sel-none'), 'click', function(){
    selItems.forEach(function(i){
      var cb=i.querySelector('input[type=checkbox]'); if(cb) cb.checked=false;
    });
  });
  on(qs('#nb-sel-apply'), 'click', function(){
    var slugs=[];
    selItems.forEach(function(i){
      var cb=i.querySelector('input[type=checkbox]');
      if(cb && cb.checked) slugs.push(i.dataset.block);
    });
    if(slugs.length===0) slugs = allSlugs.slice(0,1);
    applySelection(slugs);
    closeSel();
  });
  on(window, 'keydown', function(e){
    if(e.key==='Escape') closeSel();
    if(e.key==='ArrowRight' && !e.target.matches('input,textarea')) navigate(1);
    if(e.key==='ArrowLeft'  && !e.target.matches('input,textarea')) navigate(-1);
  });

  // ── column-filter popup (singleton) ──────────────────────────────────
  var cfpBackdrop = qs('#cfp-backdrop');
  var cfpPopup    = qs('#cfp-popup');
  var cfpTarget   = null;  // { card, colIdx, headers, rows, colFilters, applyFn }

  function closeCfp() {
    if (cfpPopup) cfpPopup.style.display = 'none';
    if (cfpBackdrop) cfpBackdrop.classList.remove('open');
    cfpTarget = null;
  }
  on(cfpBackdrop, 'click', closeCfp);
  on(qs('#cfp-close-btn'), 'click', closeCfp);

  var NUM_OPS = [
    {op:'=',   label:'= equal'},
    {op:'!=',  label:'&#x2260; not'},
    {op:'>',   label:'&gt; greater'},
    {op:'>=',  label:'&ge; at least'},
    {op:'<',   label:'&lt; less'},
    {op:'<=',  label:'&le; at most'},
    {op:'[]',  label:'&#x2208; between'},
  ];

  function openCfp(card, colIdx, headers, rows, colFilters, applyFn, thEl) {
    // gather column values
    var vals = rows.map(function(r){ return r.cells[colIdx] ? r.cells[colIdx].textContent.trim() : ''; });
    var nonEmpty = vals.filter(function(v){ return v !== ''; });
    var isNum = nonEmpty.length > 0 && nonEmpty.every(function(v){ return v !== '' && !isNaN(parseFloat(v)); });

    // position popup below the th
    var rect = thEl.getBoundingClientRect();
    cfpPopup.style.display = 'block';
    cfpPopup.style.left = Math.min(rect.left, window.innerWidth - 310) + 'px';
    cfpPopup.style.top  = (rect.bottom + 4) + 'px';
    cfpBackdrop.classList.add('open');

    // column name
    var nameEl = qs('#cfp-col-name');
    if (nameEl) nameEl.textContent = headers[colIdx].textContent.replace(/[▲▼]/g,'').trim();

    var cur = colFilters[colIdx];

    if (isNum) {
      qs('#cfp-num-mode').style.display = '';
      qs('#cfp-text-mode').style.display = 'none';

      // build op buttons
      var opsEl = qs('#cfp-num-ops');
      opsEl.innerHTML = '';
      NUM_OPS.forEach(function(item){
        var btn = document.createElement('button');
        btn.className = 'cfp-op-btn';
        btn.innerHTML = item.label;
        btn.dataset.op = item.op;
        if (cur && cur.op === item.op) btn.classList.add('active');
        on(btn, 'click', function(){
          qsa('.cfp-op-btn', opsEl).forEach(function(b){ b.classList.remove('active'); });
          btn.classList.add('active');
          var isBetween = item.op === '[]';
          qs('#cfp-num-val2').style.display   = isBetween ? '' : 'none';
          qs('#cfp-between-sep').style.display = isBetween ? '' : 'none';
        });
        opsEl.appendChild(btn);
      });
      var v1 = qs('#cfp-num-val'), v2 = qs('#cfp-num-val2'), sep = qs('#cfp-between-sep');
      v1.value = (cur && cur.val  != null) ? cur.val  : '';
      v2.value = (cur && cur.val2 != null) ? cur.val2 : '';
      var isBetween = cur && cur.op === '[]';
      v2.style.display   = isBetween ? '' : 'none';
      sep.style.display  = isBetween ? '' : 'none';

      on(qs('#cfp-apply-btn'), 'click', function(){
        var activeOp = qs('.cfp-op-btn.active', opsEl);
        if (!activeOp || v1.value === '') { cfpTarget=null; closeCfp(); return; }
        colFilters[colIdx] = {type:'num', op: activeOp.dataset.op, val: parseFloat(v1.value), val2: parseFloat(v2.value||0)};
        headers[colIdx].classList.add('col-filtered');
        applyFn(); closeCfp();
      }, {once:true});
      on(qs('#cfp-clear-btn'), 'click', function(){
        colFilters[colIdx] = null;
        headers[colIdx].classList.remove('col-filtered');
        applyFn(); closeCfp();
      }, {once:true});

    } else {
      qs('#cfp-num-mode').style.display  = 'none';
      qs('#cfp-text-mode').style.display = '';

      // count unique values
      var freq = {};
      vals.forEach(function(v){ freq[v] = (freq[v]||0) + 1; });
      var uniq = Object.keys(freq).sort();
      var curSet = (cur && cur.type==='text' && cur.selected) ? cur.selected : null;

      var listEl = qs('#cfp-val-list');
      listEl.innerHTML = '';
      qs('#cfp-text-search').value = '';

      uniq.forEach(function(v){
        var item = document.createElement('div');
        item.className = 'cfp-val-item';
        var cb = document.createElement('input'); cb.type='checkbox';
        cb.checked = !curSet || !!curSet[v];
        var span = document.createElement('span'); span.className='cfp-val-text'; span.textContent=v||'(empty)';
        var cnt  = document.createElement('span'); cnt.className='cfp-val-n'; cnt.textContent=freq[v];
        item.appendChild(cb); item.appendChild(span); item.appendChild(cnt);
        on(item, 'click', function(e){ if(e.target!==cb) cb.checked=!cb.checked; });
        listEl.appendChild(item);
      });

      on(qs('#cfp-text-search'), 'input', function(){
        var term = this.value.toLowerCase();
        qsa('.cfp-val-item', listEl).forEach(function(it){
          it.classList.toggle('hidden-row', term && it.querySelector('.cfp-val-text').textContent.toLowerCase().indexOf(term)===-1);
        });
      });
      on(qs('#cfp-val-all'),  'click', function(){ qsa('.cfp-val-item:not(.hidden-row) input', listEl).forEach(function(c){ c.checked=true; }); });
      on(qs('#cfp-val-none'), 'click', function(){ qsa('.cfp-val-item:not(.hidden-row) input', listEl).forEach(function(c){ c.checked=false; }); });

      on(qs('#cfp-apply-btn'), 'click', function(){
        var sel = {};
        qsa('.cfp-val-item', listEl).forEach(function(it){
          var cb = it.querySelector('input'); var v = it.querySelector('.cfp-val-text').textContent;
          if(v==='(empty)') v='';
          if(cb && cb.checked) sel[v]=true;
        });
        var allChecked = Object.keys(sel).length === uniq.length;
        colFilters[colIdx] = allChecked ? null : {type:'text', selected: sel};
        if (allChecked) headers[colIdx].classList.remove('col-filtered');
        else headers[colIdx].classList.add('col-filtered');
        applyFn(); closeCfp();
      }, {once:true});
      on(qs('#cfp-clear-btn'), 'click', function(){
        colFilters[colIdx] = null;
        headers[colIdx].classList.remove('col-filtered');
        applyFn(); closeCfp();
      }, {once:true});
    }
  }

  // ── table cards ───────────────────────────────────────────────────────
  function initTableCard(card) {
    var globalFilter = card.querySelector('.tbl-filter');
    var countEl      = card.querySelector('.tbl-count');
    var table        = card.querySelector('table');
    if (!table) return;
    var thead  = table.querySelector('thead');
    var tbody  = table.querySelector('tbody');
    var rows   = qsa('tr', tbody);
    var labelRow = thead.querySelector('tr');
    labelRow.classList.add('tbl-head-labels');
    var headers = qsa('th', labelRow);
    var sortCol = -1, sortDir = 1;

    // per-column filter state: null = no filter; otherwise {type, ...}
    var colFilters = headers.map(function(){ return null; });

    // inject filter-icon button into each header th
    headers.forEach(function(h, i){
      var btn = document.createElement('button');
      btn.className = 'tbl-col-filter-btn';
      btn.title = 'Filter this column';
      btn.innerHTML = '&#9660;';
      btn.dataset.col = i;
      on(btn, 'click', function(e){
        e.stopPropagation();
        openCfp(card, i, headers, rows, colFilters, applyFilter, h);
      });
      h.appendChild(btn);
    });

    // ── filter logic ──────────────────────────────────────────────────
    function updateCount() {
      var vis = rows.filter(function(r){ return !r.classList.contains('hidden-row'); }).length;
      if (countEl) countEl.textContent = vis + ' / ' + rows.length + ' rows';
    }
    function testNum(val, f) {
      var n = parseFloat(val);
      if (isNaN(n)) return false;
      switch(f.op) {
        case '=':  return n === f.val;
        case '!=': return n !== f.val;
        case '>':  return n > f.val;
        case '>=': return n >= f.val;
        case '<':  return n < f.val;
        case '<=': return n <= f.val;
        case '[]': return n >= f.val && n <= f.val2;
        default:   return true;
      }
    }
    function applyFilter() {
      var global = globalFilter ? globalFilter.value.toLowerCase() : '';
      rows.forEach(function(r) {
        var cells = qsa('td', r);
        var failG = global && r.textContent.toLowerCase().indexOf(global) === -1;
        var failC = colFilters.some(function(f, i) {
          if (!f) return false;
          var v = cells[i] ? cells[i].textContent.trim() : '';
          if (f.type === 'num')  return !testNum(v, f);
          if (f.type === 'text') return !f.selected[v] && !(v==='' && f.selected['']);
          return false;
        });
        r.classList.toggle('hidden-row', failG || failC);
      });
      updateCount();
    }
    on(globalFilter, 'input', applyFilter);

    // ── sort ──────────────────────────────────────────────────────────
    function sortBy(col) {
      if (sortCol===col) sortDir=-sortDir; else { sortCol=col; sortDir=1; }
      headers.forEach(function(h,i){
        h.classList.remove('sort-asc','sort-desc');
        if(i===col) h.classList.add(sortDir===1?'sort-asc':'sort-desc');
      });
      rows.sort(function(a,b){
        var av=a.cells[col]?a.cells[col].textContent:'';
        var bv=b.cells[col]?b.cells[col].textContent:'';
        var an=parseFloat(av), bn=parseFloat(bv);
        return ((!isNaN(an)&&!isNaN(bn))?(an-bn):av.localeCompare(bv))*sortDir;
      });
      rows.forEach(function(r){ tbody.appendChild(r); });
    }
    headers.forEach(function(h,i){
      on(h,'click',function(e){
        if(e.target.classList.contains('tbl-col-filter-btn')) return;
        sortBy(i);
      });
    });

    // ── multi-format download ─────────────────────────────────────────
    function getColNames() {
      return headers.map(function(h){ return h.textContent.replace(/[▲▼⌄]/g,'').trim(); });
    }
    function visRows() { return rows.filter(function(r){ return !r.classList.contains('hidden-row'); }); }
    function getCells(r) { return qsa('td',r).map(function(c){ return c.textContent; }); }
    function dlCSV() {
      var c=getColNames().map(function(n){ return '"'+n.replace(/"/g,'""')+'"'; });
      var o=[c.join(',')];
      visRows().forEach(function(r){ o.push(getCells(r).map(function(v){ return '"'+v.replace(/"/g,'""')+'"'; }).join(',')); });
      return {txt:o.join('\n'),mime:'text/csv',ext:'csv'};
    }
    function dlTSV() {
      var o=[getColNames().join('\t')];
      visRows().forEach(function(r){ o.push(getCells(r).join('\t')); });
      return {txt:o.join('\n'),mime:'text/tab-separated-values',ext:'tsv'};
    }
    function dlJSON() {
      var c=getColNames();
      var o=visRows().map(function(r){ var obj={}; getCells(r).forEach(function(v,i){ obj[c[i]]=v; }); return obj; });
      return {txt:JSON.stringify(o,null,2),mime:'application/json',ext:'json'};
    }
    var fmtMap={csv:dlCSV,tsv:dlTSV,json:dlJSON};
    qsa('.tbl-dl-btn',card).forEach(function(btn){
      on(btn,'click',function(){
        var d=(fmtMap[btn.dataset.fmt]||dlCSV)();
        var base=(card.dataset.filename||'table.csv').replace(/\.[^.]+$$/,'');
        var a=document.createElement('a');
        a.href=URL.createObjectURL(new Blob([d.txt],{type:d.mime}));
        a.download=base+'.'+d.ext; a.click();
      });
    });
    updateCount();
  }
  qsa('.tbl-card').forEach(initTableCard);

  // ── statistics toggle ─────────────────────────────────────────────────
  qsa('.stats-btn').forEach(function(btn){
    on(btn,'click',function(){
      var key=btn.dataset.stats;
      qsa('.stats-btn').forEach(function(b){ b.classList.toggle('active', b.dataset.stats===key); });
      qsa('.stats-block').forEach(function(b){ b.classList.toggle('active', b.dataset.stats===key); });
    });
  });


  // ── init ──────────────────────────────────────────────────────────────
  initBases();
  goLayout();
  nbLayout();
  // show the first panel by default
  if (allSlugs.length > 0) applySelection([allSlugs[0]]);
  on(window,'resize', function(){ initBases(); goLayout(); });
})();
</script>
</body>
</html>
""")


def build_neighborhood_panels(extents, block_svgs, block_tables=None,
                               combined_svg=None, combined_table=None, default_view='all'):
    """
    Build the neighborhoods section's inner fragments for the pop-up
    (icon) selection model -- there is no scrollable tab strip.

    Each neighborhood is one hidden panel containing two sub-views: its
    figure and its input-rows table (as a sortable/downloadable table
    card). A single pair of Figure/Table sub-tabs at the top switches
    which sub-view shows for whichever panel is active; the pop-up (and
    the prev/next arrows) choose which panel is active.

    Parameters
    ----------
    extents : pandas.DataFrame
        Output of `compute_block_extents`; drives order and labels.
    block_svgs : dict[str, str]
        Slug -> annotated SVG markup.
    block_tables : dict[str, str] or None
        Slug -> table-card HTML for that block's input rows.
    combined_svg : str or None
        If given, an "All neighborhoods" panel (data-block `__all__`) is
        added first, showing every block at once.
    combined_table : str or None
        Table-card HTML for the full input table, shown in the "All"
        panel's table sub-view (this is what lets the standalone data
        section be dropped -- the full table lives here instead).
    default_view : str
        Which panel is active on load: 'all' (the combined panel, if
        present) or 'first' (the first block).

    Returns
    -------
    (str, str)
        (panels_html, selector_items_html).
    """
    if extents.empty:
        return '<div class="nb-empty">No blocks to show.</div>', ''

    block_tables = block_tables or {}
    nucleotides = list(dict.fromkeys(extents['nucleotide']))
    has_combined = combined_svg is not None
    selected = '__all__' if (default_view == 'all' and has_combined) else _slug(extents.iloc[0]['ID'])

    panels = []
    selector_items = []

    def panel(slug, label, head, svg, table_card):
        # Every panel starts hidden; JS activates the right set on load.
        # Each panel owns its own Figure / Table sub-tabs so multiple
        # panels can be shown simultaneously with independent view state.
        subtabs = (
            '<div class="nb-panel-subtabs">'
            '<button type="button" class="nb-subtab active" data-view="nb-view-figure">&#9654; Figure</button>'
            '<button type="button" class="nb-subtab" data-view="nb-view-table">&#9776; Table</button>'
            '</div>'
        )
        fig = f'<div class="nb-view nb-view-figure active"><div class="nb-fig-scroll"><div class="nb-fig">{svg}</div></div></div>'
        tab = f'<div class="nb-view nb-view-table">{table_card or "<p class=\'nb-empty\'>No table.</p>"}</div>'
        return (
            f'<div class="nb-panel" data-block="{slug}" data-label="{html.escape(str(label))}">'
            f'{head}{subtabs}{fig}{tab}</div>'
        )

    if has_combined:
        head = (
            f'<div class="nb-panel-head"><b>All neighborhoods</b> &middot; '
            f'{len(extents)} blocks, every block stacked</div>'
        )
        panels.append(panel('__all__', 'All neighborhoods', head, combined_svg, combined_table))
        selector_items.append(
            '<div class="nb-sel-item nb-sel-all" data-block="__all__" role="button" tabindex="0">'
            '<span class="nb-sel-icon">&#9776;</span>'
            '<span class="nb-sel-name">All neighborhoods</span>'
            f'<span class="nb-sel-sub">{len(extents)} blocks</span></div>'
        )

    for nucleotide in nucleotides:
        rows = extents[extents['nucleotide'] == nucleotide]
        selector_items.append(
            f'<div class="nb-sel-group-title">{html.escape(str(nucleotide))} '
            f'<span class="nb-sel-count">{len(rows)}</span></div>'
        )
        for _, block in rows.iterrows():
            slug = _slug(block['ID'])
            label = block['query_pid'] if block['query_pid'] is not None else block['ID']
            domain = block['query_domain'] if block['query_domain'] is not None else '-'
            head = (
                f'<div class="nb-panel-head"><b>{html.escape(str(label))}</b> &middot; '
                f'{html.escape(str(block["ID"]))} &middot; {html.escape(str(domain))} &middot; '
                f'{html.escape(str(block["org_name"]))} &middot; '
                f'{block["block_start"]:,.0f}&ndash;{block["block_end"]:,.0f} bp &middot; '
                f'{int(block["n_genes"])} genes</div>'
            )
            panels.append(panel(slug, label, head, block_svgs.get(slug, ''), block_tables.get(slug, '')))
            selector_items.append(
                f'<div class="nb-sel-item" data-block="{slug}" role="button" tabindex="0" '
                f'title="{html.escape(str(block["ID"]))}">'
                f'<span class="nb-sel-icon">&#9673;</span>'
                f'<span class="nb-sel-name">{html.escape(str(label))}</span>'
                f'<span class="nb-sel-sub">{int(block["n_genes"])} genes</span></div>'
            )

    return ''.join(panels), ''.join(selector_items)


def read_svg_logo(path):
    """
    Read an SVG file and prepare it for inline embedding in `build_html_report`.

    Strips the optional XML declaration and removes the top-level `width`/
    `height` attributes so the browser scales the SVG to fit whatever CSS
    container it is placed in (the container sizes are controlled by the
    `.top-logo-wrap` / `.footer-logo-wrap` CSS rules in the report template).

    Parameters
    ----------
    path : str
        Path to an SVG file.

    Returns
    -------
    str
        SVG markup ready to pass as `header_logo` or `footer_logo` to
        `build_html_report`.
    """
    import re as _re
    svg = open(path).read()
    svg = _re.sub(r'<\?xml[^?]*\?>\s*', '', svg)
    svg = _re.sub(r'\s+width="[^"]*"', '', svg, count=1)
    svg = _re.sub(r'\s+height="[^"]*"', '', svg, count=1)
    return svg.strip()


def build_html_report(df, output_file='operon_report.html', title='Gene Neighborhood Report',
                       group_col='block_id', org_col='organism', label_col='pfam',
                       rename_map=None, custom_colors=None, max_colors=5, ignore_domains=None,
                       nucleotide_col='nucleotide', start_col='start', end_col='end',
                       length_col='nlen', operon_kwargs=None, max_table_rows=2000,
                       work_dir=None, include_combined=True, default_view='all',
                       software_name='S(H)ARP',
                       header_logo=SHARP_HEADER_LOGO, footer_logo=SHARP_FOOTER_LOGO):
    """
    Build one self-contained, interactive HTML page for a single input
    table: a zoomable genome-wide overview with hover tooltips, a
    tabbed, sidebar-free neighborhoods section (one tab per result, plus
    an "All neighborhoods" tab), and the raw data table with a filter
    box.

    Neighborhoods section behavior:

      * each figure spans the full width of the panel and does NOT shrink
        with the number of neighbors; its own zoom controls stretch it;
      * a "Filter neighborhoods" pop-up lets the user mark/unmark which
        results appear as tabs;
      * selecting a tab shows that neighborhood's figure together with the
        slice of the input table used to draw it;
      * every protein has a hover info window (id, coordinates, strand,
        length, product; the query protein is flagged as the query in
        place of a domain line).

    A single shared domain -> color map is computed once (from the whole
    table, honoring `rename_map`/`custom_colors`/`max_colors`/
    `ignore_domains`) and reused for the genome overview and every
    per-block figure, so a given domain is the same color everywhere.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw input table; shown as-is in the data section and, per block,
        under each neighborhood figure.
    output_file : str
        Path to write the HTML report to.
    title : str
    group_col, org_col, label_col : str
        Block, organism and architecture columns (as elsewhere).
    rename_map, custom_colors, max_colors, ignore_domains :
        Color/label controls, identical in meaning to
        `neighborhood_figure`; applied once, globally.
    nucleotide_col, start_col, end_col, length_col : str
        Genomic-coordinate columns for the overview (see
        `compute_block_extents`).
    operon_kwargs : dict or None
        Extra per-figure options forwarded to `neighborhood_figure` for
        each block (e.g. `collapse_opposite_strand=True`, `font_size`,
        or `align_query_center=True` to re-enable centering, which is
        off by default). Color/label and `group_col`/`org_col`/
        `label_col` are handled centrally here, so don't pass those.
    max_table_rows : int or None
        Row cap for the main embedded table; `None` embeds every row.
        (Per-neighborhood sub-tables are never capped.)
    work_dir : str or None
        Where intermediate SVGs are written. If None, a temporary
        directory is used and cleaned up afterwards.
    include_combined : bool, default True
        Add the "All neighborhoods" tab (every block stacked in one
        figure). Turning this off saves time/size on very large inputs.
    default_view : str, default 'all'
        Which neighborhoods tab is open on load: 'all' (the combined
        tab) or 'first' (the first block). Ignored toward 'first' if
        `include_combined` is False.
    software_name : str, default 'S(H)ARP'
        Name shown in the report footer ("Made by ...").
    header_logo : str or None
        SVG string (already prepared for inline embedding) to show in the
        top-nav brand slot. Defaults to `SHARP_HEADER_LOGO` (bundled).
        Use `read_svg_logo(path)` to load a custom SVG from disk, or
        pass `None` to leave the slot empty.
    footer_logo : str or None
        SVG string for the footer logo. Defaults to `SHARP_FOOTER_LOGO`.
        Same conventions as `header_logo`.

    Returns
    -------
    str
        `output_file` (the path written).
    """
    operon_kwargs = dict(operon_kwargs) if operon_kwargs else {}

    # One global color map, shared by every figure for consistency.
    working = prepare_dataframe(df, group_col=group_col, org_col=org_col,
                                label_col=label_col, rename_map=rename_map)
    color_map = build_color_map(working, max_colors=max_colors,
                                ignore_domains=ignore_domains, custom_colors=custom_colors)
    extents = compute_block_extents(
        working, nucleotide_col=nucleotide_col, start_col=start_col,
        end_col=end_col, length_col=length_col,
    )

    genome_overview = build_genome_overview_interactive_html(extents, color_map=color_map)

    # Per-block table cards (figure sub-tab in each neighborhood panel)
    block_tables = {
        _slug(block_id): render_table_card(
            block_df, filename=f'{_slug(block_id)}.csv',
        )
        for block_id, block_df in df.groupby(group_col, sort=False)
    }

    own_tmp = work_dir is None
    tmp_dir = work_dir or tempfile.mkdtemp(prefix='operon_report_')
    try:
        per_block_operon_kwargs = dict(operon_kwargs)
        per_block_operon_kwargs.update(
            org_col=org_col, label_col=label_col, rename_map=rename_map,
        )
        block_svgs = render_neighborhood_svgs_by_block(
            df, group_col=group_col, color_map=color_map,
            operon_kwargs=per_block_operon_kwargs, tmp_dir=tmp_dir,
        )

        combined_svg = None
        if include_combined:
            combined_path = os.path.join(tmp_dir, 'combined.svg')
            _working, combined_meta = neighborhood_figure(
                df, group_col=group_col, output_file=combined_path,
                color_map=color_map, collect_node_meta=True, **per_block_operon_kwargs,
            )
            with open(combined_path) as f:
                combined_svg = annotate_neighborhood_svg(_strip_svg_prolog(f.read()), combined_meta)
    finally:
        if own_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # Combined table card lives in the "All neighborhoods" table sub-view
    combined_table = render_table_card(df, filename='all_neighborhoods.csv') if include_combined else None

    nb_panels, nb_selector = build_neighborhood_panels(
        extents, block_svgs, block_tables=block_tables,
        combined_svg=combined_svg, combined_table=combined_table,
        default_view=default_view,
    )

    # Statistics
    ref_counts, arch_counts = compute_domain_stats(working, extents, ignore_domains=ignore_domains)
    stats_html = build_stats_section_html(ref_counts, arch_counts, color_map=color_map)

    # Main data card (Data tab — full table with sort/filter/CSV)
    data_card = render_table_card(df, filename='full_data.csv', max_rows=max_table_rows)

    n_blocks = df[group_col].nunique() if group_col in df.columns else 'NA'
    html_doc = HTML_REPORT_TEMPLATE.substitute(
        title=html.escape(title),
        n_genes=f'{len(df):,}',
        n_blocks=n_blocks,
        genome_overview=genome_overview,
        nb_panels=nb_panels,
        nb_selector=nb_selector,
        stats_html=stats_html,
        data_card=data_card,
        software_name=html.escape(software_name),
        header_logo_html=header_logo or '',
        footer_logo_html=footer_logo or '',
    )

    with open(output_file, 'w') as f:
        f.write(html_doc)

    return output_file