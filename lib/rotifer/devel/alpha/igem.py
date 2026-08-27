"""
operon_fig.py
=============

Draw gene-neighborhood ("operon") figures from a long-format table of
genes/proteins that have already been grouped into genomic blocks (e.g.
one block per hit of interest plus its surrounding genomic context).

Each block is drawn as one row of the figure:

    [ text label ]  [gene] [gene] [QUERY] [gene] [gene] ...

The text label is a 3-line box with the reference query protein id, the
block id and the organism name. Every gene is drawn as an arrow pointing
the way its strand says, after the whole block has been mirrored (when
needed) so the reference query reads left-to-right -- so neighbors
transcribed the other way point back at it. Genes are colored by a
domain/annotation label (Pfam, Aravind, or free-text product, depending
on `label_col`).

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
    build_scaled_block_svg        one block drawn to real genomic scale
    render_scaled_svgs_by_block   one to-scale SVG per block
    build_gene_tooltip_html       per-protein hover "info window" body
    annotate_neighborhood_svg     inject those tooltips into a graphviz SVG
    build_neighborhood_panels     pop-up selector + merged figure/to-scale stacks
    render_table_card             sortable/filterable/downloadable table widget
    render_neighborhood_table_card single merged, block-tagged table (all blocks)
    compute_domain_stats          reference-query + full-architecture domain counts
    build_bar_list_html           interactive HTML/CSS bar list (all domains, filterable/sortable)
    build_stats_section_html      granularity x scope toggleable statistics section HTML
    build_html_report             combine everything into one HTML page

How the drawing actually works (no GUI toolkit involved)
--------------------------------------------------------
There is no Qt, matplotlib, or other plotting/GUI library here. Two
very different rendering paths are used:

  * The neighborhood figures (the gene-arrow rows) are laid out by
    **Graphviz** -- the same C graph-layout engine behind `dot` -- which
    we drive from Python through the **pygraphviz** binding. Each gene
    is a Graphviz node (`shape=cds`/`triangle`), each row is
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
import math
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

# Default header branding logo: the SVG file loaded into the report's
# top-nav brand slot. `build_html_report` resolves it through
# `resolve_logo`, so `header_logo` accepts either a path (with `~`
# expanded) or ready-to-embed SVG markup.
SHARP_HEADER_LOGO_PATH = '~/projects/igem/2026/data/logo.svg'



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


# Graphviz shape used for gene arrows: `cds` is the synthetic-biology
# CDS pentagon -- a rectangle whose leading end tapers to a point. It
# points right by default and `orientation=180` turns it around, which
# is how minus-strand genes are drawn.
GENE_ARROW_SHAPE = 'cds'
GENE_ARROW_LEFT_ORIENTATION = 180

# Graphviz 'orientation' (degrees) that makes a `shape=triangle` node
# point left. Collapsed neighbors are by definition on the strand
# opposite the query's, and blocks are normalized so the query reads
# left-to-right, so their triangle always points left. (Verified
# empirically by rendering test shapes with Graphviz 2.43: a plain
# triangle points up by default, and `orientation=90` rotates it to
# point left.)
COLLAPSED_TRIANGLE_ORIENTATION = 90


# ---------------------------------------------------------------------------
# Per-gene node styling
# ---------------------------------------------------------------------------

def gene_node_style(row, query_canonical_strand, color_map, highlight_query=True,
                     collapse_opposite_strand=False, font_size=10):
    """
    Decide the Graphviz node attributes (shape/color/label/size) for one
    gene.

    Every gene is drawn as a `cds` arrow (a rectangle tapering to a
    point at its leading end) whose direction follows its strand: right
    for +1, left for -1, and a plain box for anything else. Since
    `normalize_block_strand` has usually already mirrored the block onto
    the reference query's strand, a row reads left-to-right along the
    query's transcription direction, and genes transcribed the other way
    visibly point back at it. Every gene is filled per `color_map` and
    labeled with its domain/annotation; the query additionally gets a
    red, thicker outline when `highlight_query` is True.

    When `collapse_opposite_strand` is True, neighbors on the strand
    *opposite* the query's (`query_canonical_strand`) are drawn instead
    as small, unlabeled, grey triangles pointing left (back against the
    query's direction) -- a lightweight
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

    # Every gene -- query or neighbor -- points the way its strand says.
    # After `normalize_block_strand` that strand is the block-normalized
    # one, so a row reads left-to-right along the query's transcription
    # direction and opposite-strand neighbors point back at it.
    if strand_val == 1:
        shape_attrs = dict(shape=GENE_ARROW_SHAPE)
    elif strand_val == -1:
        shape_attrs = dict(shape=GENE_ARROW_SHAPE,
                           orientation=GENE_ARROW_LEFT_ORIENTATION)
    else:
        shape_attrs = dict(shape='box')

    return dict(
        **shape_attrs,
        label=str(row['domain']),
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
                        font_size=10, left_pad=0, right_pad=0, spacer_width=0.6,
                        show_row_label=True):
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
    show_row_label : bool, default True
        If False, the left-hand 3-line label (query protein id / block id
        / organism) is not drawn -- an invisible zero-size placeholder
        takes its place instead, so the alignment machinery (which
        references `label_node`) keeps working unchanged. Useful when
        something else already identifies the block (e.g. the report's
        neighborhood selector), so repeating it inside every stacked
        figure would just be noise.

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
    if show_row_label:
        label_html = build_row_label_html(
            query_pid=query_pid,
            block_id=block_df['ID'].iloc[0],
            org_name=block_df['org_name'].iloc[0],
            width=label_width,
            font_size=font_size,
        )
        graph.add_node(label_node_id, label=label_html, shape='none', margin=0.1)
    else:
        graph.add_node(label_node_id, label='', shape='point', style='invis', width=0.01, height=0.01)

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
                        spacer_width=0.6, color_map=None, collect_node_meta=False,
                        show_row_label=True):
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
    Every gene is drawn pointing the way its (block-normalized) strand
    says -- see `gene_node_style`. A block with more than one query gene uses the
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
    show_row_label : bool, default True
        If False, don't draw the left-hand 3-line label (query protein id
        / block id / organism) on any row. Off by default in the HTML
        report's per-block figures (the neighborhood selector already
        shows that info), so stacking several selected figures doesn't
        repeat it before every one. Left on by default here since
        `neighborhood_figure` is also used standalone.
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
            show_row_label=show_row_label,
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

    # NOTE: the shared tooltip div lives once at <body> level in the report
    # template (outside every page-section), not here -- a page-section is
    # display:none when inactive, which would hide an embedded tooltip too.
    return f'<div class="go-wrap">{controls}{"".join(rows)}</div>'


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


# ---------------------------------------------------------------------------
# To-scale ("biological scale") neighborhood view
# ---------------------------------------------------------------------------

def _nice_tick_step(span, target_ticks=6):
    """
    Pick a round tick interval (1/2/5 x 10^n) that puts roughly
    `target_ticks` ticks across `span` base pairs.
    """
    if span <= 0:
        return 1
    raw = span / max(target_ticks, 1)
    magnitude = 10 ** math.floor(math.log10(raw))
    for multiple in (1, 2, 5):
        if raw <= multiple * magnitude:
            return multiple * magnitude
    return 10 * magnitude


def _format_bp_tick(value, step):
    """
    Format a genomic coordinate for an axis tick: plain bp for small
    steps, kb once the ticks are 1 kb or more apart.
    """
    if step >= 1000:
        text = f'{value / 1000:,.1f}'.rstrip('0').rstrip('.')
        return f'{text} kb'
    return f'{value:,.0f}'


def build_scaled_block_svg(block_df, color_map=None, nucleotide_col='nucleotide',
                            start_col='start', end_col='end',
                            normalize_orientation=True, highlight_query=True,
                            track_width=900, left_margin=210, right_margin=30,
                            gene_height=26, font_size=11, min_gene_width=2.0,
                            show_row_label=True):
    """
    Draw one block to *biological scale*: genes placed by their real
    genomic coordinates, so arrow widths are proportional to gene
    lengths and the gaps between arrows are the real intergenic
    distances.

    This is the counterpart to `neighborhood_figure`, which lays genes
    out by Graphviz in even, text-sized boxes -- great for reading
    domain labels across rows, but it says nothing about how long a
    gene is or how far apart two genes sit. Here the x axis *is* the
    contig, in base pairs, with a labeled ruler underneath.

    Genes keep the same fill colors (`color_map`) and the same red query
    outline as the Graphviz figure, and each one is wrapped in the same
    `class="node nb-gene" data-tip="..."` group the report's JavaScript
    uses for its per-protein info window -- so hovering and click-to-pin
    work here exactly as they do in the Figure view.

    Orientation follows the same rule as `normalize_block_strand`: when
    `normalize_orientation` is True and the block's reference query is
    on the minus strand, the whole block is drawn reverse-complemented
    (coordinates mirrored, arrows flipped) so the query reads
    left-to-right. The ruler still shows real coordinates -- they simply
    count down from left to right -- and the header says so.

    Parameters
    ----------
    block_df : pandas.DataFrame
        One block, already through `prepare_dataframe` (needs 'pid',
        'domain', 'is_query', 'strand', plus `start_col`/`end_col`).
    color_map : dict[str, str] or None
        Domain -> fill color, as everywhere else (see `build_color_map`).
    nucleotide_col, start_col, end_col : str
        Per-gene contig and genomic span columns.
    normalize_orientation : bool, default True
        Mirror the block so its reference query points right.
    highlight_query : bool, default True
        Outline query genes in red.
    track_width, left_margin, right_margin, gene_height, font_size : float
        Layout, in SVG user units (effectively pixels at 100% zoom).
    min_gene_width : float
        Floor on how narrow an arrow may get, so a very short gene in a
        very wide block stays visible (and hoverable).
    show_row_label : bool, default True
        Draw the left-hand query-id / block-id / organism label column.

    Returns
    -------
    str
        A self-contained `<svg>...</svg>` fragment.
    """
    color_map = color_map or {}
    block = block_df.reset_index(drop=True)

    starts = pd.to_numeric(block.get(start_col), errors='coerce') if start_col in block.columns else None
    ends = pd.to_numeric(block.get(end_col), errors='coerce') if end_col in block.columns else None
    if starts is None or ends is None:
        spans = None
    else:
        lows = np.fmin(starts, ends)
        highs = np.fmax(starts, ends)
        spans = [(lo, hi) for lo, hi in zip(lows, highs)]

    valid = [s for s in (spans or []) if not (pd.isna(s[0]) or pd.isna(s[1]))]
    if not valid:
        return ('<svg viewBox="0 0 420 40" xmlns="http://www.w3.org/2000/svg" '
                'font-family="Consolas, \'SF Mono\', Menlo, monospace" font-size="11">'
                '<text x="8" y="24" fill="#888">No genomic coordinates for this block.</text></svg>')

    lo = min(s[0] for s in valid)
    hi = max(s[1] for s in valid)
    span = max(hi - lo, 1)

    ref_idx = select_reference_query_index(block)
    flip = bool(
        normalize_orientation and ref_idx is not None
        and block.loc[ref_idx, 'strand'] == -1
    )

    if not show_row_label:
        left_margin = 20

    def to_x(pos):
        frac = (hi - pos) / span if flip else (pos - lo) / span
        return left_margin + frac * track_width

    header_y = 16
    band_top = 30
    band_bottom = band_top + gene_height
    band_mid = (band_top + band_bottom) / 2
    axis_y = band_bottom + 22
    fig_width = left_margin + track_width + right_margin
    fig_height = axis_y + 26

    parts = [
        f'<svg viewBox="0 0 {fig_width:.0f} {fig_height:.0f}" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="Consolas, \'SF Mono\', Menlo, monospace" font-size="{font_size}">',
        f'<rect x="0" y="0" width="{fig_width:.0f}" height="{fig_height:.0f}" fill="white"/>',
    ]

    # header: contig, span, and whether we mirrored the block
    nucleotide = block[nucleotide_col].iloc[0] if nucleotide_col in block.columns else ''
    header = f'{nucleotide}  {lo:,.0f}-{hi:,.0f}  ({span:,.0f} bp)'
    if flip:
        header += '  · reverse-complemented'
    parts.append(
        f'<text x="{left_margin:.0f}" y="{header_y}" fill="#888" font-size="{font_size - 1}">'
        f'{html.escape(header)}</text>'
    )

    # left-hand label column, same three lines as the Graphviz figure
    if show_row_label:
        query_pid = block.loc[ref_idx, 'pid'] if ref_idx is not None else ''
        label_lines = [
            (str(query_pid), 'bold', '#111'),
            (str(block['ID'].iloc[0]) if 'ID' in block.columns else '', 'normal', '#333'),
            (str(block['org_name'].iloc[0]) if 'org_name' in block.columns else '', 'italic', '#555'),
        ]
        for i, (text, style, fill) in enumerate(label_lines):
            style_attr = ' font-weight="700"' if style == 'bold' else (
                ' font-style="italic"' if style == 'italic' else '')
            parts.append(
                f'<text x="8" y="{band_top - 2 + i * (font_size + 3):.0f}" fill="{fill}"'
                f'{style_attr}>{html.escape(text)}</text>'
            )

    # the contig line every gene sits on
    parts.append(
        f'<line x1="{left_margin:.1f}" y1="{band_mid:.1f}" '
        f'x2="{left_margin + track_width:.1f}" y2="{band_mid:.1f}" '
        f'stroke="#d0d0d0" stroke-width="1.5"/>'
    )

    # ---- genes ----
    for i, row in block.iterrows():
        g_lo, g_hi = spans[i]
        if pd.isna(g_lo) or pd.isna(g_hi):
            continue
        x_a, x_b = sorted((to_x(g_lo), to_x(g_hi)))
        width = max(x_b - x_a, min_gene_width)
        x_b = x_a + width

        strand_val = row.get('strand', 1)
        drawn_strand = -strand_val if (flip and strand_val in (1, -1)) else strand_val
        head = min(9.0, width * 0.45)
        if drawn_strand == -1:
            pts = [(x_b, band_top), (x_a + head, band_top), (x_a, band_mid),
                   (x_a + head, band_bottom), (x_b, band_bottom)]
        elif drawn_strand == 1:
            pts = [(x_a, band_top), (x_b - head, band_top), (x_b, band_mid),
                   (x_b - head, band_bottom), (x_a, band_bottom)]
        else:
            pts = [(x_a, band_top), (x_b, band_top), (x_b, band_bottom), (x_a, band_bottom)]
        points = ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts)

        is_target = bool(row['is_query'])
        stroke = 'red' if (highlight_query and is_target) else '#333'
        stroke_width = '2.4' if (highlight_query and is_target) else '1'
        fill = color_map.get(row.get('domain'), '#ffffff')

        meta = dict(
            pid=row.get('pid'),
            start=row.get(start_col),
            end=row.get(end_col),
            strand=strand_val,
            domain=row.get('domain'),
            product=row.get('product'),
            plen=row.get('plen'),
            is_query=is_target,
        )
        tip = html.escape(build_gene_tooltip_html(meta), quote=True)

        parts.append(f'<g class="node nb-gene" data-tip="{tip}">')
        parts.append(
            f'<polygon points="{points}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{stroke_width}"/>'
        )
        # the domain label only fits inside wide-enough arrows; the rest
        # rely on the hover/click info window
        label = str(row.get('domain', ''))
        if label and width >= len(label) * font_size * 0.62 + 8:
            parts.append(
                f'<text x="{(x_a + x_b) / 2:.1f}" y="{band_mid + font_size / 3:.1f}" '
                f'text-anchor="middle" fill="#111" pointer-events="none">'
                f'{html.escape(label)}</text>'
            )
        parts.append('</g>')

    # ---- ruler ----
    parts.append(
        f'<line x1="{left_margin:.1f}" y1="{axis_y:.1f}" '
        f'x2="{left_margin + track_width:.1f}" y2="{axis_y:.1f}" '
        f'stroke="#999" stroke-width="1"/>'
    )
    step = _nice_tick_step(span)
    first_tick = math.ceil(lo / step) * step
    tick = first_tick
    while tick <= hi:
        x = to_x(tick)
        parts.append(
            f'<line x1="{x:.1f}" y1="{axis_y:.1f}" x2="{x:.1f}" y2="{axis_y + 5:.1f}" '
            f'stroke="#999" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{axis_y + 17:.1f}" text-anchor="middle" fill="#777" '
            f'font-size="{font_size - 1}">{_format_bp_tick(tick, step)}</text>'
        )
        tick += step

    parts.append('</svg>')
    return '\n'.join(parts)


def render_scaled_svgs_by_block(working, color_map=None, nucleotide_col='nucleotide',
                                 start_col='start', end_col='end', **kwargs):
    """
    Render one to-scale figure per block (see `build_scaled_block_svg`).

    Parameters
    ----------
    working : pandas.DataFrame
        The prepared table from `prepare_dataframe` (blocks are taken
        from its 'pid_order'/'ID' columns, so block order matches
        `compute_block_extents`).
    color_map : dict[str, str] or None
        Shared domain -> color mapping, so a domain is the same color
        here, in the Graphviz figures and in the genome overview.
    nucleotide_col, start_col, end_col : str
        Coordinate columns, passed through.
    **kwargs
        Any other `build_scaled_block_svg` option.

    Returns
    -------
    dict[str, str]
        Block slug (see `_slug`) -> SVG markup.
    """
    svgs = {}
    for _, block_df in working.groupby('pid_order', sort=True):
        slug = _slug(block_df['ID'].iloc[0])
        svgs[slug] = build_scaled_block_svg(
            block_df, color_map=color_map, nucleotide_col=nucleotide_col,
            start_col=start_col, end_col=end_col, **kwargs,
        )
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
        '<div class="tbl-dl-wrap">''<button type="button" class="tbl-dl-btn" data-fmt="csv">&#8681; CSV</button>''<button type="button" class="tbl-dl-btn" data-fmt="tsv">TSV</button>''<button type="button" class="tbl-dl-btn" data-fmt="json">JSON</button>''</div>'
        '</div>'
        f'<div class="tbl-scroll">{table}</div>'
        '</div>'
    )


def compute_domain_stats(working, scope='all', ignore_domains=None):
    """
    Count how often domains appear, two ways, over a chosen subset of
    proteins:

      * "domains" -- unique ATOMIC domains: every '+'-joined architecture
        is split into its individual components, so `GntR+FCD` contributes
        one count each to `GntR` and to `FCD`. This answers "how often does
        this single domain occur".
      * "architectures" -- the full COMPOSITE domain string counted as one
        unit, so `GntR+FCD` is its own category distinct from `GntR` or
        `FCD` alone. This answers "how often does this exact domain
        combination occur".

    Both are computed over whichever `scope` is requested:

      * 'all'       -- every protein (query genes and neighbors alike).
      * 'query'     -- only each block's reference query gene. Since every
        block contributes exactly one query, this is the count "how many
        neighborhoods are built around this domain" -- and keeping it
        separate matters, because query genes are often deliberately
        similar (that's how the neighborhoods were assembled), so mixing
        them into 'all' can make a domain look far more common in the
        data than it actually is among the neighbors.
      * 'neighbors' -- every protein EXCEPT the query genes -- what
        actually surrounds the queries, with the query's own domain
        removed from the count entirely.

    Generic/uninformative values (`ignore_domains`, plus 'unk'/'-'/'?' and
    anything containing "hypothetical") are dropped. For architectures a
    component is dropped from the composite; if nothing meaningful remains,
    the whole architecture is skipped.

    Parameters
    ----------
    working : pandas.DataFrame
        Prepared table (needs the 'domain' and 'is_query' columns); see
        `prepare_dataframe`.
    scope : str
        'all', 'query', or 'neighbors' (see above).
    ignore_domains : list[str] or None
        Extra values to drop. Defaults to `DEFAULT_IGNORE_DOMAINS`.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        (domain_counts, architecture_counts), each with columns
        'domain' and 'count', sorted by count descending. (The first is
        atomic single domains; the second is composite architectures.)
    """
    ignore = {d.lower() for d in (ignore_domains or DEFAULT_IGNORE_DOMAINS)} | {'unk', '-', '?', ''}

    if scope == 'query' and 'is_query' in working.columns:
        subset = working[working['is_query']]
    elif scope == 'neighbors' and 'is_query' in working.columns:
        subset = working[~working['is_query']]
    else:
        subset = working

    def keep(value):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return False
        s = str(value).lower()
        return s not in ignore and 'hypothetical' not in s

    atomic = []        # individual domains (split on '+')
    architectures = []  # full composite strings, kept intact
    for value in subset['domain'].tolist():
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        parts = [p for p in str(value).split('+') if keep(p)]
        if not parts:
            continue
        atomic.extend(parts)
        # rebuild the architecture from the kept parts so ignored
        # components (e.g. signal peptides) don't fragment the label
        architectures.append('+'.join(parts))

    domain_counts = (
        pd.Series(atomic, dtype=object).value_counts()
        .rename_axis('domain').reset_index(name='count')
    )
    arch_counts = (
        pd.Series(architectures, dtype=object).value_counts()
        .rename_axis('domain').reset_index(name='count')
    )

    return domain_counts, arch_counts



def build_bar_list_html(counts_df, color_map=None, marker_color='#2a6f77'):
    """
    Render an interactive HTML/CSS horizontal bar list for a
    `{domain, count}` table, colored to match the figures.

    Unlike a static image, this shows EVERY row (not just a top-N slice)
    inside a scrollable container, and each row carries a `data-domain`
    attribute so the report's JS can filter (search box) and re-sort
    (by count or alphabetically) it client-side without re-rendering --
    hidden rows simply drop out of the flex flow, so the remaining bars
    reflow with no gaps.

    Parameters
    ----------
    counts_df : pandas.DataFrame
        Columns 'domain' and 'count', as returned by `compute_domain_stats`.
    color_map : dict[str, str] or None
        Domain -> color, to match the neighborhood figures. Domains not
        in the map fall back to `marker_color`.
    marker_color : str
        Fallback bar color.

    Returns
    -------
    str
        HTML for one `.bar-list` (or a "nothing to show" message).
    """
    color_map = color_map or {}
    if counts_df.empty:
        return '<p class="nb-empty">No domains to summarize.</p>'

    max_count = int(counts_df['count'].max()) if len(counts_df) else 0
    rows = []
    for _, row in counts_df.iterrows():
        domain = str(row['domain'])
        count = int(row['count'])
        pct = (count / max_count * 100) if max_count else 0
        fill = color_map.get(domain, marker_color)
        rows.append(
            '<div class="bar-row" data-domain="' + html.escape(domain.lower()) + '" '
            'data-count="' + str(count) + '" data-name="' + html.escape(domain) + '">'
            '<span class="bar-label" title="' + html.escape(domain) + '">' + html.escape(domain) + '</span>'
            '<span class="bar-track"><span class="bar-fill" style="width:' + f'{pct:.1f}' + '%;'
            'background:' + fill + ';"></span></span>'
            '<span class="bar-count">' + str(count) + '</span>'
            '</div>'
        )
    return '<div class="bar-list">' + ''.join(rows) + '</div>'


def build_stats_section_html(working, color_map=None, ignore_domains=None):
    """
    Build the statistics section's inner HTML.

    Two independent toggles combine (2 x 3 = 6 panels, one shown at a
    time), plus a sort toggle that applies to whichever panel is visible:

      * granularity -- "Domains" (atomic, split on '+') vs
        "Architectures" (full composite string kept intact);
      * scope -- "All proteins", "Query only" (one count per
        neighborhood, from its reference query gene), or "Neighbors
        only" (every protein except the queries). Keeping query and
        neighbor counts separate matters: query genes are often
        deliberately similar to each other (that's how the neighborhoods
        were assembled in the first place), so folding them into "All"
        can make a domain look far more common among the neighbors than
        it actually is.

    Each panel shows only an interactive bar list (every domain, not just
    a top-N slice -- scrollable, re-sortable). Which domains/architectures
    actually appear in that list is controlled the same way the
    Neighborhoods section controls which blocks show: a "Select" pop-up
    with a checklist (All/None + "Show selected"), one per panel (see
    `build_html_report`/the report template for the matching modal).

    Parameters
    ----------
    working : pandas.DataFrame
        Prepared table (needs 'domain' and 'is_query'); see
        `prepare_dataframe`.
    color_map : dict[str, str] or None
        Domain -> color, to match the figures.
    ignore_domains : list[str] or None
        Forwarded to `compute_domain_stats`.

    Returns
    -------
    (str, str)
        (stats_html, stats_selector_html) -- the panels' HTML and the
        HTML for all 6 selector-checklist groups (one per gran/scope
        combination), each shown/hidden by the report's JS to match
        whichever panel is currently active.
    """
    scopes = [
        ('all', 'All proteins'),
        ('query', 'Query only'),
        ('neighbors', 'Neighbors only'),
    ]
    granularities = [
        ('domain', 'Domains'),
        ('arch', 'Architectures'),
    ]
    notes = {
        'domain': (
            'Each unique single domain &mdash; a composite like <code>GntR+FCD</code> '
            'adds one to <code>GntR</code> and one to <code>FCD</code>.'
        ),
        'arch': (
            'Each full domain architecture counted as one unit &mdash; '
            '<code>GntR+FCD</code> is its own category, separate from '
            '<code>GntR</code> or <code>FCD</code> alone.'
        ),
    }
    scope_notes = {
        'all': 'Counts every protein in the data, query genes and neighbors alike.',
        'query': (
            'Counts only each neighborhood\'s reference query gene &mdash; one per '
            'neighborhood, i.e. how many neighborhoods are built around each domain.'
        ),
        'neighbors': 'Counts every protein EXCEPT the query genes -- what actually surrounds the queries.',
    }

    panels = []
    selector_groups = []
    toggle_gran = ''.join(
        f'<button type="button" class="stats-btn{" active" if i == 0 else ""}" '
        f'data-gran="{key}">{label}</button>'
        for i, (key, label) in enumerate(granularities)
    )
    toggle_scope = ''.join(
        f'<button type="button" class="stats-btn{" active" if i == 0 else ""}" '
        f'data-scope="{key}">{label}</button>'
        for i, (key, label) in enumerate(scopes)
    )

    for scope_key, scope_label in scopes:
        domain_counts, arch_counts = compute_domain_stats(
            working, scope=scope_key, ignore_domains=ignore_domains,
        )
        for gran_key, counts in (('domain', domain_counts), ('arch', arch_counts)):
            bar_list = build_bar_list_html(counts, color_map=color_map)
            active = ' active' if (scope_key == 'all' and gran_key == 'domain') else ''
            panels.append(
                f'<div class="stats-block{active}" data-gran="{gran_key}" data-scope="{scope_key}">'
                f'<p class="stats-note">{notes[gran_key]} {scope_notes[scope_key]} '
                f'{len(counts)} distinct {"domains" if gran_key == "domain" else "architectures"}.</p>'
                f'<div class="stats-chart">{bar_list}</div></div>'
            )

            # matching selector checklist for the "Select domains" pop-up
            sel_display = '' if active else ' style="display:none"'
            items = ''.join(
                f'<div class="nb-sel-item" data-domain="{html.escape(str(row["domain"]).lower())}" '
                f'role="button" tabindex="0">'
                f'<span class="nb-sel-icon">&#9673;</span>'
                f'<span class="nb-sel-name">{html.escape(str(row["domain"]))}</span>'
                f'<span class="nb-sel-sub">{int(row["count"])}</span></div>'
                for _, row in counts.iterrows()
            )
            selector_groups.append(
                f'<div class="stats-sel-group" data-gran="{gran_key}" data-scope="{scope_key}"{sel_display}>'
                f'{items}</div>'
            )

    stats_html = (
        '<div class="stats-controls">'
        '<div class="stats-toggle-wrap">'
        '<span class="stats-toggle-label">View</span>'
        f'<div class="stats-toggle" id="stats-gran-toggle">{toggle_gran}</div>'
        '</div>'
        '<div class="stats-toggle-wrap">'
        '<span class="stats-toggle-label">Scope</span>'
        f'<div class="stats-toggle" id="stats-scope-toggle">{toggle_scope}</div>'
        '</div>'
        '<button type="button" class="nb-icon-btn" id="stats-sel-open">'
        '<svg viewBox="0 0 24 24"><circle cx="12" cy="5" r="1.5"/><circle cx="12" cy="12" r="1.5"/>'
        '<circle cx="12" cy="19" r="1.5"/><line x1="5" y1="5" x2="19" y2="5"/>'
        '<line x1="5" y1="12" x2="19" y2="12"/><line x1="5" y1="19" x2="19" y2="19"/></svg>'
        'Select</button>'
        '<div class="stats-sort-group" id="stats-sort-toggle">'
        '<button type="button" class="stats-btn active" data-sort="count">Sort: count</button>'
        '<button type="button" class="stats-btn" data-sort="alpha">Sort: A&rarr;Z</button>'
        '</div>'
        '<div class="tbl-dl-wrap" id="stats-dl-group" title="Download the distribution currently shown">'
        '<button type="button" class="tbl-dl-btn stats-dl-btn" data-fmt="csv">&#8681; CSV</button>'
        '<button type="button" class="tbl-dl-btn stats-dl-btn" data-fmt="tsv">TSV</button>'
        '<button type="button" class="tbl-dl-btn stats-dl-btn" data-fmt="json">JSON</button>'
        '</div>'
        '</div>'
        f'<div id="stats-panels">{"".join(panels)}</div>'
    )
    return stats_html, ''.join(selector_groups)


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

  .top-tabs{display:flex;flex:1;gap:0;}
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
  .sec-inner{max-width:none;margin:0;padding:28px 40px 72px;}
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
  /* pinned (click-to-keep) card: stays put, takes the mouse so its text
     can be selected, and grows a close button */
  .go-tooltip.pinned{
    pointer-events:auto;user-select:text;padding-right:26px;
    border:1px solid #4a5468;box-shadow:0 10px 30px rgba(0,0,0,.42);
  }
  .go-tooltip .tip-close{
    position:absolute;top:4px;right:7px;cursor:pointer;display:none;
    color:#9aa4b4;font-size:15px;line-height:1;padding:2px 3px;
  }
  .go-tooltip.pinned .tip-close{display:block;}
  .go-tooltip .tip-close:hover{color:#fff;}
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

  /* single global Figure/Table toggle (shared across every selected block) */
  .nb-view-tabs{display:flex;gap:0;border-bottom:2px solid var(--line);margin-bottom:14px;}
  .nb-subtab{
    padding:6px 16px;font-size:13px;cursor:pointer;border:none;background:none;
    color:var(--muted);border-bottom:3px solid transparent;margin-bottom:-2px;
  }
  .nb-subtab.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:600;}
  .nb-subtab:hover{color:var(--ink);}

  /* ONE window holds every currently-selected neighborhood -- no separate
     boxed "windows" per block, just a divider between stacked blocks */
  .nb-window{
    border:1px solid var(--line);border-radius:10px;background:#fff;padding:20px;
  }
  .nb-view{display:none;}
  /* which sub-view shows is controlled centrally via a class on .nb-window,
     so figure vs table applies to the whole merged selection at once */
  .nb-window.view-figure .nb-view-figure{display:block;}
  .nb-window.view-table  .nb-view-table{display:block;}
  .nb-window.view-scale  .nb-view-scale{display:block;}

  /* one block's row inside the single merged figure stack */
  .nb-fig-block{display:none;}
  .nb-fig-block.active{display:block;}
  .nb-fig-block.active ~ .nb-fig-block.active{margin-top:22px;}

  /* same, for the to-scale stack (its own class so block selection can
     toggle both stacks without the slug list being counted twice) */
  .nb-scale-block{display:none;}
  .nb-scale-block.active{display:block;}
  .nb-scale-block.active ~ .nb-scale-block.active{
    margin-top:18px;border-top:1px solid var(--line);padding-top:14px;
  }
  .nb-scale-note{color:var(--muted);font-size:12px;margin:0 0 12px;}

  /* figure wrapper: full width, zoom stretches it (chrome lives on .nb-window now) */
  .nb-fig-scroll{overflow-x:auto;}
  /* figure at natural size; zoom scales it via JS-driven inline width */
  .nb-fig{display:inline-block;width:auto;}
  .nb-fig svg{display:block;width:100%;height:auto;}
  .nb-empty{color:var(--muted);font-size:13px;padding:20px;}
  /* the merged table nested in the shared window doesn't need its own border */
  .nb-window .tbl-card{border:none;border-radius:0;}
  .nb-window .nb-view-table+.nb-view-table,
  .nb-window .tbl-card+.tbl-card{margin-top:12px;}

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
  .stats-controls{display:flex;flex-wrap:wrap;align-items:flex-end;gap:10px 24px;margin-bottom:14px;}
  .stats-toggle-wrap{display:flex;flex-direction:column;gap:4px;}
  .stats-toggle-label{
    font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;
    color:var(--muted);font-weight:600;
  }
  .stats-toggle{display:flex;gap:0;border-bottom:2px solid var(--line);}
  .stats-btn{
    padding:9px 18px;font-size:13px;cursor:pointer;
    border:none;background:none;color:var(--muted);
    border-bottom:3px solid transparent;margin-bottom:-2px;white-space:nowrap;
  }
  .stats-btn.active{color:var(--accent);border-bottom-color:var(--accent);font-weight:600;}
  .stats-btn:hover{color:var(--ink);}
  .stats-sort-group{display:flex;gap:0;border:1px solid var(--line);border-radius:6px;overflow:hidden;}
  .stats-sort-group .stats-btn{padding:6px 12px;border-bottom:none;margin-bottom:0;}
  .stats-sort-group .stats-btn.active{background:var(--accent-soft);}
  .stats-block{display:none;}
  .stats-block.active{display:block;}
  .stats-note{color:var(--muted);font-size:13px;margin:0 0 16px;}
  .stats-chart{margin-bottom:20px;}

  /* interactive HTML/CSS bar list (replaces the old static top-15 SVG) */
  .bar-list{display:flex;flex-direction:column;gap:3px;max-height:480px;overflow-y:auto;padding-right:6px;}
  .bar-row{
    display:grid;grid-template-columns:190px 1fr 48px;align-items:center;gap:10px;
    font-size:12px;font-family:Consolas,"SF Mono",Menlo,monospace;
  }
  .bar-row.hidden-row{display:none;}
  .bar-label{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;text-align:right;color:#222;}
  .bar-track{background:#eee;border-radius:3px;height:16px;overflow:hidden;}
  .bar-fill{display:block;height:100%;border-radius:3px 0 0 3px;}
  .bar-count{color:var(--muted);font-size:11px;}

  /* ---- column filter popup ---- */
  .cfp-backdrop{
    display:none;position:fixed;inset:0;z-index:800;
  }
  .cfp-backdrop.open{display:block;}
  .cfp-popup{
    position:fixed;z-index:801;background:#fff;border:1px solid var(--line);
    border-radius:10px;box-shadow:0 8px 32px rgba(0,0,0,.18);
    min-width:240px;max-width:300px;max-height:80vh;overflow-y:auto;
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
  <p class="sec-desc">Use the <b>&#9776; Select</b> button to choose which neighborhoods are in view -- pick as many as you like, they all show together in one window below. The <b>Figure</b> / <b>To scale</b> / <b>Table</b> toggle switches every visible block at once -- <b>Figure</b> spaces genes evenly so the domain labels read across rows, <b>To scale</b> places them at their real genomic coordinates. Hover any gene arrow for its info window, or click it to keep the window open.</p>

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

  <div class="nb-view-tabs" id="nb-view-tabs">
    <button type="button" class="nb-subtab active" data-view="nb-view-figure">&#9654; Figure</button>
    <button type="button" class="nb-subtab" data-view="nb-view-scale">&#8596; To scale</button>
    <button type="button" class="nb-subtab" data-view="nb-view-table">&#9776; Table</button>
  </div>

  <div class="nb-window view-figure" id="nb-window">
    <div class="nb-view nb-view-figure">
$nb_fig_stack
    </div>
    <div class="nb-view nb-view-scale">
      <p class="nb-scale-note">Genes drawn to <b>biological scale</b>: arrow width is the real gene length and the gaps are the real intergenic distances, along a base-pair ruler. Blocks whose query sits on the minus strand are shown reverse-complemented, so the ruler counts down.</p>
$nb_scale_stack
    </div>
    <div class="nb-view nb-view-table">
$nb_table_card
    </div>
  </div>

</div>
</div>

<!-- ===================== 03 STATISTICS ===================== -->
<div class="page-section" data-page="statistics">
<div class="sec-inner">
  <h1 class="sec-title">Domain statistics</h1>
  <p class="sec-desc"><b>View</b> switches Domains (atomic, split on &lsquo;+&rsquo;) vs Architectures (full composite string); <b>Scope</b> switches All proteins / Query only / Neighbors only. Use <b>&#9776; Select</b> to choose which ones appear in the list, the sort toggle to reorder, and <b>&#8681; CSV/TSV/JSON</b> to download the distribution currently shown.</p>
  <div class="panel">
$stats_html
  </div>
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

<div class="nb-sel-backdrop" id="stats-sel-modal">
  <div class="nb-sel-modal">
    <div class="nb-sel-header">
      <h3 id="stats-sel-title">Select domains</h3>
      <button type="button" class="nb-sel-close" id="stats-sel-close">&times;</button>
    </div>
    <input type="text" class="nb-sel-search" id="stats-sel-search" placeholder="Search domains&hellip;">
    <div class="nb-sel-actions">
      <button type="button" id="stats-sel-all">&#9745; All</button>
      <button type="button" id="stats-sel-none">&#9744; None</button>
      <button type="button" class="nb-sel-show" id="stats-sel-apply">Show selected</button>
    </div>
    <div class="nb-sel-list" id="stats-sel-list">
$stats_selector
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

<footer>Made by <b>S(H)ARP</b> &mdash; Biosynthetic Gene Cluster Analysis</footer>

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

  // ── shared tooltip ─────────────────────────────────────
  // One card, two modes: it follows the mouse while hovering, and a
  // click *pins* it -- frozen where it was opened and ignoring every
  // hover -- so its text can be read and selected without keeping the
  // pointer on the gene. A pinned card closes on its ×, on a click
  // anywhere outside it, or on Escape.
  var tipPinned = false;
  function moveTip(e) {
    var x=e.clientX+14, y=e.clientY+14, r=tip.getBoundingClientRect();
    if (x+r.width>window.innerWidth)  x=e.clientX-r.width-14;
    if (y+r.height>window.innerHeight) y=e.clientY-r.height-14;
    tip.style.left=x+'px'; tip.style.top=y+'px';
  }
  function showTip(content, e) {
    tip.innerHTML = content + '<span class="tip-close" title="Close">&times;</span>';
    tip.style.display = 'block';
    moveTip(e);
  }
  function hideTip()  { tip.style.display='none'; }
  function unpinTip() { tipPinned=false; tip.classList.remove('pinned'); hideTip(); }
  function pinTip(el, e) {
    tipPinned = false;             // let showTip/moveTip place the fresh card
    showTip(el.dataset.tip, e);
    tipPinned = true;
    tip.classList.add('pinned');
  }
  function attachTip(el) {
    on(el, 'mouseenter', function(e){ if (!tipPinned) showTip(el.dataset.tip, e); });
    on(el, 'mousemove',  function(e){ if (!tipPinned) moveTip(e); });
    on(el, 'mouseleave', function(){ if (!tipPinned) hideTip(); });
  }
  on(tip, 'click', function (e) {
    if (e.target.classList.contains('tip-close')) { unpinTip(); return; }
    e.stopPropagation();           // clicks inside a pinned card keep it open
  });
  on(document, 'click', function () { if (tipPinned) unpinTip(); });
  var goMarkers = qsa('.go-marker');
  goMarkers.forEach(function(m){ attachTip(m); on(m,'click',function(){ unpinTip(); openNeighborhood(m.dataset.block); }); });
  // ── per-protein info window (same card as the genome overview) ───────
  // Hover → the card follows the mouse; click → it stays put until closed.
  qsa('.nb-gene').forEach(function (el) {
    attachTip(el);
    on(el, 'click', function (e) {
      e.stopPropagation();         // don't let the document handler unpin it
      pinTip(el, e);
    });
  });

  // ── neighborhood selection (single merged figure + single merged table) ──
  var figBlocks   = qsa('.nb-fig-block');
  var scaleBlocks = qsa('.nb-scale-block');
  var selItems    = qsa('.nb-sel-item');
  var allSlugs  = figBlocks.map(function(f){ return f.dataset.block; });

  // Which slugs are currently shown (set by applySelection)
  var activeSet = {};

  // Single global Figure/Table toggle -- applies to the whole merged
  // window at once via a class on .nb-window.
  var nbWindowEl = qs('#nb-window');
  function setNbView(name) {
    if (nbWindowEl) nbWindowEl.className = 'nb-window ' + name.replace('nb-view-', 'view-');
    qsa('#nb-view-tabs .nb-subtab').forEach(function (t) { t.classList.toggle('active', t.dataset.view === name); });
  }
  qsa('#nb-view-tabs .nb-subtab').forEach(function (t) {
    on(t, 'click', function () { setNbView(t.dataset.view); });
  });

  function updateCrumb() {
    var n = Object.keys(activeSet).length;
    var el = qs('#nb-crumb-count');
    if (el) {
      el.innerHTML = n > 0 ?
        (n + '&nbsp;neighborhood' + (n>1?'s':'') + '&nbsp;shown') : '';
    }
  }

  function applySelection(slugs) {
    activeSet = {};
    slugs.forEach(function(s){ activeSet[s]=true; });
    figBlocks.forEach(function(f){ f.classList.toggle('active', !!activeSet[f.dataset.block]); });
    scaleBlocks.forEach(function(f){ f.classList.toggle('active', !!activeSet[f.dataset.block]); });
    selItems.forEach(function(i){
      i.classList.toggle('active', !!activeSet[i.dataset.block]);
      var cb = i.querySelector('input[type=checkbox]');
      if (cb) cb.checked = !!activeSet[i.dataset.block];
    });
    goMarkers.forEach(function(m){ m.classList.toggle('selected', !!activeSet[m.dataset.block]); });
    updateCrumb();
    // re-filter the merged table's rows to match the new selection
    var nbTable = qs('#nb-table-card');
    if (nbTable && nbTable._applyFilter) nbTable._applyFilter();
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

  // prev / next (cycle through *all* blocks, showing only that one)
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
  function nbNaturalWidth(fig) {
    if (fig._natW) return fig._natW;
    var svg = fig.querySelector('svg');
    var w = 0;
    if (svg) {
      var vb = svg.getAttribute('viewBox');
      if (vb) { var p = vb.split(/[ ,]+/); w = parseFloat(p[2]) || 0; }
      if (!w) w = svg.getBoundingClientRect().width || 0;
    }
    fig._natW = w || 600;
    return fig._natW;
  }
  function nbLayout() {
    qsa('.nb-fig').forEach(function(f){
      var w = nbNaturalWidth(f) * nbZoom;
      f.style.width = Math.max(60, Math.round(w)) + 'px';
    });
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
    if(e.key==='Escape') { closeSel(); closeStatsSel(); unpinTip(); }
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

    // position popup near the th; final clamp happens after content is
    // built below, once the popup's real (content-dependent) size is known
    var rect = thEl.getBoundingClientRect();
    cfpPopup.style.left = Math.max(8, Math.min(rect.left, window.innerWidth - 310)) + 'px';
    cfpPopup.style.top  = (rect.bottom + 4) + 'px';
    cfpPopup.style.display = 'block';
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

    // final clamp: now that the popup's real content (num controls or the
    // value checklist) is in place, measure it and make sure it's fully
    // inside the viewport -- flip above the header if there's no room
    // below, and never let it hang off the left/right/bottom edges.
    requestAnimationFrame(function () {
      var pr = cfpPopup.getBoundingClientRect();
      var left = Math.max(8, Math.min(rect.left, window.innerWidth - pr.width - 8));
      var top = rect.bottom + 4;
      if (top + pr.height > window.innerHeight - 8) {
        var above = rect.top - pr.height - 4;
        top = above >= 8 ? above : Math.max(8, window.innerHeight - pr.height - 8);
      }
      cfpPopup.style.left = left + 'px';
      cfpPopup.style.top = top + 'px';
    });
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
        // rows in the merged neighborhoods table carry data-block; only
        // show rows whose block is currently selected (other tables, e.g.
        // statistics, don't set data-block so this check is a no-op there)
        var failB = r.dataset.block && typeof activeSet !== 'undefined' && !activeSet[r.dataset.block];
        r.classList.toggle('hidden-row', failG || failC || failB);
      });
      updateCount();
    }
    on(globalFilter, 'input', applyFilter);
    card._applyFilter = applyFilter;

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

  // ── statistics: granularity x scope toggle + domain selection + sort ───
  var statsGran = 'domain', statsScope = 'all', statsSort = 'count';
  // per-panel (gran|scope) selected-domain sets; a panel with no entry here
  // yet means "everything selected" (the default, until the user customizes it)
  var statsSelections = {};
  function statsKey() { return statsGran + '|' + statsScope; }

  function applyStatsSelection() {
    var active = qs('.stats-block.active');
    if (!active) return;
    var sel = statsSelections[statsKey()];
    qsa('.bar-row', active).forEach(function (r) {
      var shown = !sel || !!sel[r.dataset.domain];
      r.classList.toggle('hidden-row', !shown);
    });
  }
  function applyStatsSort() {
    var active = qs('.stats-block.active');
    if (!active) return;
    var list = active.querySelector('.bar-list');
    if (!list) return;
    var rows = qsa('.bar-row', list);
    rows.sort(function (a, b) {
      if (statsSort === 'alpha') return a.dataset.name.localeCompare(b.dataset.name);
      return parseFloat(b.dataset.count) - parseFloat(a.dataset.count);
    });
    rows.forEach(function (r) { list.appendChild(r); });
  }
  function showStatsPanel() {
    qsa('.stats-block').forEach(function (b) {
      b.classList.toggle('active', b.dataset.gran === statsGran && b.dataset.scope === statsScope);
    });
    applyStatsSelection();
    applyStatsSort();
  }
  qsa('#stats-gran-toggle .stats-btn').forEach(function (btn) {
    on(btn, 'click', function () {
      statsGran = btn.dataset.gran;
      qsa('#stats-gran-toggle .stats-btn').forEach(function (b) { b.classList.toggle('active', b === btn); });
      showStatsPanel();
    });
  });
  qsa('#stats-scope-toggle .stats-btn').forEach(function (btn) {
    on(btn, 'click', function () {
      statsScope = btn.dataset.scope;
      qsa('#stats-scope-toggle .stats-btn').forEach(function (b) { b.classList.toggle('active', b === btn); });
      showStatsPanel();
    });
  });
  qsa('#stats-sort-toggle .stats-btn').forEach(function (btn) {
    on(btn, 'click', function () {
      statsSort = btn.dataset.sort;
      qsa('#stats-sort-toggle .stats-btn').forEach(function (b) { b.classList.toggle('active', b === btn); });
      applyStatsSort();
    });
  });
  showStatsPanel();

  // ── "Select domains" pop-up (mirrors the neighborhoods Select pop-up) ──
  var statsSelModal = qs('#stats-sel-modal');
  var statsSelGroups = qsa('.stats-sel-group');

  // inject checkboxes into every group's items, once
  statsSelGroups.forEach(function (g) {
    qsa('.nb-sel-item', g).forEach(function (i) {
      var cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = true;
      on(cb, 'click', function (e) { e.stopPropagation(); });
      on(cb, 'change', function (e) { e.stopPropagation(); });
      i.appendChild(cb);
    });
  });

  function currentStatsSelGroup() {
    return statsSelGroups.filter(function (g) {
      return g.dataset.gran === statsGran && g.dataset.scope === statsScope;
    })[0];
  }
  function openStatsSel() {
    var group = currentStatsSelGroup();
    statsSelGroups.forEach(function (g) { g.style.display = (g === group) ? '' : 'none'; });
    if (group) {
      var sel = statsSelections[statsKey()];
      qsa('.nb-sel-item', group).forEach(function (i) {
        var cb = i.querySelector('input[type=checkbox]');
        if (cb) cb.checked = !sel || !!sel[i.dataset.domain];
      });
    }
    var titleEl = qs('#stats-sel-title');
    if (titleEl) titleEl.textContent = 'Select ' + (statsGran === 'domain' ? 'domains' : 'architectures');
    if (statsSelModal) statsSelModal.classList.add('open');
    var s = qs('#stats-sel-search');
    if (s) { s.value = ''; s.dispatchEvent(new Event('input')); s.focus(); }
  }
  function closeStatsSel() { if (statsSelModal) statsSelModal.classList.remove('open'); }
  on(qs('#stats-sel-open'), 'click', openStatsSel);
  on(qs('#stats-sel-close'), 'click', closeStatsSel);
  on(statsSelModal, 'click', function (e) { if (e.target === statsSelModal) closeStatsSel(); });
  on(qs('#stats-sel-search'), 'input', function () {
    var term = this.value.toLowerCase();
    var group = currentStatsSelGroup();
    if (!group) return;
    qsa('.nb-sel-item', group).forEach(function (i) {
      i.classList.toggle('hidden-row', term && i.textContent.toLowerCase().indexOf(term) === -1);
    });
  });
  on(qs('#stats-sel-all'), 'click', function () {
    var group = currentStatsSelGroup();
    if (!group) return;
    qsa('.nb-sel-item', group).forEach(function (i) {
      if (!i.classList.contains('hidden-row')) {
        var cb = i.querySelector('input[type=checkbox]'); if (cb) cb.checked = true;
      }
    });
  });
  on(qs('#stats-sel-none'), 'click', function () {
    var group = currentStatsSelGroup();
    if (!group) return;
    qsa('.nb-sel-item', group).forEach(function (i) {
      var cb = i.querySelector('input[type=checkbox]'); if (cb) cb.checked = false;
    });
  });
  on(qs('#stats-sel-apply'), 'click', function () {
    var group = currentStatsSelGroup();
    if (!group) { closeStatsSel(); return; }
    var sel = {};
    qsa('.nb-sel-item', group).forEach(function (i) {
      var cb = i.querySelector('input[type=checkbox]');
      if (cb && cb.checked) sel[i.dataset.domain] = true;
    });
    statsSelections[statsKey()] = sel;
    applyStatsSelection();
    closeStatsSel();
  });

  // ── download the distribution currently shown (respects the domain
  // selection and sort order, same "download what's visible" pattern as
  // the neighborhoods table) ──────────────────────────────────────────
  function statsVisibleRows() {
    var active = qs('.stats-block.active');
    if (!active) return [];
    return qsa('.bar-row', active)
      .filter(function (r) { return !r.classList.contains('hidden-row'); })
      .map(function (r) { return { domain: r.dataset.name, count: r.dataset.count }; });
  }
  function statsDownload(fmt) {
    var rows = statsVisibleRows();
    var base = 'distribution_' + statsGran + '_' + statsScope;
    var d;
    if (fmt === 'tsv') {
      var t = ['domain\tcount'];
      rows.forEach(function (r) { t.push(r.domain + '\t' + r.count); });
      d = { txt: t.join('\n'), mime: 'text/tab-separated-values', ext: 'tsv' };
    } else if (fmt === 'json') {
      d = { txt: JSON.stringify(rows, null, 2), mime: 'application/json', ext: 'json' };
    } else {
      var c = ['"domain","count"'];
      rows.forEach(function (r) { c.push('"' + r.domain.replace(/"/g, '""') + '","' + r.count + '"'); });
      d = { txt: c.join('\n'), mime: 'text/csv', ext: 'csv' };
    }
    var a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([d.txt], { type: d.mime }));
    a.download = base + '.' + d.ext;
    a.click();
  }
  qsa('.stats-dl-btn').forEach(function (btn) {
    on(btn, 'click', function () { statsDownload(btn.dataset.fmt); });
  });


  // ── init ──────────────────────────────────────────────────────────────
  initBases();
  goLayout();
  nbLayout();
  // Respect whichever blocks the report was built with as active
  // (default_view='all' bakes every block active, 'first' bakes just one).
  var initialSlugs = figBlocks
    .filter(function (f) { return f.classList.contains('active'); })
    .map(function (f) { return f.dataset.block; });
  if (initialSlugs.length === 0 && allSlugs.length > 0) initialSlugs = [allSlugs[0]];
  applySelection(initialSlugs);
  on(window,'resize', function(){ initBases(); goLayout(); });
})();
</script>
</body>
</html>
""")


def render_neighborhood_table_card(df, group_col='block_id', filename='neighborhoods.csv',
                                    max_rows=None):
    """
    Build ONE sortable/filterable/downloadable table-card containing every
    row of `df`, with each `<tr>` tagged `data-block="<slug>"` (see
    `_slug`). This is the single merged table shown in the Neighborhoods
    section's Table view: the report's JS shows/hides rows to match
    whichever neighborhoods are currently selected, on top of the table's
    own text/column filters, rather than swapping between separate
    per-block tables.

    Parameters
    ----------
    df : pandas.DataFrame
        Full input table (all blocks).
    group_col : str
        Column identifying each block; used to compute each row's slug.
    filename : str
        Suggested name for downloads.
    max_rows : int or None
        Optional row cap. `None` (the default) embeds every row.

    Returns
    -------
    str
        HTML for one `.tbl-card`, id'd `nb-table-card` so the report's JS
        can target it directly.
    """
    shown = df if max_rows is None else df.head(max_rows)
    slugs = shown[group_col].map(_slug) if group_col in shown.columns else [''] * len(shown)

    header_cells = ''.join(f'<th>{html.escape(str(c))}</th>' for c in shown.columns)
    body_rows = []
    for slug, row in zip(slugs, shown.itertuples(index=False, name=None)):
        cells = ''.join('<td>' + ('' if pd.isna(v) else html.escape(str(v))) + '</td>' for v in row)
        body_rows.append(f'<tr data-block="{slug}">{cells}</tr>')

    table_html = (
        '<table class="data-tbl">'
        f'<thead><tr>{header_cells}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody>'
        '</table>'
    )
    note = ''
    if max_rows is not None and len(df) > max_rows:
        note = f'<p class="table-note">Showing the first {max_rows:,} of {len(df):,} rows.</p>'

    return (
        f'<div class="tbl-card" id="nb-table-card" data-filename="{html.escape(filename)}">'
        '<div class="tbl-controls">'
        '<input type="text" class="tbl-filter" placeholder="Filter rows...">'
        '<span class="tbl-count"></span>'
        '<div class="tbl-dl-wrap">'
        '<button type="button" class="tbl-dl-btn" data-fmt="csv">&#8681; CSV</button>'
        '<button type="button" class="tbl-dl-btn" data-fmt="tsv">TSV</button>'
        '<button type="button" class="tbl-dl-btn" data-fmt="json">JSON</button>'
        '</div></div>'
        f'<div class="tbl-scroll">{table_html}{note}</div>'
        '</div>'
    )


def build_neighborhood_panels(extents, block_svgs, table_card, default_view='all',
                               scale_svgs=None):
    """
    Build the neighborhoods section's inner fragments for the pop-up
    (icon) multi-select model.

    Every currently-selected neighborhood renders together inside ONE
    shared window: a single merged figure stack (each selected block's
    SVG, one below another, all inside one bordered/scrollable box) and
    ONE merged, sortable/filterable/downloadable table (every row of the
    full input table, tagged by block; the report's JS shows only the
    rows belonging to selected blocks). A single Figure/Table toggle
    switches the whole window's sub-view at once. There's no separate
    "All neighborhoods" entry -- checking every box in the pop-up (via
    its "All" button) achieves the same result through the same
    mechanism.

    Parameters
    ----------
    extents : pandas.DataFrame
        Output of `compute_block_extents`; drives order and labels.
    block_svgs : dict[str, str]
        Slug -> annotated SVG markup.
    table_card : str
        The single merged table-card HTML, from
        `render_neighborhood_table_card`.
    default_view : str
        Which blocks are active on load: 'all' (every block) or 'first'
        (just the first block).
    scale_svgs : dict[str, str] or None
        Slug -> to-scale SVG markup (see `render_scaled_svgs_by_block`),
        for the window's "To scale" sub-view. None leaves that stack
        empty.

    Returns
    -------
    (str, str, str)
        (fig_stack_html, scale_stack_html, selector_items_html).
    """
    if extents.empty:
        empty = '<div class="nb-empty">No blocks to show.</div>'
        return empty, empty, ''

    scale_svgs = scale_svgs or {}

    nucleotides = list(dict.fromkeys(extents['nucleotide']))
    first_slug = _slug(extents.iloc[0]['ID'])

    fig_blocks = []
    scale_blocks = []
    selector_items = []

    for nucleotide in nucleotides:
        rows = extents[extents['nucleotide'] == nucleotide]
        selector_items.append(
            f'<div class="nb-sel-group-title">{html.escape(str(nucleotide))} '
            f'<span class="nb-sel-count">{len(rows)}</span></div>'
        )
        for _, block in rows.iterrows():
            slug = _slug(block['ID'])
            label = block['query_pid'] if block['query_pid'] is not None else block['ID']
            active = ' active' if (default_view == 'all' or slug == first_slug) else ''

            fig_blocks.append(
                f'<div class="nb-fig-block{active}" data-block="{slug}">'
                f'<div class="nb-fig">{block_svgs.get(slug, "")}</div></div>'
            )
            scale_blocks.append(
                f'<div class="nb-scale-block{active}" data-block="{slug}">'
                f'<div class="nb-fig">{scale_svgs.get(slug, "")}</div></div>'
            )
            selector_items.append(
                f'<div class="nb-sel-item" data-block="{slug}" role="button" tabindex="0" '
                f'title="{html.escape(str(block["ID"]))}">'
                f'<span class="nb-sel-icon">&#9673;</span>'
                f'<span class="nb-sel-name">{html.escape(str(label))}</span>'
                f'<span class="nb-sel-sub">{int(block["n_genes"])} genes</span></div>'
            )

    fig_stack = (
        f'<div class="nb-fig-scroll"><div id="nb-fig-stack">{"".join(fig_blocks)}</div></div>'
    )
    scale_stack = (
        f'<div class="nb-fig-scroll"><div id="nb-scale-stack">{"".join(scale_blocks)}</div></div>'
    )
    return fig_stack, scale_stack, ''.join(selector_items)


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


def resolve_logo(value):
    """
    Turn a `header_logo`/`footer_logo` argument into inline SVG markup.

    Accepts either ready-to-embed SVG markup (returned unchanged) or the
    path of an SVG file (`~` expanded), which is read through
    `read_svg_logo`. A path that does not exist gives an empty string,
    so a report still builds on a machine without the branding file.

    Parameters
    ----------
    value : str or None
        SVG markup, a path to an SVG file, or None/'' for "no logo".

    Returns
    -------
    str
        SVG markup, or '' when there is nothing to embed.
    """
    if not value:
        return ''
    if value.lstrip().startswith('<'):
        return value
    path = os.path.expanduser(value)
    if not os.path.exists(path):
        return ''
    return read_svg_logo(path)


def build_html_report(df, output_file='operon_report.html', title='Gene Neighborhood Report',
                       group_col='block_id', org_col='organism', label_col='pfam',
                       rename_map=None, custom_colors=None, max_colors=5, ignore_domains=None,
                       nucleotide_col='nucleotide', start_col='start', end_col='end',
                       length_col='nlen', operon_kwargs=None, max_table_rows=None,
                       work_dir=None, default_view='all',
                       software_name='S(H)ARP',
                       header_logo=SHARP_HEADER_LOGO_PATH, footer_logo=None):
    """
    Build one self-contained, interactive HTML page for a single input
    table: a zoomable genome-wide overview with hover tooltips, a
    multi-select neighborhoods section, and a statistics section.

    Neighborhoods section behavior:

      * a "Select" pop-up lets the user check any number of neighborhoods
        (with "All"/"None" bulk buttons); every checked one renders
        together in a single window;
      * that window holds ONE merged figure (each selected block's SVG
        stacked in one bordered/scrollable box, not separate boxes per
        block) and ONE merged, sortable/filterable/downloadable table
        (every row of the full input table, restricted to the rows of
        whichever blocks are selected); a single Figure/Table toggle
        switches the whole window's sub-view at once;
      * each figure spans the full width of the panel and does NOT shrink
        with the number of neighbors in it; the zoom controls stretch it;
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
        Raw input table.
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
        Row cap for the merged neighborhoods table. `None` (the default)
        embeds every row -- unrestricted.
    work_dir : str or None
        Where intermediate SVGs are written. If None, a temporary
        directory is used and cleaned up afterwards.
    default_view : str, default 'all'
        Which neighborhoods are selected on load: 'all' (every block) or
        'first' (just the first one).
    software_name : str, default 'S(H)ARP'
        Name shown in the report footer ("Made by ...").
    header_logo : str or None
        Logo for the top-nav brand slot: either the path of an SVG file
        (`~` expanded) or ready-to-embed SVG markup -- see
        `resolve_logo`. Defaults to `SHARP_HEADER_LOGO_PATH`; a path
        that does not exist, or `None`, leaves the slot empty.
    footer_logo : str or None
        Logo for the footer. Defaults to None (no footer logo). Same
        conventions as `header_logo`.

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

    own_tmp = work_dir is None
    tmp_dir = work_dir or tempfile.mkdtemp(prefix='operon_report_')
    try:
        # The left-hand id/block/organism label (see `show_row_label` on
        # `neighborhood_figure`) stays on by default here too -- it's the
        # per-figure identifier. Pass `operon_kwargs=dict(show_row_label=False)`
        # to turn it off if it's ever not wanted.
        per_block_operon_kwargs = dict(operon_kwargs)
        per_block_operon_kwargs.update(
            org_col=org_col, label_col=label_col, rename_map=rename_map,
        )
        block_svgs = render_neighborhood_svgs_by_block(
            df, group_col=group_col, color_map=color_map,
            operon_kwargs=per_block_operon_kwargs, tmp_dir=tmp_dir,
        )
    finally:
        if own_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # ONE merged table for every block, tagged per row by block -- the
    # report's JS shows only the rows for whichever blocks are selected.
    nb_table_card = render_neighborhood_table_card(
        df, group_col=group_col, filename='neighborhoods.csv', max_rows=max_table_rows,
    )

    # Same blocks, drawn to real genomic scale for the "To scale" sub-view.
    scale_svgs = render_scaled_svgs_by_block(
        working, color_map=color_map, nucleotide_col=nucleotide_col,
        start_col=start_col, end_col=end_col,
        normalize_orientation=operon_kwargs.get('normalize_orientation', True),
    )

    nb_fig_stack, nb_scale_stack, nb_selector = build_neighborhood_panels(
        extents, block_svgs, nb_table_card, default_view=default_view,
        scale_svgs=scale_svgs,
    )

    # Statistics (granularity x scope, all computed inside build_stats_section_html)
    stats_html, stats_selector = build_stats_section_html(
        working, color_map=color_map, ignore_domains=ignore_domains,
    )

    n_blocks = df[group_col].nunique() if group_col in df.columns else 'NA'
    html_doc = HTML_REPORT_TEMPLATE.substitute(
        title=html.escape(title),
        n_genes=f'{len(df):,}',
        n_blocks=n_blocks,
        genome_overview=genome_overview,
        nb_fig_stack=nb_fig_stack,
        nb_scale_stack=nb_scale_stack,
        nb_table_card=nb_table_card,
        nb_selector=nb_selector,
        stats_html=stats_html,
        stats_selector=stats_selector,
        software_name=html.escape(software_name),
        header_logo_html=resolve_logo(header_logo),
    )

    with open(output_file, 'w') as f:
        f.write(html_doc)

    return output_file