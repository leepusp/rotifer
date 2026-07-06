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
    build_tabs_and_panels_html    tabs + per-block panels (figure + table)
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
    --bg: #faf9f6;
    --panel: #ffffff;
    --ink: #1f2430;
    --muted: #6b7280;
    --line: #e4e1d8;
    --accent: #2a6f77;
    --accent-soft: #e4f0f1;
    --selected: #fdecdc;
    --selected-line: #e07b39;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    background: var(--bg);
    color: var(--ink);
    font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.45;
  }
  header {
    padding: 32px 40px 22px;
    border-bottom: 1px solid var(--line);
    background: var(--panel);
  }
  header h1 { margin: 0 0 6px; font-size: 25px; letter-spacing: -0.01em; }
  header .meta { color: var(--muted); font-size: 13px; font-family: Consolas, "SF Mono", Menlo, monospace; }
  header .meta b { color: var(--accent); }
  main { max-width: 1280px; margin: 0 auto; padding: 30px 40px 80px; }
  section { margin-bottom: 52px; }
  section .eyebrow {
    font-family: Consolas, "SF Mono", Menlo, monospace;
    font-size: 12px; letter-spacing: 0.08em; color: var(--accent); margin: 0 0 4px;
  }
  section h2 { margin: 0 0 4px; font-size: 19px; }
  section p.desc { color: var(--muted); margin: 0 0 16px; font-size: 14px; max-width: 760px; }
  .panel {
    background: var(--panel); border: 1px solid var(--line);
    border-radius: 8px; padding: 18px;
  }

  /* ---- genome overview ---- */
  .go-controls { display: flex; align-items: center; gap: 8px; margin-bottom: 14px; flex-wrap: wrap; }
  .go-controls button {
    border: 1px solid var(--line); background: #fff; color: var(--ink);
    border-radius: 6px; padding: 4px 11px; font-size: 14px; cursor: pointer; line-height: 1;
  }
  .go-controls button:hover { background: var(--accent-soft); border-color: var(--accent); }
  .go-controls button:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
  #go-zoom-val { font-family: Consolas, "SF Mono", Menlo, monospace; font-size: 13px; min-width: 38px; text-align: center; }
  .go-hint { color: var(--muted); font-size: 12px; margin-left: 6px; }
  .go-track { display: grid; grid-template-columns: 150px 1fr; align-items: center; margin: 16px 0; gap: 12px; }
  .go-track-label { text-align: right; font-family: Consolas, "SF Mono", Menlo, monospace; overflow: hidden; }
  .go-contig { font-size: 13px; color: #222; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
  .go-len { font-size: 11px; color: var(--muted); }
  .go-viewport { overflow-x: auto; overflow-y: hidden; cursor: grab; padding: 18px 0; }
  .go-viewport.grabbing { cursor: grabbing; }
  .go-inner { position: relative; height: var(--track-h, 14px); }
  .go-axis {
    position: absolute; left: 0; right: 0; top: 50%; transform: translateY(-50%);
    height: 6px; background: #e3e3e3; border: 0.5px solid #bbb; border-radius: 3px;
  }
  .go-marker {
    position: absolute; top: 0; height: 100%;
    border: 1.2px solid; border-radius: 2px; cursor: pointer; box-sizing: border-box;
    transition: transform 0.06s ease;
  }
  .go-marker:hover { transform: scaleY(1.45); z-index: 3; }
  .go-marker.selected { box-shadow: 0 0 0 2px var(--selected-line); z-index: 2; }
  .go-tooltip {
    position: fixed; display: none; z-index: 1000; pointer-events: none;
    background: #1f2430; color: #fff; border-radius: 7px; padding: 9px 11px;
    font-size: 12px; max-width: 320px; box-shadow: 0 6px 22px rgba(0,0,0,0.25);
    font-family: Consolas, "SF Mono", Menlo, monospace; line-height: 1.5;
  }
  .go-tooltip b { display: block; margin-bottom: 3px; font-size: 12.5px; }
  .go-tooltip .t-row { display: block; color: #cdd3da; }
  .go-tooltip .t-query { color: #ffd28a; font-weight: bold; }

  /* hoverable proteins inside neighborhood figures */
  .nb-gene { cursor: pointer; }
  .nb-gene:hover { filter: brightness(0.92); }
  .nb-gene:hover polygon, .nb-gene:hover ellipse { stroke-width: 2.4px; }

  /* ---- neighborhoods: toolbar + tabs + full-width figure ---- */
  .nb-toolbar { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; margin-bottom: 12px; }
  .nb-btn {
    border: 1px solid var(--line); background: #fff; color: var(--ink);
    border-radius: 6px; padding: 6px 12px; font-size: 13px; cursor: pointer; line-height: 1;
  }
  .nb-btn:hover { background: var(--accent-soft); border-color: var(--accent); }
  .nb-btn:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
  .nb-zoom { display: flex; align-items: center; gap: 6px; margin-left: auto; }
  .nb-zoom button {
    border: 1px solid var(--line); background: #fff; border-radius: 6px;
    padding: 4px 10px; font-size: 14px; cursor: pointer; line-height: 1;
  }
  .nb-zoom button:hover { background: var(--accent-soft); border-color: var(--accent); }
  #nb-zoom-val { font-family: Consolas, "SF Mono", Menlo, monospace; font-size: 13px; min-width: 46px; text-align: center; }

  .nb-tabs {
    display: flex; gap: 4px; overflow-x: auto; border-bottom: 2px solid var(--line);
    margin-bottom: 14px; padding-bottom: 0;
  }
  .nb-tab {
    border: 1px solid var(--line); border-bottom: none; background: #f4f2ec; color: #555;
    border-radius: 7px 7px 0 0; padding: 7px 14px; font-size: 12.5px; cursor: pointer;
    white-space: nowrap; font-family: Consolas, "SF Mono", Menlo, monospace;
    margin-bottom: -2px; transition: background 0.1s ease;
  }
  .nb-tab:hover { background: var(--accent-soft); }
  .nb-tab.active { background: #fff; color: var(--accent); border-color: var(--line); border-bottom: 2px solid #fff; font-weight: bold; }
  .nb-tab.nb-tab-all { color: var(--accent); }
  .nb-tab.nb-hidden { display: none; }

  .nb-panel { display: none; }
  .nb-panel.active { display: block; }
  .nb-panel-head { font-family: Consolas, "SF Mono", Menlo, monospace; font-size: 13px; color: var(--muted); margin-bottom: 10px; }
  .nb-panel-head b { color: var(--ink); }
  /* the figure occupies the full window width and does NOT shrink with
     gene count; the zoom controls stretch .nb-fig beyond 100% */
  .nb-fig-scroll { overflow-x: auto; border: 1px solid var(--line); border-radius: 8px; padding: 14px; background: var(--panel); }
  .nb-fig { width: 100%; }
  .nb-fig svg { display: block; width: 100%; height: auto; }
  .nb-subtable-wrap { margin-top: 16px; }
  .nb-subtable-title { font-family: Consolas, "SF Mono", Menlo, monospace; font-size: 12px; color: var(--muted); margin-bottom: 6px; }
  .nb-subtable { max-height: 340px; overflow: auto; border: 1px solid var(--line); border-radius: 8px; }
  .nb-empty { color: var(--muted); font-size: 13px; padding: 20px; }

  /* ---- filter pop-up ---- */
  .nb-modal-backdrop {
    display: none; position: fixed; inset: 0; background: rgba(20,24,32,0.42);
    z-index: 900; align-items: center; justify-content: center;
  }
  .nb-modal-backdrop.open { display: flex; }
  .nb-modal {
    background: #fff; border-radius: 12px; width: min(460px, 92vw); max-height: 78vh;
    display: flex; flex-direction: column; box-shadow: 0 18px 50px rgba(0,0,0,0.3); overflow: hidden;
  }
  .nb-modal header { display: flex; align-items: center; padding: 16px 18px; border-bottom: 1px solid var(--line); }
  .nb-modal header h3 { margin: 0; font-size: 15px; }
  .nb-modal header button { margin-left: auto; border: none; background: none; font-size: 22px; line-height: 1; cursor: pointer; color: var(--muted); }
  .nb-modal-search { margin: 12px 18px 0; padding: 8px 11px; border: 1px solid var(--line); border-radius: 6px; font-size: 13px; font-family: Consolas, "SF Mono", Menlo, monospace; }
  .nb-modal-actions { display: flex; gap: 8px; padding: 10px 18px; }
  .nb-modal-actions button { font-size: 12px; padding: 4px 10px; border: 1px solid var(--line); background: #fff; border-radius: 6px; cursor: pointer; }
  .nb-modal-actions button:hover { background: var(--accent-soft); }
  .nb-modal-list { overflow-y: auto; padding: 4px 18px 18px; }
  .nb-filter-group-title { font-family: Consolas, "SF Mono", Menlo, monospace; font-size: 12px; color: var(--accent); margin: 12px 0 4px; }
  .nb-filter-count { color: var(--muted); }
  .nb-filter-item { display: flex; align-items: center; gap: 8px; padding: 5px 4px; font-size: 12.5px; font-family: Consolas, "SF Mono", Menlo, monospace; cursor: pointer; border-radius: 5px; }
  .nb-filter-item:hover { background: var(--accent-soft); }
  .nb-filter-item.hidden-row { display: none; }
  .nb-filter-item .nb-filter-sub { color: var(--muted); font-size: 11px; margin-left: auto; }

  /* ---- data table ---- */
  .search-row { display: flex; align-items: center; gap: 10px; margin-bottom: 12px; }
  .search-row input {
    flex: 0 1 320px; padding: 8px 12px; border: 1px solid var(--line);
    border-radius: 6px; font-size: 13px; font-family: Consolas, "SF Mono", Menlo, monospace;
  }
  .search-row input:focus { outline: 2px solid var(--accent); outline-offset: 1px; }
  .search-row .count { color: var(--muted); font-size: 12px; }
  .table-scroll { max-height: 560px; overflow: auto; border: 1px solid var(--line); border-radius: 8px; }
  table { border-collapse: collapse; width: 100%; font-size: 12.5px; font-family: Consolas, "SF Mono", Menlo, monospace; }
  thead th {
    position: sticky; top: 0; background: var(--accent-soft); color: var(--ink);
    text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--line); white-space: nowrap;
  }
  tbody td { padding: 6px 10px; border-bottom: 1px solid var(--line); white-space: nowrap; }
  tbody tr:nth-child(even) { background: #fbfbf9; }
  tbody tr.hidden-row { display: none; }
  .table-note { color: var(--muted); font-size: 12px; margin-top: 10px; }
  footer { text-align: center; color: var(--muted); font-size: 12px; padding: 24px 0 48px; }

  @media (max-width: 820px) {
    .go-track { grid-template-columns: 110px 1fr; }
    .nb-zoom { margin-left: 0; width: 100%; }
  }
</style>
</head>
<body>
<header>
  <h1>$title</h1>
  <div class="meta"><b>$n_genes</b> genes &middot; <b>$n_blocks</b> neighborhoods</div>
</header>
<main>

  <section>
    <p class="eyebrow">01 &middot; GENOME-WIDE VIEW</p>
    <h2>Where each neighborhood sits in the genome</h2>
    <p class="desc">One track per contig, scaled to its own length. Each marker is one neighborhood. Hover a marker for its details, zoom in to separate crowded hits, and click a marker to open that neighborhood below.</p>
    <div class="panel">
$genome_overview
    </div>
  </section>

  <section>
    <p class="eyebrow">02 &middot; NEIGHBORHOODS</p>
    <h2>Gene neighborhood detail</h2>
    <p class="desc">Each figure spans the full width of the panel. Use <b>Filter neighborhoods</b> to choose which results appear as tabs, switch tabs to inspect one at a time (its figure and the input rows that built it), and zoom the figure as needed. Hover any protein for its info window; the query protein is flagged as such. Query genes are outlined in red, neighbors always point right, and opposite-strand neighbors collapse to small triangles.</p>
    <div class="nb-toolbar">
      <button type="button" class="nb-btn" id="nb-filter-open">&#9776;&nbsp;Filter neighborhoods</button>
      <div class="nb-zoom">
        <button type="button" id="nb-zoom-out" title="Zoom out">&minus;</button>
        <span id="nb-zoom-val">100%</span>
        <button type="button" id="nb-zoom-in" title="Zoom in">+</button>
        <button type="button" id="nb-zoom-reset" title="Reset zoom">reset</button>
        <button type="button" id="nb-zoom-fit" title="Fit to width">fit width</button>
      </div>
    </div>
    <div class="nb-tabs">
$nb_tabs
    </div>
    <div class="nb-content">
$nb_panels
    </div>
  </section>

  <section>
    <p class="eyebrow">03 &middot; DATA</p>
    <h2>Input table</h2>
    <p class="desc">The raw table this report was built from.</p>
    <div class="search-row">
      <input type="text" id="table-search" placeholder="Filter rows...">
      <span class="count" id="table-count"></span>
    </div>
    <div class="table-scroll">
$table_html
    </div>
  </section>

</main>

<div class="nb-modal-backdrop" id="nb-modal">
  <div class="nb-modal">
    <header><h3>Show neighborhoods</h3><button type="button" id="nb-modal-close" title="Close">&times;</button></header>
    <input type="text" class="nb-modal-search" id="nb-modal-search" placeholder="Search...">
    <div class="nb-modal-actions">
      <button type="button" id="nb-check-all">Select all</button>
      <button type="button" id="nb-check-none">Clear</button>
    </div>
    <div class="nb-modal-list">
$nb_filter_items
    </div>
  </div>
</div>

<footer>Made by <b>$software_name</b> &middot; gene-neighborhood analysis</footer>
<script>
(function () {
  // ---------- genome overview: zoom / pan / tooltip ----------
  var zoom = 1;
  var tracks = Array.prototype.slice.call(document.querySelectorAll('.go-track'));
  var tip = document.getElementById('go-tooltip');
  var zoomVal = document.getElementById('go-zoom-val');

  function layout() {
    tracks.forEach(function (tr) {
      var vp = tr.querySelector('.go-viewport');
      var inner = tr.querySelector('.go-inner');
      var base = tr._base || vp.clientWidth || 600;
      var w = base * zoom;
      inner.style.width = w + 'px';
      var markers = inner.querySelectorAll('.go-marker');
      for (var i = 0; i < markers.length; i++) {
        var m = markers[i];
        var fs = parseFloat(m.dataset.fs), fe = parseFloat(m.dataset.fe);
        var left = fs * w;
        var width = (fe - fs) * w;
        if (width < 5) width = 5;
        m.style.left = left + 'px';
        m.style.width = width + 'px';
      }
    });
  }

  function setZoom(z) {
    zoom = Math.min(500, Math.max(1, z));
    if (zoomVal) zoomVal.textContent = (zoom < 10 ? zoom.toFixed(1) : Math.round(zoom)) + 'x';
    layout();
  }

  function initBases() {
    tracks.forEach(function (tr) {
      var vp = tr.querySelector('.go-viewport');
      tr._base = vp.clientWidth || 600;
    });
  }

  var zin = document.getElementById('go-zoom-in');
  var zout = document.getElementById('go-zoom-out');
  var zreset = document.getElementById('go-zoom-reset');
  if (zin) zin.addEventListener('click', function () { setZoom(zoom * 1.6); });
  if (zout) zout.addEventListener('click', function () { setZoom(zoom / 1.6); });
  if (zreset) zreset.addEventListener('click', function () {
    setZoom(1);
    tracks.forEach(function (tr) { tr.querySelector('.go-viewport').scrollLeft = 0; });
  });

  tracks.forEach(function (tr) {
    var vp = tr.querySelector('.go-viewport');

    vp.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        var rect = vp.getBoundingClientRect();
        var base = tr._base || vp.clientWidth;
        var frac = (vp.scrollLeft + (e.clientX - rect.left)) / (base * zoom);
        setZoom(zoom * (e.deltaY < 0 ? 1.2 : 1 / 1.2));
        vp.scrollLeft = frac * (base * zoom) - (e.clientX - rect.left);
      } else {
        e.preventDefault();
        vp.scrollLeft += (e.deltaY + e.deltaX);
      }
    }, { passive: false });

    var dragging = false, startX = 0, startScroll = 0;
    vp.addEventListener('mousedown', function (e) {
      if (e.target.classList.contains('go-marker')) return;
      dragging = true; startX = e.clientX; startScroll = vp.scrollLeft;
      vp.classList.add('grabbing');
    });
    window.addEventListener('mousemove', function (e) {
      if (dragging) vp.scrollLeft = startScroll - (e.clientX - startX);
    });
    window.addEventListener('mouseup', function () {
      dragging = false; vp.classList.remove('grabbing');
    });
  });

  function moveTip(e) {
    var x = e.clientX + 14, y = e.clientY + 14;
    var r = tip.getBoundingClientRect();
    if (x + r.width > window.innerWidth) x = e.clientX - r.width - 14;
    if (y + r.height > window.innerHeight) y = e.clientY - r.height - 14;
    tip.style.left = x + 'px';
    tip.style.top = y + 'px';
  }
  function attachTip(el) {
    el.addEventListener('mouseenter', function (e) {
      tip.innerHTML = el.dataset.tip;
      tip.style.display = 'block';
      moveTip(e);
    });
    el.addEventListener('mousemove', moveTip);
    el.addEventListener('mouseleave', function () { tip.style.display = 'none'; });
  }

  var markers = Array.prototype.slice.call(document.querySelectorAll('.go-marker'));
  markers.forEach(function (m) {
    attachTip(m);
    m.addEventListener('click', function () { selectBlock(m.dataset.block); });
  });

  // per-protein info windows inside the neighborhood figures
  Array.prototype.slice.call(document.querySelectorAll('.nb-gene')).forEach(attachTip);

  // ---------- neighborhoods: tabs, zoom, filter pop-up ----------
  var panels = Array.prototype.slice.call(document.querySelectorAll('.nb-panel'));
  var tabs = Array.prototype.slice.call(document.querySelectorAll('.nb-tab'));

  function selectBlock(slug) {
    var tab = tabs.filter(function (t) { return t.dataset.block === slug; })[0];
    if (tab && tab.classList.contains('nb-hidden')) {
      // a marker click can target a filtered-out neighborhood: reveal it
      tab.classList.remove('nb-hidden');
      var cb = document.querySelector('.nb-filter-item input[data-block="' + slug + '"]');
      if (cb) cb.checked = true;
    }
    panels.forEach(function (p) { p.classList.toggle('active', p.dataset.block === slug); });
    tabs.forEach(function (t) { t.classList.toggle('active', t.dataset.block === slug); });
    markers.forEach(function (m) { m.classList.toggle('selected', m.dataset.block === slug); });
  }
  window._selectBlock = selectBlock;

  tabs.forEach(function (t) {
    t.addEventListener('click', function () { selectBlock(t.dataset.block); });
  });

  // neighborhood figure zoom: figures always fill the panel width (they
  // don't shrink with gene count); zoom stretches the figure wrapper and
  // the panel scrolls horizontally.
  var nbZoom = 1;
  var nbZoomVal = document.getElementById('nb-zoom-val');
  function nbLayout() {
    document.querySelectorAll('.nb-fig').forEach(function (f) { f.style.width = (nbZoom * 100) + '%'; });
    if (nbZoomVal) nbZoomVal.textContent = Math.round(nbZoom * 100) + '%';
  }
  function nbSetZoom(z) { nbZoom = Math.min(8, Math.max(1, z)); nbLayout(); }
  var b;
  if ((b = document.getElementById('nb-zoom-in'))) b.addEventListener('click', function () { nbSetZoom(nbZoom * 1.3); });
  if ((b = document.getElementById('nb-zoom-out'))) b.addEventListener('click', function () { nbSetZoom(nbZoom / 1.3); });
  if ((b = document.getElementById('nb-zoom-reset'))) b.addEventListener('click', function () { nbSetZoom(1); });
  if ((b = document.getElementById('nb-zoom-fit'))) b.addEventListener('click', function () { nbSetZoom(1); });

  // filter pop-up
  var modal = document.getElementById('nb-modal');
  var filterInputs = Array.prototype.slice.call(document.querySelectorAll('.nb-filter-item input'));
  function openModal() { if (modal) modal.classList.add('open'); }
  function closeModal() { if (modal) modal.classList.remove('open'); }
  if ((b = document.getElementById('nb-filter-open'))) b.addEventListener('click', openModal);
  if ((b = document.getElementById('nb-modal-close'))) b.addEventListener('click', closeModal);
  if (modal) modal.addEventListener('click', function (e) { if (e.target === modal) closeModal(); });

  function applyFilter() {
    filterInputs.forEach(function (cb) {
      var slug = cb.dataset.block;
      var tab = tabs.filter(function (t) { return t.dataset.block === slug; })[0];
      if (tab) tab.classList.toggle('nb-hidden', !cb.checked);
    });
    // if the active tab was hidden, jump to the first visible one
    var active = tabs.filter(function (t) { return t.classList.contains('active'); })[0];
    if (!active || active.classList.contains('nb-hidden')) {
      var firstVisible = tabs.filter(function (t) { return !t.classList.contains('nb-hidden'); })[0];
      if (firstVisible) selectBlock(firstVisible.dataset.block);
    }
  }
  filterInputs.forEach(function (cb) { cb.addEventListener('change', applyFilter); });
  if ((b = document.getElementById('nb-check-all'))) b.addEventListener('click', function () {
    filterInputs.forEach(function (cb) { cb.checked = true; }); applyFilter();
  });
  if ((b = document.getElementById('nb-check-none'))) b.addEventListener('click', function () {
    filterInputs.forEach(function (cb) { cb.checked = false; }); applyFilter();
  });
  var modalSearch = document.getElementById('nb-modal-search');
  if (modalSearch) modalSearch.addEventListener('input', function () {
    var term = modalSearch.value.toLowerCase();
    Array.prototype.slice.call(document.querySelectorAll('.nb-filter-item')).forEach(function (it) {
      it.classList.toggle('hidden-row', it.textContent.toLowerCase().indexOf(term) === -1);
    });
  });
  window.addEventListener('keydown', function (e) { if (e.key === 'Escape') closeModal(); });

  // ---------- data table filter ----------
  var input = document.getElementById('table-search');
  var table = document.getElementById('data-table');
  if (input && table) {
    var rows = Array.prototype.slice.call(table.querySelectorAll('tbody tr'));
    var countEl = document.getElementById('table-count');
    var update = function () {
      var visible = rows.filter(function (r) { return !r.classList.contains('hidden-row'); }).length;
      if (countEl) countEl.textContent = visible + ' of ' + rows.length + ' rows';
    };
    input.addEventListener('input', function () {
      var term = input.value.toLowerCase();
      rows.forEach(function (r) {
        r.classList.toggle('hidden-row', r.textContent.toLowerCase().indexOf(term) === -1);
      });
      update();
    });
    update();
  }

  // ---------- init ----------
  initBases();
  setZoom(1);
  nbLayout();
  window.addEventListener('resize', function () { initBases(); layout(); });
})();
</script>
</body>
</html>
""")


def build_tabs_and_panels_html(extents, block_svgs, block_tables=None,
                                combined_svg=None, default_view='all'):
    """
    Build the fragments the (sidebar-free) neighborhoods section needs:
    a single horizontal row of tabs, the stack of panels (only the
    active one shown), and the checklist that populates the "filter"
    pop-up. Selecting a tab shows that neighborhood's figure *and* the
    slice of the input table used to draw it.

    Parameters
    ----------
    extents : pandas.DataFrame
        Output of `compute_block_extents`; drives order and labels.
    block_svgs : dict[str, str]
        Slug -> annotated SVG markup, from
        `render_neighborhood_svgs_by_block`.
    block_tables : dict[str, str] or None
        Slug -> HTML `<table>` of that block's input rows. Shown under
        the figure in each per-block panel.
    combined_svg : str or None
        If given, an extra "All neighborhoods" tab (data-block
        `__all__`) is added first, showing every block at once.
    default_view : str
        Which tab is active on load: 'all' (the combined tab, if present)
        or 'first' (the first block).

    Returns
    -------
    (str, str, str)
        (tabs_html, panels_html, filter_items_html).
    """
    if extents.empty:
        empty = '<div class="nb-empty">No blocks to show.</div>'
        return '', empty, ''

    block_tables = block_tables or {}
    nucleotides = list(dict.fromkeys(extents['nucleotide']))
    has_combined = combined_svg is not None
    selected = '__all__' if (default_view == 'all' and has_combined) else _slug(extents.iloc[0]['ID'])

    tabs = []
    panels = []
    filter_items = []

    def fig_block(svg):
        # The figure fills the full width of the panel (it does NOT scale
        # with the number of neighbors); the inner wrapper's width is
        # what the neighborhood zoom controls stretch.
        return f'<div class="nb-fig-scroll"><div class="nb-fig">{svg}</div></div>'

    if has_combined:
        act = ' active' if selected == '__all__' else ''
        tabs.append(
            f'<button type="button" class="nb-tab nb-tab-all{act}" data-block="__all__">'
            f'&#9776; All neighborhoods</button>'
        )
        panels.append(
            f'<div class="nb-panel{act}" data-block="__all__">'
            f'<div class="nb-panel-head"><b>All neighborhoods</b> &middot; {len(extents)} blocks, every block stacked</div>'
            f'{fig_block(combined_svg)}</div>'
        )

    for nucleotide in nucleotides:
        rows = extents[extents['nucleotide'] == nucleotide]
        filter_group = [
            f'<div class="nb-filter-group-title">{html.escape(str(nucleotide))} '
            f'<span class="nb-filter-count">{len(rows)}</span></div>'
        ]
        for _, block in rows.iterrows():
            slug = _slug(block['ID'])
            label = block['query_pid'] if block['query_pid'] is not None else block['ID']
            domain = block['query_domain'] if block['query_domain'] is not None else '-'
            act = ' active' if slug == selected else ''

            tabs.append(
                f'<button type="button" class="nb-tab{act}" data-block="{slug}" '
                f'title="{html.escape(str(block["ID"]))}">{html.escape(str(label))}</button>'
            )

            head = (
                f'<div class="nb-panel-head"><b>{html.escape(str(label))}</b> &middot; '
                f'{html.escape(str(block["ID"]))} &middot; {html.escape(str(domain))} &middot; '
                f'{html.escape(str(block["org_name"]))} &middot; '
                f'{block["block_start"]:,.0f}&ndash;{block["block_end"]:,.0f} bp &middot; '
                f'{int(block["n_genes"])} genes</div>'
            )
            table_html = block_tables.get(slug, '')
            table_block = (
                f'<div class="nb-subtable-wrap"><div class="nb-subtable-title">Input rows for this neighborhood</div>'
                f'<div class="nb-subtable">{table_html}</div></div>'
            ) if table_html else ''
            panels.append(
                f'<div class="nb-panel{act}" data-block="{slug}">{head}'
                f'{fig_block(block_svgs.get(slug, ""))}{table_block}</div>'
            )

            checked = ' checked' if slug != '__never__' else ''
            filter_group.append(
                f'<label class="nb-filter-item"><input type="checkbox" data-block="{slug}"{checked}> '
                f'<span>{html.escape(str(label))}</span>'
                f'<span class="nb-filter-sub">{int(block["n_genes"])} genes</span></label>'
            )
        filter_items.append(f'<div class="nb-filter-group">{"".join(filter_group)}</div>')

    return ''.join(tabs), ''.join(panels), ''.join(filter_items)


def build_html_report(df, output_file='operon_report.html', title='Gene Neighborhood Report',
                       group_col='block_id', org_col='organism', label_col='pfam',
                       rename_map=None, custom_colors=None, max_colors=5, ignore_domains=None,
                       nucleotide_col='nucleotide', start_col='start', end_col='end',
                       length_col='nlen', operon_kwargs=None, max_table_rows=2000,
                       work_dir=None, include_combined=True, default_view='all',
                       software_name='SHARP'):
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
    software_name : str, default 'SHARP'
        Name shown in the report footer ("Made by ...").

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

    # One input-table slice per block, shown under its figure.
    block_tables = {
        _slug(block_id): render_dataframe_html(block_df, table_id=None, css_class='nb-subtable-table')
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

    nb_tabs, nb_panels, nb_filter_items = build_tabs_and_panels_html(
        extents, block_svgs, block_tables=block_tables,
        combined_svg=combined_svg, default_view=default_view,
    )

    n_blocks = df[group_col].nunique() if group_col in df.columns else 'NA'
    html_doc = HTML_REPORT_TEMPLATE.substitute(
        title=html.escape(title),
        n_genes=f'{len(df):,}',
        n_blocks=n_blocks,
        genome_overview=genome_overview,
        nb_tabs=nb_tabs,
        nb_panels=nb_panels,
        nb_filter_items=nb_filter_items,
        table_html=render_dataframe_html(df, max_rows=max_table_rows),
        software_name=html.escape(software_name),
    )

    with open(output_file, 'w') as f:
        f.write(html_doc)

    return output_file