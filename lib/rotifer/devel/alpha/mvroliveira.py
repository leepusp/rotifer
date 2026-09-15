#!/usr/bin/env python3
"""
mvroliveira.py - Sequence similarity networks (SSNs) from all-vs-all BLAST.

    all_vs_all_search(seqobj) -> raw hit table (BLAST, DIAMOND or MMseqs2 --
                                  your choice, `tool=`)
    build_ssn(hits)       -> networkx.Graph filtered by any criteria
    sequence_ssn(seqobj)  -> does both steps at once
    export_ssn(G, path)   -> save graph to disk (any tool that reads
                                  GraphML/TSV)
    to_cytoscape(G)           -> push straight into a running Cytoscape
                                  Desktop (via py4cytoscape/cyREST)
    get_selected_ids()        -> IDs currently selected in Cytoscape,
                                  to filter your sequence object

    build_tree(seqobj)        -> tree, default engine='fasttree' (fast,
                                  approximate ML). engine='iqtree2' also
                                  available, with UFBoot + SH-aLRT (1000/1000)
                                  already set as soon as you switch to it.
                                  Returns plain Newick text -- no extra
                                  parsing dependency needed.

    build_ssn_and_tree(seqobj) -> runs sequence_ssn + build_tree in one
                                  call; {'graph', 'hits', 'tree'}

    closeness_scan(hits)      -> scans cutoffs, average closeness centrality
                                  per cutoff (Hornung & Terrapon 2023) -- an
                                  objective way to pick a subfamily cutoff
    plot_closeness_scan(scan) -> PNG of the curve above, peaks marked

Split BLAST from ssn-building in two steps so BLAST (the expensive
part) only runs once, while you try different edge criteria freely on
the same hit table.

Visualization deliberately does NOT reimplement layout/styling/selection
in Python -- Cytoscape already does all of that, well, and your team is
already using it. `to_cytoscape` just hands the graph over; use
Cytoscape's own Layout menu, Style panel and Select panel from there.

Importing this module also attaches `.to_ssn()`, `.to_tree()` and
`.build_ssn_and_tree()` to every rotifer `sequence` object
(monkey-patch, see bottom of this file):

    seqs = rdbs.sequence("teste.fasta")
    G = seqs.to_ssn('rede.graphml')
    tree = seqs.align().to_tree(cpu=8)
    result = seqs.align().build_ssn_and_tree(cpu=8)   # both at once: {'graph','hits','tree'}

NOTE: these only exist AFTER this module has been imported at least
once -- the patch happens at import time, not before.
"""

import os
import json
import tempfile
from subprocess import Popen, PIPE, STDOUT

import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm

import rotifer
import rotifer.devel.beta.sequence as rdbs
logger = rotifer.logging.getLogger(__name__)

HIT_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
]

# Columns where LOWER = better hit. Everything else is higher = better.
LOWER_IS_BETTER = {'evalue', 'mismatch', 'gapopen'}

TOOLS = ('blast', 'diamond', 'mmseqs')


def _search_commands(tool, program, moltype, fasta, db, out, evalue, cpu, max_target_seqs, extra_args):
    """
    Build the (db_command, search_command) shell strings for one of the
    three supported engines. All three are asked to emit the exact same
    columns, in the same order (HIT_COLUMNS), so the rest of the
    pipeline never needs to know which engine actually ran.

    Returns db_command as None when the engine doesn't need a separate
    database-building step (mmseqs easy-search builds its own).
    """
    cols = ' '.join(HIT_COLUMNS)

    if tool == 'blast':
        dbtype = 'prot' if program in ('blastp', 'blastx') else 'nucl'
        db_cmd = f'makeblastdb -in {fasta} -dbtype {dbtype} -out {db}'
        cmd = (f'{program} -query {fasta} -db {db} -evalue {evalue} '
               f'-num_threads {cpu} -outfmt "6 {cols}" -out {out}')
        if max_target_seqs:
            cmd += f' -max_target_seqs {max_target_seqs}'

    elif tool == 'diamond':
        if moltype == 'nucl':
            raise ValueError(
                "DIAMOND doesn't do nucleotide-vs-nucleotide search. "
                "Use tool='blast' or tool='mmseqs' for that, or moltype='prot' "
                "with program='blastx' for a translated search."
            )
        db_cmd = f'diamond makedb --in {fasta} -d {db}'
        cmd = f'diamond {program} -q {fasta} -d {db} -e {evalue} -p {cpu} -f 6 {cols} -o {out}'
        if max_target_seqs:
            cmd += f' --max-target-seqs {max_target_seqs}'

    elif tool == 'mmseqs':
        search_type = {'prot': 1, 'nucl': 3}[moltype]
        db_cmd = None  # easy-search builds its own temp databases
        cmd = (f'mmseqs easy-search {fasta} {fasta} {out} {db}_tmp '
               f'--search-type {search_type} -e {evalue} --threads {cpu} '
               f'--format-output "{",".join(HIT_COLUMNS)}"')
        if max_target_seqs:
            cmd += f' --max-seqs {max_target_seqs}'

    else:
        raise ValueError(f"tool must be one of {TOOLS}, got {tool!r}")

    if extra_args:
        cmd += f' {extra_args}'
    return db_cmd, cmd


def all_vs_all_search(seqobj, tool='blast', program=None, moltype='prot', cpu=8,
                       evalue=10, max_target_seqs=None, keep_self_hits=False,
                       extra_args='', progress=True):
    """
    Search every sequence in seqobj against all the others -- with
    BLAST, DIAMOND, or MMseqs2, your choice.

    All three engines are asked for the same output columns in the
    same order, so the returned table looks identical regardless of
    which one ran; nothing downstream (build_ssn, etc.) needs to
    know or care.

    Parameters
    ----------
    seqobj : rotifer.devel.beta.sequence.sequence
        Alignments are handled fine -- gaps are stripped before searching.
    tool : str, default 'blast'
        'blast', 'diamond', or 'mmseqs'.
    program : str or None, default None
        Only meaningful for tool='blast' or tool='diamond' (both have
        'blastp'/'blastx' subcommands/programs). Defaults to 'blastp'
        when left as None. Ignored for tool='mmseqs' (use `moltype`
        instead).
    moltype : str, default 'prot'
        'prot' or 'nucl'. Used by mmseqs (--search-type) and to sanity
        -check diamond (which can't do pure nucleotide-vs-nucleotide
        search). For tool='blast', `program` still decides this
        directly (blastp/blastn/blastx) -- moltype is only consulted
        there to pick makeblastdb's -dbtype.
    cpu : int, default 8
        Threads (mapped to each engine's own flag).
    evalue : float, default 10
        Each engine's own reporting threshold. Keep this loose; filter
        for real in build_ssn's `criteria` instead (no need to
        search again).
    max_target_seqs : int or None, default None
        Caps hits per query (mapped to each engine's own flag).
    keep_self_hits : bool, default False
        Drop each sequence's hit against itself.
    extra_args : str, default ''
        Any extra flags for the underlying engine, appended verbatim.
    progress : bool, default True
        Show a tqdm progress bar with the current pipeline stage.

    Returns
    -------
    pandas.DataFrame
        One row per pairwise hit: qseqid, sseqid, pident, length,
        mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore.

    Examples
    --------
    >>> hits = rdao.all_vs_all_search(seqs, tool='diamond', cpu=16)
    >>> hits = rdao.all_vs_all_search(seqs, tool='mmseqs', cpu=16)
    """
    if len(seqobj) < 2:
        logger.warning('all_vs_all_search: need at least 2 sequences')
        return pd.DataFrame(columns=HIT_COLUMNS)
    if tool not in TOOLS:
        raise ValueError(f"tool must be one of {TOOLS}, got {tool!r}")
    if program is None:
        program = 'blastp'

    steps = 4 if tool in ('blast', 'diamond') else 3  # mmseqs skips the db-build step
    pbar = tqdm(total=steps, desc='all_vs_all_search', unit='step', disable=not progress)

    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)

        pbar.set_description('writing fasta')
        fasta = f'{tmp}/seqs.fasta'
        seqobj.to_file(fasta, remove_gaps=True)
        pbar.update(1)

        out = f'{tmp}/hits.tsv'
        db_cmd, search_cmd = _search_commands(
            tool, program, moltype, fasta, f'{tmp}/db', out,
            evalue, cpu, max_target_seqs, extra_args,
        )

        if db_cmd:
            pbar.set_description(f'{tool} makedb')
            Popen(db_cmd, stdout=PIPE, stderr=STDOUT, shell=True).communicate()
            pbar.update(1)

        pbar.set_description(f'running {tool}')
        Popen(search_cmd, stdout=PIPE, stderr=STDOUT, shell=True).communicate()
        pbar.update(1)

        pbar.set_description('parsing hits')
        if not os.path.exists(out) or not os.path.getsize(out):
            logger.warning(f'all_vs_all_search: {tool} produced no hits')
            df = pd.DataFrame(columns=HIT_COLUMNS)
        else:
            df = pd.read_csv(out, sep='\t', names=HIT_COLUMNS)
        pbar.update(1)
    os.chdir(cwd)
    pbar.close()

    if not keep_self_hits and not df.empty:
        df = df.query('qseqid != sseqid').reset_index(drop=True)
    return df


def build_ssn(hits, criteria='bitscore > 70', keep_best='bitscore',
                   directed=False, add_community=False, resolution_parameter=0.05,
                   auto_cutoff=False, cutoff_column='bitscore', cutoff_steps=50,
                   min_component_size=2, max_singleton_fraction=0.5,
                   auto_start=True, peak_selection='first', progress=True):
    """
    Turn a BLAST hit table into a networkx graph.

    Parameters
    ----------
    hits : pandas.DataFrame
        Output of all_vs_all_search.
    criteria : str, default 'bitscore > 70'
        Pandas .query() expression deciding which pairs become edges.
        Any column works, e.g. 'evalue < 1e-3' or 'bitscore > 70 and pident > 30'.
        Ignored entirely when `auto_cutoff=True`.
    keep_best : str or None, default 'bitscore'
        BLAST isn't symmetric (A-vs-B != B-vs-A). This column breaks
        ties between the two directions of a pair. Direction (higher
        or lower wins) is picked automatically.
    add_community : bool, default False
        Also run Louvain/Leiden clustering (rotifer.devel.alpha.net_functions)
        and store it as node attributes 'Louvain' / 'Leiden' -- handy to
        color by in Cytoscape afterwards (see to_cytoscape's color_by).
    auto_cutoff : bool, default False
        If True, `criteria` is ignored -- the SSN is built AFTER
        running closeness_scan on `cutoff_column` (Hornung & Terrapon
        2023, doi:10.1371/journal.pcbi.1010881) and picking a cutoff
        from its local peaks, instead of a value chosen by hand. The
        chosen cutoff and the full scan are stashed on the returned
        graph so you can check what was picked and why:

            G.graph['auto_cutoff']       # e.g. 'bitscore > 84.3'
            G.graph['closeness_scan']    # the full scan DataFrame

        Raises ValueError if the scan finds no usable peak (can mean
        the family genuinely has no clean subfamilies at this scale --
        see closeness_scan/plot_closeness_scan to inspect the curve
        yourself and fall back to an explicit `criteria` instead).
    cutoff_column : str, default 'bitscore'
        Column scanned when `auto_cutoff=True`.
    cutoff_steps : int, default 50
        Passed to closeness_scan's `n_steps` when `auto_cutoff=True`.
    min_component_size, max_singleton_fraction
        Passed straight to closeness_scan -- see there. These are what
        keep `auto_cutoff` from landing on an over-fragmented "many
        tiny blobs" cutoff (a known limitation of raw closeness
        centrality, discussed explicitly in the paper).
    peak_selection : str, default 'first'
        Which local peak to use when `auto_cutoff=True` and several
        exist: 'first' (broadest/most permissive clean split -- the
        conservative default), 'last' (strictest split found), or
        'max' (the single highest closeness value among all peaks).
    progress : bool, default True
        Show a tqdm progress bar with the current pipeline stage.

    Returns
    -------
    networkx.Graph (or DiGraph). Every BLAST column is kept as an edge
    attribute, regardless of what `criteria` filtered on.

    Examples
    --------
    >>> G = rdao.build_ssn(hits, criteria='bitscore > 70')       # manual, as before
    >>> G = rdao.build_ssn(hits, auto_cutoff=True)               # cutoff picked for you
    >>> G.graph['auto_cutoff'], G.graph['closeness_scan']            # inspect the choice
    """
    if auto_cutoff:
        scan = closeness_scan(hits, column=cutoff_column, n_steps=cutoff_steps,
                               min_component_size=min_component_size,
                               max_singleton_fraction=max_singleton_fraction,
                               auto_start=auto_start, progress=progress)
        # Callers treat this ValueError as "no real peak, fall back to a
        # single community", so every failure mode routes through it.
        if 'avg_closeness' not in scan.columns:
            raise ValueError(
                "auto_cutoff: closeness_scan returned no 'avg_closeness' column -- "
                "the scan itself failed before producing usable data. Inspect "
                "closeness_scan(hits) directly to see why."
            )
        valid = scan.dropna(subset=['avg_closeness']).reset_index(drop=True)
        if len(valid) < 3:
            raise ValueError(
                'auto_cutoff: not enough valid points in the closeness scan to find a '
                'peak -- check hits/cutoff_column, or inspect closeness_scan(hits) yourself.'
            )
        v = valid['avg_closeness'].to_numpy()
        is_peak = np.zeros(len(v), dtype=bool)
        for i in range(1, len(v) - 1):
            if v[i] >= v[i - 1] and v[i] >= v[i + 1] and (v[i] > v[i - 1] or v[i] > v[i + 1]):
                is_peak[i] = True
        if v[0] > v[1]:
            is_peak[0] = True
        if v[-1] > v[-2]:
            is_peak[-1] = True
        # A high closeness value with n_components == 1 just means "still
        # one connected blob" -- not a split into subfamilies. A single
        # dense blob scores well on closeness too (see Hornung & Terrapon
        # 2023, Fig 1: the fully-connected end of their artificial data
        # scores the same as the fully-separated end). Without this
        # filter, the loosest cutoff tested is very often mistaken for
        # "the peak", handing back the whole family as one undivided blob.
        #
        # over_fragmented excludes the OTHER extreme: past a point, most
        # sequences have splintered into singletons, and whatever remains
        # connected scores well on closeness for an equally uninformative
        # reason (a tiny leftover clique, not a real subfamily). The scan
        # itself still reports these points -- for the curve to be worth
        # looking at, someone has to be able to see the whole shape and
        # judge for themselves -- this filter only keeps them out of the
        # AUTOMATIC choice.
        not_over_frag = ~valid['over_fragmented'].to_numpy() if 'over_fragmented' in valid.columns \
                        else np.ones(len(v), dtype=bool)
        peaks = valid.loc[is_peak & (valid['n_components'] >= 2) & not_over_frag]
        if peaks.empty:
            raise ValueError(
                'auto_cutoff: no local peak found -- this family may not have clean '
                'subfamilies at this scale. Inspect plot_closeness_scan(closeness_scan(hits)) '
                'and fall back to an explicit `criteria` if so.'
            )
        if peak_selection == 'first':
            chosen = peaks.iloc[0]
        elif peak_selection == 'last':
            chosen = peaks.iloc[-1]
        elif peak_selection == 'max':
            chosen = peaks.loc[peaks['avg_closeness'].idxmax()]
        else:
            raise ValueError("peak_selection must be 'first', 'last', or 'max'")
        criteria = f'{cutoff_column} > {chosen["cutoff"]}'
        logger.info(f'auto_cutoff: picked "{criteria}" (avg_closeness={chosen["avg_closeness"]:.3f})')

    steps = 3 + int(add_community)
    pbar = tqdm(total=steps, desc='build_ssn', unit='step', disable=not progress)

    pbar.set_description('filtering by criteria')
    df = hits.query(criteria).copy() if len(hits) else hits.copy()
    if df.empty:
        logger.warning('build_ssn: no pairs passed the criteria')
    pbar.update(1)

    pbar.set_description('resolving symmetric pairs')
    if not directed and keep_best and not df.empty:
        ascending = keep_best in LOWER_IS_BETTER
        df['source'], df['target'] = np.where(
            df['qseqid'] > df['sseqid'],
            [df['qseqid'], df['sseqid']],
            [df['sseqid'], df['qseqid']],
        )
        df = df.sort_values(keep_best, ascending=ascending).drop_duplicates(['source', 'target'])
    else:
        df = df.rename({'qseqid': 'source', 'sseqid': 'target'}, axis=1)
    pbar.update(1)

    pbar.set_description('building graph')
    edge_attrs = [c for c in df.columns if c not in ('source', 'target')]
    G = nx.from_pandas_edgelist(
        df, source='source', target='target',
        edge_attr=edge_attrs or None,
        create_using=(nx.DiGraph if directed else nx.Graph)(),
    )
    pbar.update(1)

    if auto_cutoff:
        G.graph['auto_cutoff'] = criteria
        G.graph['closeness_scan'] = scan

    if add_community:
        pbar.set_description('detecting communities')
        if G.number_of_edges() > 0:
            from rotifer.devel.alpha import net_functions
            weight = 'bitscore' if 'bitscore' in edge_attrs else (edge_attrs[0] if edge_attrs else None)
            c = net_functions.get_network_community(G, weight=weight, resolution_parameter=resolution_parameter)
            nx.set_node_attributes(G, c['Louvain'][0], name='Louvain')
            nx.set_node_attributes(G, c['Leidein'][0], name='Leiden')
        pbar.update(1)

    pbar.close()
    return G


def sequence_ssn(seqobj, tool='blast', program=None, moltype='prot', cpu=8,
                      search_evalue=10, max_target_seqs=None, keep_self_hits=False,
                      criteria='bitscore > 70', keep_best='bitscore',
                      directed=False, add_community=False, resolution_parameter=0.05,
                      auto_cutoff=False, cutoff_column='bitscore', cutoff_steps=50,
                      min_component_size=2, max_singleton_fraction=0.5,
                      auto_start=True, peak_selection='first', progress=True):
    """
    Convenience wrapper: all_vs_all_search + build_ssn in one call.
    Returns (G, hits) -- keep `hits` to rebuild the graph with a
    different `criteria` without searching again.

    `tool` picks the engine ('blast', 'diamond', or 'mmseqs') -- see
    all_vs_all_search for what each parameter maps to per engine.
    `auto_cutoff`/`cutoff_column`/`cutoff_steps`/`min_component_size`/
    `max_singleton_fraction`/`auto_start`/`peak_selection` are forwarded
    to build_ssn -- see there for details.
    """
    hits = all_vs_all_search(seqobj, tool=tool, program=program, moltype=moltype,
                              cpu=cpu, evalue=search_evalue, max_target_seqs=max_target_seqs,
                              keep_self_hits=keep_self_hits, progress=progress)
    G = build_ssn(hits, criteria=criteria, keep_best=keep_best, directed=directed,
                       add_community=add_community, resolution_parameter=resolution_parameter,
                       auto_cutoff=auto_cutoff, cutoff_column=cutoff_column,
                       cutoff_steps=cutoff_steps, min_component_size=min_component_size,
                       max_singleton_fraction=max_singleton_fraction,
                       auto_start=auto_start,
                       peak_selection=peak_selection, progress=progress)
    return G, hits


def export_ssn(G, path, fmt='graphml'):
    """
    Save an SSN to disk -- for Cytoscape, Gephi, or any other tool.

    Parameters
    ----------
    G : networkx.Graph
    path : str
    fmt : str, default 'graphml'
        'graphml' -- single file, keeps every node/edge attribute.
                     In Cytoscape: File > Import > Network from File.
        'tsv'     -- plain edge table (source, target, + attributes).

    Notes
    -----
    GraphML only supports primitive attribute values (str/int/float/bool).
    If `G.graph['closeness_scan']` exists (set by build_ssn when
    auto_cutoff=True), it's a pandas.DataFrame and gets serialized to a
    JSON string before writing -- only in the file, not in your `G`
    (which is left untouched, so you can still call plot_closeness_scan
    on G.graph['closeness_scan'] afterwards as a real DataFrame).

    Examples
    --------
    >>> export_ssn(G, 'ssn.graphml')
    """
    if fmt == 'graphml':
        if 'closeness_scan' in G.graph:
            G = G.copy()
            G.graph['closeness_scan'] = G.graph['closeness_scan'].to_json()
        nx.write_graphml(G, path)
    elif fmt == 'tsv':
        nx.to_pandas_edgelist(G).to_csv(path, sep='\t', index=False)
    else:
        raise ValueError("fmt must be 'graphml' or 'tsv'")
    logger.info(f'ssn saved to {path} ({fmt})')


def to_cytoscape(G, title='rotifer_ssn', layout='force-directed', color_by=None):
    """
    Push an SSN straight into a running Cytoscape Desktop, via
    py4cytoscape/cyREST -- instead of reimplementing layout, styling
    and selection in Python.

    Requires:
      - Cytoscape open on the same machine (cyREST listens on
        localhost:1234 by default; this is automatic once Cytoscape
        is running, nothing extra to start).
      - `pip install py4cytoscape` in this environment.

    Once the SSN is in Cytoscape, use its native tools:
      - Layout menu: any algorithm, re-run anytime, no need to call
        this function again.
      - Style panel: map any column (bitscore, evalue, Louvain, Leiden,
        pident, ...) to color/size/width with a couple of clicks --
        far more flexible than a hand-coded stylesheet.
      - Select panel (or click/lasso the canvas): select nodes by
        clicking, dragging, or filtering on a column (e.g. Louvain == 3
        selects a whole community). Get the selection back into
        Python any time with get_selected_ids().

    Parameters
    ----------
    G : networkx.Graph
        Output of build_ssn / sequence_ssn / to_ssn.
    title : str, default 'rotifer_ssn'
        Name the SSN gets inside Cytoscape.
    layout : str, default 'force-directed'
        Applied once on load. Any name Cytoscape recognizes (see
        py4cytoscape.get_layout_names() for the full list) -- you can
        always pick a different one from Cytoscape's own Layout menu
        afterwards.
    color_by : str or None, default None
        If given (e.g. 'Louvain' or 'Leiden', from
        build_ssn(add_community=True)), immediately maps that
        node column to color as a discrete mapping.

    Returns
    -------
    int : the new SSN's SUID in Cytoscape.

    Examples
    --------
    >>> G, hits = rdao.sequence_ssn(seqs, add_community=True)
    >>> rdao.to_cytoscape(G, color_by='Louvain')
    """
    if G.number_of_nodes() == 0:
        raise ValueError(
            'to_cytoscape: G has 0 nodes -- nothing to push to Cytoscape. '
            'This means build_ssn\'s `criteria` matched no pairs at all '
            '(nodes only exist because an edge brought them in). Try a looser '
            'criteria, e.g. build_ssn(hits, criteria="bitscore > 0"), to confirm.'
        )
    import py4cytoscape as p4c
    suid = p4c.create_network_from_networkx(G, title=title)
    p4c.layout_network(layout)
    if color_by:
        p4c.set_node_color_mapping(**p4c.gen_node_color_map(color_by, mapping_type='d'))
    return suid


def get_selected_ids():
    """
    IDs currently selected in Cytoscape (click, lasso-select, or the
    Select panel's column filter -- e.g. filter Louvain == 3 to grab a
    whole community). Use this to collect or eliminate sequences from
    your original sequence object:

        ids = rdao.get_selected_ids()
        mantidos = seqs.df[seqs.df.id.isin(ids)]     # keep just these
        restante = seqs.df[~seqs.df.id.isin(ids)]    # drop these

    Returns
    -------
    list of str
    """
    import py4cytoscape as p4c
    return p4c.get_selected_nodes()


def closeness_scan(hits, column='bitscore', cutoffs=None, n_steps=50,
                    min_component_size=2, max_singleton_fraction=0.5,
                    auto_start=True, progress=True):
    """
    Scan a range of edge cutoffs and compute average closeness
    centrality at each one -- an objective, data-driven way to pick a
    subfamily-defining cutoff instead of guessing a number (Hornung &
    Terrapon 2023, PLOS Comput Biol, doi:10.1371/journal.pcbi.1010881).

    As the cutoff gets stricter, edges disappear and the ssn
    gradually splits into disconnected subnetworks (candidate
    subfamilies). Average closeness centrality peaks exactly when a
    split happens cleanly (no sparse straggling connections left
    between the pieces) -- so LOCAL MAXIMA in the resulting curve mark
    good candidate cutoffs, not necessarily a single "best" value.
    Multiple peaks can be legitimate (different levels of granularity,
    e.g. taxonomic vs functional splits); a curve that never rises
    much above its starting point suggests the family may not have
    real subfamilies at all -- see plot_closeness_scan.

    Two corrections applied here that the paper flags explicitly and
    that a naive implementation misses:

    1. Closeness centrality alone does NOT penalize excessive
       fragmentation -- the paper's own words: "two perfect
       subnetworks of size 50 will give the same outcome as 50 perfect
       subnetworks of size 2". Left unchecked, a cutoff that shatters
       the family into dozens of tiny cliques can score just as well
       (or better) than one clean, meaningful split. `min_component_size`
       excludes components smaller than this from the average -- they
       don't count as "good" or "bad", they're just not counted.
    2. Cutoffs past `max_singleton_fraction` fragmentation are still
       scanned and returned -- the scan always covers the matrix's real
       min-to-max range -- but they're excluded from automatic peak
       selection (see `over_fragmented` in the returned columns).
       A domain with a tight conserved core surrounded by weak,
       divergent peripheral sequences fragments early at the edges
       while the core can still hold a genuine peak much higher up; an
       early stop used to hide that peak entirely. This mirrors the
       paper's stated goal -- not handing back an over-fragmented split
       as the answer -- without the side effect of never testing the
       rest of the range.

    Parameters
    ----------
    hits : pandas.DataFrame
        Output of all_vs_all_search -- the RAW table, not yet filtered.
    column : str, default 'bitscore'
        Column to scan. Higher values are assumed stricter.
    cutoffs : iterable or None, default None
        Explicit cutoff values to test. If None, uses `n_steps` values
        evenly spaced between the column's min and max.
    n_steps : int, default 50
        Only used when `cutoffs` is None.
    min_component_size : int, default 2
        Connected components smaller than this are excluded from the
        average closeness calculation at every cutoff (see point 1
        above). Raise it if your family's meaningful subfamilies are
        never smaller than, say, 5 members, to also ignore small but
        non-singleton noise.
    max_singleton_fraction : float, default 0.5
        Cutoffs where more than this fraction of nodes have become
        singletons are excluded from automatic peak selection, but still
        scanned and returned (see point 2 above). Lower it (e.g. 0.2) to
        be stricter about what counts as "too fragmented"; raise it
        (e.g. 0.8) to allow the automatic choice to land on a stricter
        cutoff even with many singletons around it.
    auto_start : bool, default True
        Start the scan at the lowest cutoff that actually fragments the
        network, instead of at the column minimum. This is the mirror
        image of `max_singleton_fraction`: that one stops the scan once
        the network has shattered, this one skips the head of the range
        where nothing has split yet.

        The reasoning is structural, not heuristic. While the network is
        still a SINGLE connected component, no subfamily split has
        occurred, so by construction that cutoff cannot be the answer to
        the question this scan asks -- and `build_ssn`/`plot_closeness_scan`
        already refuse to treat those points as candidate peaks. Scanning
        them anyway burns `n_steps` on a region whose outcome is known in
        advance, leaving fewer points for the region that decides the
        result. The fragmentation point is located by bisection (about
        log2(n_steps) extra network builds), and all `n_steps` are then
        spent between it and the column maximum.

        Set to False to scan the raw full range (e.g. to plot the whole
        curve for a figure, including the uninformative flat head).
    progress : bool, default True
        Show a tqdm progress bar (one step per cutoff tested).

    Returns
    -------
    pandas.DataFrame
        Columns: cutoff, n_nodes, n_edges, n_components, n_unclustered,
        avg_closeness. Rows after the stopping point are not included.
        The frame carries `.attrs['scan_start']` and
        `.attrs['scan_start_reason']` describing where the scan began.

    Examples
    --------
    >>> hits = pd.read_csv('LIC_11568_ssn_hits.tsv', sep='\\t')
    >>> scan = rdao.closeness_scan(hits)
    >>> rdao.plot_closeness_scan(scan, path='LIC_11568_closeness.png')
    """
    if len(hits) == 0:
        logger.warning('closeness_scan: hits table is empty')
        return pd.DataFrame(columns=['cutoff', 'n_nodes', 'n_edges', 'n_components',
                                      'n_unclustered', 'avg_closeness'])

    total_nodes = len(set(hits['qseqid']) | set(hits['sseqid']))

    lo_raw, hi_raw = float(hits[column].min()), float(hits[column].max())
    scan_start, start_reason = lo_raw, 'column minimum'

    if cutoffs is None:
        if auto_start and hi_raw > lo_raw:
            scan_start, start_reason = _first_fragmenting_cutoff(
                hits, column, lo_raw, hi_raw, n_steps=n_steps)
            if scan_start > lo_raw:
                logger.info(
                    f'closeness_scan: starting at {column} > {scan_start:.4g} instead of '
                    f'{lo_raw:.4g} -- below that the network is still one connected '
                    f'component, so no split can be detected there ({start_reason}). '
                    f'Pass auto_start=False to scan the full raw range.'
                )
        cutoffs = np.linspace(scan_start, hi_raw, n_steps)

    rows = []
    for c in tqdm(cutoffs, desc='closeness_scan', unit='cutoff', disable=not progress):
        G = build_ssn(hits, criteria=f'{column} > {c}', progress=False)

        # Sequences with zero surviving hits simply aren't added as
        # nodes by build_ssn -- but they're just as "unclustered"
        # as a true singleton, and need to count as such here (the
        # paper explicitly adds these as singleton nodes for the same
        # bookkeeping reason). n_unclustered below counts both.
        if G.number_of_nodes() == 0:
            rows.append({'cutoff': c, 'n_nodes': 0, 'n_edges': 0, 'n_components': 0,
                         'n_unclustered': total_nodes, 'avg_closeness': np.nan,
                         'over_fragmented': True})
            continue

        components = list(nx.connected_components(G))
        nodes_in_big_components = sum(len(comp) for comp in components if len(comp) >= min_component_size)
        n_unclustered = total_nodes - nodes_in_big_components

        # Flag over-fragmented cutoffs instead of stopping the scan there:
        # a dense conserved core ringed by weak divergent sequences
        # fragments early at the edges while still holding a genuine peak
        # much higher up. Excluded from automatic choice, still reported.
        over_frag = bool(total_nodes and (n_unclustered / total_nodes) > max_singleton_fraction)

        scored_nodes = [n for comp in components if len(comp) >= min_component_size for n in comp]
        if scored_nodes:
            sub = G.subgraph(scored_nodes)
            closeness = nx.closeness_centrality(sub, wf_improved=False)
            avg_closeness = float(np.mean(list(closeness.values())))
        else:
            avg_closeness = np.nan

        rows.append({
            'cutoff': c,
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'n_components': len(components),
            'n_unclustered': n_unclustered,
            'avg_closeness': avg_closeness,
            'over_fragmented': over_frag,

        })
    out = pd.DataFrame(rows)
    out.attrs['scan_start'] = scan_start
    out.attrs['scan_start_reason'] = start_reason
    return out


def _first_fragmenting_cutoff(hits, column, lo, hi, n_steps=50, max_probes=None):
    """
    Lowest cutoff at which the network stops being a single connected
    component, found by bisection.

    Monotonicity is what makes bisection valid here: raising the cutoff
    only ever REMOVES edges, and removing edges can never merge two
    components. So "is the network still in one piece?" flips from True
    to False exactly once along the cutoff axis, and never flips back.

    Returns (cutoff, reason). Falls back to `lo` when the network is
    already fragmented at the bottom of the range (nothing to skip), and
    to `lo` as well when it never fragments (degenerate -- the caller
    should scan normally and let the usual peak filtering handle it).
    """
    import math

    def _n_components(c):
        G = build_ssn(hits, criteria=f'{column} > {c}', progress=False)
        if G.number_of_nodes() == 0:
            return 0
        return nx.number_connected_components(G)

    if _n_components(lo) != 1:
        return lo, 'already fragmented at the range minimum'
    if _n_components(hi) == 1:
        return lo, 'never fragments across the range'

    budget = max_probes if max_probes else max(6, int(math.log2(max(n_steps, 2))) + 3)
    a, b = lo, hi
    for _ in range(budget):
        mid = (a + b) / 2.0
        if _n_components(mid) == 1:
            a = mid          # ainda inteiro -> fragmentacao esta acima
        else:
            b = mid          # ja quebrou -> fragmentacao esta em ou abaixo
    return a, f'bisection over {budget} probes'


def scan_diagnostics(scan, good_closeness=0.7):
    """
    Turn a closeness scan into a short, checkable verdict about whether
    the cutoff it suggests can be trusted.

    Every field below is read straight off the scan table -- nothing is
    inferred or estimated. The point is to catch, automatically, the two
    failure modes that are easy to miss when eyeballing a single curve:

    `truncated`
        The curve was still climbing when the scan ended, and its
        maximum sits at the very edge. The real optimum is then somewhere
        ABOVE the range that was tested, so the chosen cutoff is only the
        best of what was sampled. Re-run over a wider/stricter range.
    `flat`
        The curve barely moves across the whole range. There is no peak
        to speak of, so whichever cutoff comes out is close to arbitrary
        -- often a sign the family simply has no clean subfamily
        structure at this resolution.

    `max_closeness` is worth reading alongside both: in the reference
    work a clean split lands near 0.9-1.0, so a maximum far below that
    means the pieces never separated tidily at any tested cutoff.

    Returns
    -------
    dict with keys: n_points, max_closeness, cutoff_at_max, truncated,
    flat, span, peaks, verdict (a short human-readable summary).
    """
    if scan is None or len(scan) == 0 or 'avg_closeness' not in scan:
        return {'n_points': 0, 'verdict': 'empty scan'}

    valid = scan.dropna(subset=['avg_closeness'])
    if len(valid) < 3:
        return {'n_points': len(valid), 'verdict': 'too few points to judge'}

    v = valid['avg_closeness'].to_numpy()
    c = valid['cutoff'].to_numpy()
    vmax, vmin = float(v.max()), float(v.min())
    imax = int(v.argmax())
    span = (vmax - vmin) / vmax if vmax > 0 else 0.0

    tail_rising = bool(len(v) >= 3 and v[-1] > v[-2] > v[-3])
    max_at_edge = imax >= len(v) - 2
    truncated = bool(tail_rising and max_at_edge)
    flat = bool(span < 0.10)

    peaks = []
    for i in range(1, len(v) - 1):
        if v[i] >= v[i-1] and v[i] >= v[i+1] and (v[i] > v[i-1] or v[i] > v[i+1]):
            if 'n_components' in valid.columns:
                if int(valid['n_components'].iloc[i]) < 2:
                    continue
            peaks.append(float(c[i]))

    if truncated:
        verdict = ('curve still rising at the end of the scanned range -- '
                   'the real optimum is likely above it; widen the range')
    elif flat:
        verdict = ('curve is essentially flat -- no meaningful peak, the '
                   'chosen cutoff is close to arbitrary')
    elif vmax < good_closeness:
        verdict = (f'best closeness {vmax:.2f} stays below {good_closeness:.2f} -- '
                   'components never separated cleanly at any tested cutoff')
    else:
        verdict = f'{len(peaks)} candidate peak(s), max closeness {vmax:.2f}'

    return {'n_points': len(valid), 'max_closeness': vmax, 'cutoff_at_max': float(c[imax]),
            'truncated': truncated, 'flat': flat, 'span': float(span),
            'peaks': peaks, 'verdict': verdict}


def plot_closeness_scan(scan, path=None, title=None):
    """
    Plot the closeness_scan curve (cutoff vs average closeness
    centrality) and mark its local maxima -- your candidate cutoffs.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of closeness_scan.
    path : str or None, default None
        If given, saves the figure there as PNG (e.g.
        'LIC_11568_closeness.png'). If None, the figure is returned
        without saving (handy in a notebook).
    title : str or None, default None
        Plot title -- pass the alignment's own name/prefix to tell
        multiple plots apart later.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> scan = rdao.closeness_scan(hits)
    >>> rdao.plot_closeness_scan(scan, path='LIC_11568_closeness.png', title='LIC_11568')
    """
    import matplotlib
    matplotlib.use('Agg')  # headless-safe: no X server/display needed
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))

    if scan is None or len(scan) == 0 or 'avg_closeness' not in scan.columns:
        ax.text(0.5, 0.5, 'sem dados de scan (rede degenerada ou N pequeno demais)',
                ha='center', va='center', transform=ax.transAxes, color='#888')
        ax.set_xlabel('corte')
        ax.set_ylabel('closeness centrality média')
        ax.set_title(title or 'Closeness centrality scan')
        fig.tight_layout()
        if path:
            fig.savefig(path, dpi=150)
            plt.close(fig)
            logger.info(f'plot saved to {path} (empty scan placeholder)')
        return fig

    ax.plot(scan['cutoff'], scan['avg_closeness'], '-o', markersize=3, color='#332288')

    valid = scan.dropna(subset=['avg_closeness']).reset_index(drop=True)
    if len(valid) >= 3:
        v = valid['avg_closeness'].to_numpy()
        # Local maxima, tolerant of plateaus (closeness often sits at
        # exactly 1.0 across several consecutive cutoffs once a split
        # is clean) -- marks the point where the curve rises into a
        # flat top, not just single-point spikes.
        is_peak = np.zeros(len(v), dtype=bool)
        for i in range(1, len(v) - 1):
            if v[i] >= v[i - 1] and v[i] >= v[i + 1] and (v[i] > v[i - 1] or v[i] > v[i + 1]):
                is_peak[i] = True
        if v[0] > v[1]:
            is_peak[0] = True
        if v[-1] > v[-2]:
            is_peak[-1] = True
        # Same filter as build_ssn's auto_cutoff: a high closeness
        # value with n_components == 1 is just "still one blob", not a
        # real subfamily split -- don't mark it as a candidate cutoff.
        if 'n_components' in valid.columns:
            is_peak &= (valid['n_components'] >= 2).to_numpy()
        if is_peak.any():
            ax.scatter(valid.loc[is_peak, 'cutoff'], valid.loc[is_peak, 'avg_closeness'],
                       color='#CC6677', s=60, zorder=5, label='picos locais (cortes candidatos)')
            ax.legend()

    ax.set_xlabel('corte')
    ax.set_ylabel('closeness centrality média')
    ax.set_title(title or 'Closeness centrality scan')
    fig.tight_layout()

    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
        logger.info(f'plot saved to {path}')
    return fig



# =============================================================================
# interactive SSN panel
# =============================================================================

def _dash_scan_msa_files(msa_dir):
    """
    Read every alignment file in `msa_dir`, keeping headers in original
    order (needed for display -- alignments read better in whatever order
    they were built in, not resorted). Returns [(path, [(header, seq), ...]), ...].

    Matching to a unit happens by content (does this file contain that
    unit's query header), not by filename -- these files may have been
    produced by hand, outside this pipeline's own naming convention.
    """
    import glob as _glob
    out = []
    paths = []
    for ext in ("*.msa", "*.aln", "*.afa", "*.fasta", "*.fa"):
        paths.extend(_glob.glob(os.path.join(msa_dir, ext)))
    for fp in sorted(set(paths)):
        records, header, buf = [], None, []
        try:
            with open(fp) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if header is not None:
                            records.append((header, "".join(buf)))
                        header, buf = line[1:].split()[0], []
                    else:
                        buf.append(line.strip())
            if header is not None:
                records.append((header, "".join(buf)))
        except Exception as e:
            logger.warning(f"_dash_scan_msa_files: could not read {fp}: {e}")
            continue
        if len(records) >= 2:  # a 1-sequence "alignment" is of no use here
            out.append((fp, records))
    return out


def _dash_find_reference_nodes(G, unit, aliases=None):
    """
    Locate the reference ("query") node(s) of one network.

    Returns a list, not a single id -- a unit can have more than one known
    reference sequence landing in the same network (e.g. several named
    seeds from an HMM-based search, where there's no single literal query
    sequence the way there is for a sequence-seeded search). Every
    alias/unit-name that actually exists as a node in G is kept, not just
    the first match.

    `aliases` is an optional {unit: [id, ...]} map. When a unit has exactly
    one reference (the common case), this returns a single-item list, so
    nothing downstream needs to special-case "one vs many".
    """
    candidates = list((aliases or {}).get(unit, [])) + [unit]
    found = [c for c in candidates if c in G]
    # nothing matched by name: accept a unique substring match
    if not found:
        hits = [n for n in G.nodes if unit in str(n)]
        if len(hits) == 1:
            found = hits
    # de-duplicate, preserving order
    seen = set()
    out = []
    for f in found:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def _dash_community_attr(G):
    for _, d in G.nodes(data=True):
        for a in ("Louvain", "Leiden", "leidein", "community"):
            if a in d:
                return a
        break
    return None


def _dash_build_unit(gml_path, matrix_path, aliases, max_nodes, layout_iter):
    import networkx as nx
    import base64
    name = os.path.splitext(os.path.basename(gml_path))[0]
    try:
        G = nx.read_graphml(gml_path)
    except Exception as e:
        return {"name": name, "error": f"graphml ilegivel: {e}"}

    n_total = G.number_of_nodes()
    if n_total == 0:
        return {"name": name, "error": "grafo vazio"}

    attr = _dash_community_attr(G)
    query_nodes = _dash_find_reference_nodes(G, name, aliases)
    # Every reference is marked as a query node, but the highlighted GROUP
    # follows the first one only: if a reference diverges into another
    # community that stays visible instead of being absorbed into a union.
    primary_query = query_nodes[0] if query_nodes else None
    query_comms = ({G.nodes[primary_query].get(attr)}
                   if (primary_query and attr) else set())
    query_comms.discard(None)

    # numeric value of the chosen cutoff (e.g. "bitscore > 56.2")
    auto_cutoff_raw = G.graph.get("auto_cutoff")
    auto_cutoff_val = None
    if auto_cutoff_raw:
        try:
            auto_cutoff_val = float(str(auto_cutoff_raw).split(">")[-1].strip())
        except (ValueError, IndexError):
            pass

    # Curve still rising at the right edge means the real peak lies beyond
    # the tested range (Hornung & Terrapon 2023: a good cutoff reaches
    # ~0.9-1.0; low values with a rising edge mean the scan was too short).
    scan_truncated = None
    scan_max = None
    scan_json_raw = G.graph.get("closeness_scan")
    if scan_json_raw:
        try:
            import io, pandas as _pd
            _scan = _pd.read_json(io.StringIO(scan_json_raw))
            _col = next((c for c in _scan.columns if "closeness" in c.lower()), None)
            if _col is not None and len(_scan) >= 4:
                _v = _scan[_col].dropna().tolist()
                scan_max = float(max(_v))
                # rising over the last steps AND the global max sits at the end
                tail_rising = _v[-1] > _v[-2] > _v[-3]
                max_at_edge = _v.index(max(_v)) >= len(_v) - 2
                scan_truncated = bool(tail_rising and max_at_edge)
        except Exception:
            pass

    # closeness_scan PNG, written next to the .graphml by the 'ssn' stage
    scan_png_b64 = None
    # famflow writes "{name}_closeness_scan.png"; accept "{name}.png" too
    base = os.path.splitext(gml_path)[0]
    png_path = base + "_closeness_scan.png"
    if not os.path.exists(png_path):
        png_path = base + ".png"
    if os.path.exists(png_path):
        try:
            with open(png_path, "rb") as fh:
                scan_png_b64 = base64.b64encode(fh.read()).decode("ascii")
        except Exception:
            pass

    # optional subsampling; the query's own community is always kept whole
    subsampled = False
    if max_nodes and n_total > max_nodes:
        subsampled = True
        keep = set()
        if query_comms:
            keep |= {n for n, d in G.nodes(data=True) if d.get(attr) in query_comms}
        # fill up with the most connected nodes from other communities
        resto = sorted((n for n in G.nodes if n not in keep),
                       key=lambda n: -G.degree(n))
        for n in resto:
            if len(keep) >= max_nodes:
                break
            keep.add(n)
        G = G.subgraph(keep).copy()

    # deterministic layout (fixed seed -> same figure every time)
    try:
        # Dense networks collapse into a flat blob under plain
        # spring_layout: a larger 'k' pushes nodes apart and more
        # iterations let the layout actually converge.
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        densidade = (2 * n_edges / (n_nodes * (n_nodes - 1))) if n_nodes > 1 else 0

        import math
        # networkx default k is 1/sqrt(n); denser graphs need more spread
        k = (1.0 / math.sqrt(n_nodes)) * (1.0 + 2.0 * min(densidade, 1.0))

        pos = nx.spring_layout(G, seed=42, iterations=max(layout_iter, 100), k=k)
    except Exception:
        pos = nx.random_layout(G, seed=42)

    # normalise coordinates to [0,1]
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    rangex = (maxx - minx) or 1
    rangey = (maxy - miny) or 1

    # Connected component per node. In Hornung & Terrapon (2023) a
    # "subfamily" IS the connected component -- the closeness criterion
    # looks for the cutoff at which groups come apart. Louvain partitions
    # by modularity even inside a still-connected graph.
    comp_of = {}
    for ci, comp in enumerate(nx.connected_components(G)):
        for n in comp:
            comp_of[n] = ci
    query_components = ({comp_of[primary_query]} if primary_query in comp_of else set())
    query_nodes_set = set(query_nodes)

    node_list, node_index = [], {}
    for i, n in enumerate(G.nodes):
        d = G.nodes[n]
        comm = d.get(attr) if attr else None
        node_index[n] = i
        node_list.append({
            "id": str(n),
            "x": round((pos[n][0] - minx) / rangex, 4),
            "y": round((pos[n][1] - miny) / rangey, 4),
            "c": str(comm) if comm is not None else "-",
            "cp": comp_of.get(n, -1),
            "q": bool(n in query_nodes_set),
            "qc": bool(comm in query_comms),
            "qcp": bool(comp_of.get(n, -1) in query_components),
            "deg": G.degree(n),
        })

    # Slider edges come from the RAW MATRIX, not the graph: the .graphml is
    # written with the cutoff already applied, so it holds nothing below it
    # for the slider to bring back. The graph still supplies per-node
    # attributes (Louvain community), which the browser cannot recompute.
    edge_list, edge_w = [], []
    has_weights = False
    matrix_pairs = None

    if matrix_path and os.path.exists(matrix_path):
        try:
            import pandas as pd
            _df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
            _cq = next((c for c in _df.columns if c.lower() in ("qseqid", "query", "source")), None)
            _cs = next((c for c in _df.columns if c.lower() in ("sseqid", "subject", "target")), None)
            _cb = next((c for c in _df.columns
                        if "bitscore" in c.lower() or c.lower() in ("bits", "score")), None)
            if _cq and _cs and _cb:
                matrix_pairs = _df[[_cq, _cs, _cb]].dropna()
        except Exception as e:
            logger.warning(f"could not read matrix for slider edges: {e}")

    if matrix_pairs is not None and len(matrix_pairs):
        seen = set()
        for q, s, b in matrix_pairs.itertuples(index=False):
            q, s = str(q), str(s)
            if q == s or q not in node_index or s not in node_index:
                continue
            iq, is_ = node_index[q], node_index[s]
            key = (iq, is_) if iq < is_ else (is_, iq)
            if key in seen:
                continue
            seen.add(key)
            edge_list.append([key[0], key[1]])
            edge_w.append(round(float(b), 1))
        has_weights = bool(edge_w)

    if not edge_list:
        # no usable matrix: fall back to the graph's own edges (slider is
        # then limited to what survived the cutoff; the panel says so)
        _wkey = None
        for u, v, d in G.edges(data=True):
            if _wkey is None:
                _wkey = next((k for k in d
                              if "bitscore" in k.lower() or k.lower() in ("bits", "weight")), False)
            edge_list.append([node_index[u], node_index[v]])
            edge_w.append(round(float(d.get(_wkey, 0)), 1) if _wkey else 0)
        has_weights = bool(_wkey) and any(w > 0 for w in edge_w)
        edges_from_matrix = False
    else:
        edges_from_matrix = True

    # community counts on the FULL graph, not the subsampled one
    if attr:
        G_full = nx.read_graphml(gml_path)
        from collections import Counter
        comm_counts = Counter(str(d.get(attr)) for _, d in G_full.nodes(data=True))
        comm_sizes = sorted(comm_counts.items(), key=lambda kv: -kv[1])[:15]
    else:
        comm_sizes = []

    # raw matrix statistics
    matrix_stats = None
    if matrix_path and os.path.exists(matrix_path):
        try:
            import pandas as pd
            df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
            col_bit = next((c for c in df.columns
                            if "bitscore" in c.lower() or c.lower() == "bits"), None)
            if col_bit:
                vals = df[col_bit].dropna()
                # 30-bin histogram
                import numpy as np
                counts, edges = np.histogram(vals, bins=30)
                matrix_stats = {
                    "n_pairs": len(df),
                    "min": float(vals.min()), "max": float(vals.max()),
                    "mean": float(vals.mean()), "median": float(vals.median()),
                    "hist_counts": [int(c) for c in counts],
                    "hist_edges": [round(float(e), 1) for e in edges],
                }
        except Exception as e:
            matrix_stats = {"error": str(e)}

    return {
        "name": name,
        "n_nodes_total": n_total,
        "n_nodes_shown": len(node_list),
        "n_edges_shown": len(edge_list),
        "subsampled": subsampled,
        "query_nodes": [str(q) for q in query_nodes],
        # preformatted for the header, so the browser needn't join it
        "query_node": ", ".join(str(q) for q in query_nodes) if query_nodes else None,
        "query_comm": ",".join(str(c) for c in sorted(query_comms, key=str)) if query_comms else None,
        "query_comm_size": sum(1 for n in node_list if n["qc"]),
        "n_communities": len(comm_sizes),
        "comm_sizes": comm_sizes,
        "nodes": node_list,
        "edges": edge_list,
        "matrix": matrix_stats,
        "n_components": len(set(comp_of.values())) if comp_of else 0,
        "query_component_size": sum(1 for n in node_list if n["qcp"]),
        "scan_truncated": scan_truncated,
        "scan_max": scan_max,
        "edge_w": edge_w if has_weights else None,
        "has_weights": has_weights,
        "edges_from_matrix": edges_from_matrix,
        "auto_cutoff": auto_cutoff_val,
        "auto_cutoff_raw": str(auto_cutoff_raw) if auto_cutoff_raw else None,
        "scan_png": scan_png_b64,
    }


_DASH_TEMPLATE = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<title>Painel SSN — domínios T2SS</title>
<style>
  :root {
    /* -- superficies: azul-marinho profundo, nao preto puro (traco de placa-mae) -- */
    --void: #0a0e17; --panel: #0e1420; --panel-2: #121a2b; --card: #10182799;
    --border: #1f2c42; --border-bright: #2c3f5c;

    /* -- texto -- */
    --text: #f0f6ff; --dim: #a8bcd8; --faint: #7089ab;

    /* -- duotone tecnico: ciano = dado/interacao, magenta = alerta/corte -- */
    --cyan: #4fe8d6; --cyan-dim: #2c9c8f; --cyan-glow: rgba(79,232,214,0.35);
    --magenta: #ff5f9e; --magenta-dim: #b23d6d; --magenta-glow: rgba(255,95,158,0.30);
    --amber: #ffb454;

    /* -- papeis dos nos da rede -- */
    --node-query: #f4fffc; --node-comm: #33d9c4; --node-other: #56688c;

    --mono: ui-monospace, "SF Mono", "Cascadia Code", Consolas, monospace;
    --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; font-family: var(--sans);
    background: var(--void); color: var(--text); font-size: 14px; line-height: 1.5;
    font-weight: 500;
    background-image:
      radial-gradient(ellipse 900px 500px at 15% -10%, rgba(79,232,214,0.05), transparent),
      radial-gradient(ellipse 700px 500px at 100% 10%, rgba(255,95,158,0.04), transparent);
    background-attachment: fixed;
  }
  .layout { display: flex; height: 100vh; }

  /* -------- sidebar -------- */
  aside {
    width: 264px; flex-shrink: 0; background: var(--panel-2);
    border-right: 1px solid var(--border); overflow-y: auto; padding: 18px 0;
    background-image: repeating-linear-gradient(180deg, rgba(79,232,214,0.02) 0px,
      rgba(79,232,214,0.02) 1px, transparent 1px, transparent 28px);
  }
  aside h1 {
    font-size: 14px; margin: 0 16px 2px; letter-spacing: 0.02em;
    font-family: var(--mono); color: var(--cyan); text-shadow: 0 0 12px var(--cyan-glow);
  }
  aside .sub { font-size: 11.5px; color: var(--dim); margin: 0 16px 18px;
    font-family: var(--mono); font-weight: 500; }
  .protein { margin-bottom: 3px; }
  .protein-name {
    padding: 8px 16px 5px; font-weight: 700; font-size: 11.5px;
    color: var(--cyan-dim); letter-spacing: 0.08em; font-family: var(--mono);
  }
  .domain {
    padding: 7px 16px 7px 26px; cursor: pointer; font-size: 13px;
    border-left: 2px solid transparent; transition: background 0.12s, border-color 0.12s;
    display: flex; justify-content: space-between; align-items: center; gap: 8px;
    color: var(--text); font-weight: 500;
  }
  .domain:hover { background: rgba(79,232,214,0.04); color: var(--text); }
  .domain.active {
    background: linear-gradient(90deg, rgba(79,232,214,0.10), transparent 85%);
    border-left-color: var(--cyan); color: var(--text); font-weight: 600;
  }
  .domain .badge {
    font-size: 11px; color: var(--dim); font-family: var(--mono);
    font-variant-numeric: tabular-nums; font-weight: 600;
  }
  .domain.active .badge { color: var(--cyan-dim); }
  .domain.err { color: var(--magenta-dim); font-style: italic; }

  /* -------- main -------- */
  main { flex: 1; display: flex; flex-direction: column; overflow: hidden; background: var(--void); }
  header {
    padding: 20px 26px 0; border-bottom: 1px solid var(--border); background: var(--panel);
  }
  header h2 {
    margin: 0 0 5px; font-size: 19px; letter-spacing: -0.01em; font-family: var(--mono);
    color: var(--text);
  }
  header .meta { font-size: 12.5px; color: var(--text); margin-bottom: 15px;
    font-family: var(--mono); font-weight: 600; }
  .tabs { display: flex; gap: 4px; }
  .tab {
    padding: 8px 16px; cursor: pointer; font-size: 13px; border: none;
    background: none; color: var(--dim); border-bottom: 2px solid transparent;
    font-family: var(--sans);
  }
  .tab:hover { color: var(--text); }
  .tab.active {
    color: var(--cyan); border-bottom-color: var(--cyan); font-weight: 600;
    text-shadow: 0 0 10px var(--cyan-glow);
  }

  .content { flex: 1; overflow: auto; padding: 22px 26px; }
  .pane { display: none; }
  .pane.active { display: block; }

  /* -------- rede -------- */
  #net-wrap { position: relative; background: var(--panel);
    border: 1px solid var(--border); border-radius: 6px;
    box-shadow: inset 0 0 60px rgba(79,232,214,0.02); }
  canvas { display: block; width: 100%; cursor: grab; border-radius: 6px 6px 0 0;
    touch-action: none; }
  canvas:active { cursor: grabbing; }
  .legend {
    display: flex; gap: 20px; align-items: center; padding: 11px 15px;
    border-top: 1px solid var(--border); font-size: 12.5px; color: var(--text);
    flex-wrap: wrap; font-family: var(--mono); font-weight: 600;
  }
  .legend .item { display: flex; align-items: center; gap: 7px; }
  .dot { width: 9px; height: 9px; border-radius: 50%; }
  .hint { font-size: 12px; color: var(--dim); margin: 11px 2px 0; font-family: var(--mono);
    font-weight: 500; }
  #tooltip {
    position: absolute; pointer-events: none; background: #060a12; color: var(--cyan);
    border: 1px solid var(--cyan-dim); padding: 5px 10px; border-radius: 4px; font-size: 12px;
    display: none; white-space: nowrap; z-index: 10; font-family: var(--mono);
    box-shadow: 0 0 16px rgba(79,232,214,0.25);
  }
  .warn {
    background: rgba(255,180,84,0.07); border: 1px solid rgba(255,180,84,0.35); color: var(--amber);
    padding: 10px 14px; border-radius: 6px; font-size: 12.5px; margin-bottom: 15px;
    font-family: var(--mono);
  }

  /* -------- stats -------- */
  .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 12px; margin-bottom: 24px; }
  .card { background: var(--panel); border: 1px solid var(--border);
    border-radius: 6px; padding: 14px 16px; }
  .card .label { font-size: 11px; color: var(--dim); letter-spacing: 0.05em;
    margin-bottom: 6px; font-family: var(--mono); font-weight: 600; }
  .card .value { font-size: 22px; font-weight: 600; font-variant-numeric: tabular-nums;
    letter-spacing: -0.01em; font-family: var(--mono); color: var(--text); }
  .card .value.accent { color: var(--cyan); text-shadow: 0 0 14px var(--cyan-glow); }
  .card .value.warn { color: var(--magenta); text-shadow: 0 0 14px var(--magenta-glow); }
  table { border-collapse: collapse; width: 100%; font-size: 13px; }
  th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--border); }
  th { font-weight: 700; font-size: 11.5px; color: var(--text); letter-spacing: 0.05em;
    font-family: var(--mono); }
  td { font-family: var(--mono); font-size: 12.5px; color: var(--text); font-weight: 500; }
  td.num { text-align: right; font-variant-numeric: tabular-nums; }
  tr.is-query td { background: rgba(79,232,214,0.06); font-weight: 600; color: var(--cyan); }
  tr.is-cutoff td { background: rgba(255,95,158,0.08); }
  h3 { font-size: 13.5px; margin: 24px 0 11px; color: var(--text); font-family: var(--mono);
    letter-spacing: 0.03em; font-weight: 700; }
  .bar { height: 7px; background: rgba(86,104,140,0.35); border-radius: 2px; }
  .bar.q { background: linear-gradient(90deg, var(--cyan-dim), var(--cyan)); box-shadow: 0 0 8px var(--cyan-glow); }
  .bar.cutoff-bin { background: linear-gradient(90deg, var(--magenta-dim), var(--magenta));
    box-shadow: 0 0 8px var(--magenta-glow); }
  .empty { color: var(--faint); font-style: italic; padding: 34px 0; text-align: center;
    font-family: var(--mono); }
  #cutoff-bar {
    display: flex; align-items: center; gap: 14px; padding: 12px 15px;
    background: var(--panel); border: 1px solid var(--border);
    border-radius: 6px; margin-bottom: 12px; font-family: var(--mono);
    font-size: 12.5px; flex-wrap: wrap;
  }
  #cutoff-bar label { color: var(--text); font-weight: 600; white-space: nowrap; }
  #cutoff-slider {
    flex: 1; min-width: 220px; accent-color: var(--magenta); cursor: pointer;
    height: 6px; border-radius: 3px; -webkit-appearance: none; appearance: none;
  }
  /* a trilha nativa do slider ignora o "background" do elemento pai por
     padrao -- precisa herdar explicitamente, senao o sombreamento da
     faixa sem dado real (aplicado via JS inline) nunca aparece */
  #cutoff-slider::-webkit-slider-runnable-track { background: inherit; height: 6px; border-radius: 3px; }
  #cutoff-slider::-moz-range-track { background: inherit; height: 6px; border-radius: 3px; }
  #cutoff-slider::-webkit-slider-thumb {
    -webkit-appearance: none; width: 14px; height: 14px; border-radius: 50%;
    background: var(--magenta); margin-top: -4px; cursor: pointer;
    box-shadow: 0 0 6px var(--magenta-glow);
  }
  #cutoff-slider::-moz-range-thumb {
    width: 14px; height: 14px; border-radius: 50%; border: none;
    background: var(--magenta); cursor: pointer; box-shadow: 0 0 6px var(--magenta-glow);
  }
  #cutoff-val {
    color: var(--magenta); font-weight: 700; min-width: 62px;
    text-shadow: 0 0 10px var(--magenta-glow); font-variant-numeric: tabular-nums;
  }
  #cutoff-stats { color: var(--dim); font-weight: 500; }
  #cutoff-stats b { color: var(--cyan); font-weight: 700; }
  .reset-btn {
    background: transparent; border: 1px solid var(--border-bright); color: var(--dim);
    padding: 4px 10px; border-radius: 4px; font-family: var(--mono); font-size: 11px;
    cursor: pointer; font-weight: 600;
  }
  .reset-btn:hover { color: var(--text); border-color: var(--cyan-dim); }
  .use-btn {
    background: rgba(79,232,214,0.10); border: 1px solid var(--cyan-dim); color: var(--cyan);
    padding: 4px 10px; border-radius: 4px; font-family: var(--mono); font-size: 11px;
    cursor: pointer; font-weight: 700;
  }
  .use-btn:hover { background: rgba(79,232,214,0.20); }
  .domain.overridden { border-left-color: var(--magenta) !important; }
  .domain.overridden .badge::after { content: ' ✓'; color: var(--magenta); }
  .mode-btn {
    background: transparent; border: 1px solid var(--border-bright); color: var(--cyan);
    padding: 4px 11px; border-radius: 4px; font-family: var(--mono); font-size: 11.5px;
    cursor: pointer; font-weight: 600;
  }
  .mode-btn:hover { background: rgba(79,232,214,0.08); border-color: var(--cyan-dim); }
  .scan-img { max-width: 100%; border: 1px solid var(--border); border-radius: 6px;
    background: #fff; padding: 6px; margin-top: 4px; }
  .hist-wrap { background: var(--panel); border: 1px solid var(--border);
    border-radius: 6px; padding: 14px 16px 4px; }
  .hist-wrap rect { transition: opacity 0.1s; cursor: default; }
  .hist-wrap rect:hover { opacity: 0.75; }

  .msa-wrap { background: var(--panel); border: 1px solid var(--border);
    border-radius: 6px; overflow-x: auto; max-height: 65vh; overflow-y: auto; }
  .msa-row { display: flex; border-bottom: 1px solid var(--border); white-space: nowrap;
    transition: opacity 0.15s, filter 0.15s; }
  .msa-row.out-of-group { opacity: 0.22; filter: grayscale(1); }
  .msa-row.msa-query { background: rgba(79,232,214,0.06); }
  .msa-label { flex-shrink: 0; width: 200px; padding: 3px 10px; font-family: var(--mono);
    font-size: 12px; color: var(--dim); border-right: 1px solid var(--border);
    position: sticky; left: 0; background: var(--panel); overflow: hidden;
    text-overflow: ellipsis; }
  .msa-row.msa-query .msa-label { color: var(--cyan); font-weight: 700; }
  .msa-seq { font-family: var(--mono); font-size: 12px; line-height: 1.5;
    letter-spacing: 0.5px; padding: 3px 6px; }
  .msa-seq span { padding: 0 0.5px; }
</style>
</head>
<body>
<div class="layout">
  <aside>
    <h1>Painel SSN</h1>
    <p class="sub">__NDOM__ redes</p>
    <div id="nav"></div>
  </aside>
  <main>
    <header>
      <h2 id="title">—</h2>
      <div class="meta" id="meta"></div>
      <div class="tabs">
        <button class="tab active" data-pane="ssn">Rede (SSN)</button>
        <button class="tab" data-pane="matrix">Matriz</button>
        <button class="tab" data-pane="stats">Comunidades</button>
        <button class="tab" data-pane="msa" id="tab-msa" style="display:none">Alinhamento</button>
        <button id="export-cutoffs-btn" class="use-btn" style="margin-left:auto;display:none"
                onclick="exportManualCutoffs()">baixar cortes manuais</button>
      </div>
    </header>
    <div class="content">
      <div class="pane active" id="pane-ssn">
        <div id="subsample-warn"></div>
        <div id="cutoff-bar"></div>
        <div id="net-wrap">
          <canvas id="net"></canvas>
          <div id="tooltip"></div>
          <div class="legend">
            <div class="item"><span class="dot" style="background:#f4fffc;box-shadow:0 0 6px #4fe8d6"></span> query</div>
            <div class="item"><span class="dot" style="background:#33d9c4"></span> grupo da query</div>
            <div class="item"><span class="dot" style="background:#56688c"></span> outras comunidades</div>
            <div class="item" id="legend-counts"></div>
            <div class="item" style="margin-left:auto">
              <span id="louvain-status" style="color:var(--faint);font-size:11px;font-family:var(--mono)"></span>
              <button id="group-mode-btn" class="mode-btn" title="equivalente ao Shift+arraste, para telas de toque">modo: nó único</button>
              <button id="mode-btn" class="mode-btn">agrupar por: componente conectado</button>
            </div>
          </div>
        </div>
        <p class="hint">Arraste um nó para movê-lo · <b>Shift + arrastar</b> (ou o botão "modo: grupo inteiro") move a comunidade inteira · arraste o fundo para deslocar · rolagem ou pinça para zoom</p>
        <div id="scan-section"></div>
      </div>
      <div class="pane" id="pane-matrix"></div>
      <div class="pane" id="pane-stats"></div>
      <div class="pane" id="pane-msa"></div>
    </div>
  </main>
</div>

<script>
const DATA = __DATA__;

// ---------- navegacao ----------
const nav = document.getElementById('nav');
const byProtein = {};
DATA.forEach((d, i) => {
  const prot = d.name.includes('_dominio_') ? d.name.split('_dominio_')[0]
    : (d.name.includes('__') ? d.name.split('__')[0] : d.name);
  (byProtein[prot] = byProtein[prot] || []).push({d, i});
});

Object.keys(byProtein).sort().forEach(prot => {
  const wrap = document.createElement('div');
  wrap.className = 'protein';
  const h = document.createElement('div');
  h.className = 'protein-name';
  h.textContent = prot;
  wrap.appendChild(h);
  byProtein[prot].forEach(({d, i}) => {
    const el = document.createElement('div');
    el.className = 'domain' + (d.error ? ' err' : '');
    el.dataset.idx = i;
    el.dataset.name = d.name;
    const num = d.name.includes('_dominio_') ? 'domínio ' + d.name.split('_dominio_')[1]
      : (d.name.includes('__') ? d.name.split('__').slice(1).join('__') : d.name);
    el.innerHTML = '<span>' + num + '</span>' +
      (d.error ? '' : '<span class="badge">' + d.n_nodes_total + '</span>');
    if (!d.error) el.onclick = () => select(i);
    wrap.appendChild(el);
  });
  nav.appendChild(wrap);
});

// ---------- abas ----------
document.querySelectorAll('.tab').forEach(t => {
  t.onclick = () => {
    document.querySelectorAll('.tab').forEach(x => x.classList.remove('active'));
    document.querySelectorAll('.pane').forEach(x => x.classList.remove('active'));
    t.classList.add('active');
    document.getElementById('pane-' + t.dataset.pane).classList.add('active');
    if (t.dataset.pane === 'ssn') { resize(); draw(); }
    if (t.dataset.pane === 'msa') { renderMsa(); }
  };
});

// ---------- estado ----------
let cur = null;
let view = {scale: 1, ox: 0, oy: 0};
const canvas = document.getElementById('net');
const ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');

function select(i) {
  cur = DATA[i];
  view = {scale: 1, ox: 0, oy: 0};
  document.querySelectorAll('.domain').forEach(e => e.classList.remove('active'));
  const el = document.querySelector('.domain[data-idx="' + i + '"]');
  if (el) el.classList.add('active');

  document.getElementById('title').textContent = cur.name;
  document.getElementById('meta').textContent =
    cur.n_nodes_total.toLocaleString('pt-BR') + ' nós · ' +
    cur.n_edges_shown.toLocaleString('pt-BR') + ' arestas exibidas · ' +
    cur.n_communities + ' comunidades' +
    (cur.query_node ? ' · query: ' + cur.query_node : ' · query ausente do grafo');

  document.getElementById('subsample-warn').innerHTML = cur.subsampled
    ? '<div class="warn">Rede grande (' + cur.n_nodes_total.toLocaleString('pt-BR') +
      ' nós). Exibindo ' + cur.n_nodes_shown.toLocaleString('pt-BR') +
      ' — a comunidade da query está inteira; as demais foram amostradas pelos nós mais conectados.</div>'
    : '';

  updateLegendCounts();

  const msaTab = document.getElementById('tab-msa');
  if (cur.msa) {
    msaTab.style.display = 'inline-block';
  } else {
    msaTab.style.display = 'none';
    if (msaTab.classList.contains('active')) {
      // dominio sem alinhamento, mas a aba de alinhamento estava aberta
      // -- volta pra rede em vez de deixar a aba ativa vazia
      msaTab.classList.remove('active');
      document.querySelector('.tab[data-pane="ssn"]').classList.add('active');
      document.getElementById('pane-msa').classList.remove('active');
      document.getElementById('pane-ssn').classList.add('active');
    }
  }

  renderMatrix();
  renderStats();
  renderScanSection();
  renderMsa();
  startSim();
}

// ---------- simulacao de forca (viva, arrastavel) ----------
// Layout calculado no navegador em vez de estatico: permite arrastar nos
// e comunidades igual no Cytoscape. Repulsao usa quadtree (Barnes-Hut)
// pra aguentar milhares de nos sem travar.

let sim = null;
// 'louvain' = comunidade por modularidade; 'component' = componente conectado,
// que e a definicao de subfamilia usada no criterio de closeness centrality
// (Hornung & Terrapon 2023)
let groupMode = 'component';

// cortes escolhidos a mao, por dominio: {nome: valor}. Guardado tambem no
// localStorage do navegador, entao sobrevive a fechar e reabrir o painel
// (mas NAO viaja com o arquivo -- e por isso que existe o botao de
// exportar, para levar a escolha de volta ao pipeline).
let manualCutoffs = {};
try {
  manualCutoffs = JSON.parse(localStorage.getItem('ssn_manual_cutoffs') || '{}');
} catch (e) { manualCutoffs = {}; }

function persistManualCutoffs() {
  try { localStorage.setItem('ssn_manual_cutoffs', JSON.stringify(manualCutoffs)); }
  catch (e) { /* localStorage indisponivel (ex: file:// em alguns navegadores) -- tudo bem, so nao persiste entre sessoes */ }
}

function markDomainOverridden(name) {
  persistManualCutoffs();
  const el = document.querySelector('.domain[data-name="' + name + '"]');
  if (el) el.classList.add('overridden');
}

function updateExportBtn() {
  const n = Object.keys(manualCutoffs).length;
  const btn = document.getElementById('export-cutoffs-btn');
  if (!btn) return;
  btn.textContent = n > 0 ? 'baixar cortes manuais (' + n + ')' : 'baixar cortes manuais';
  btn.style.display = n > 0 ? 'inline-block' : 'none';
}

function exportManualCutoffs() {
  const blob = new Blob([JSON.stringify(manualCutoffs, null, 2)], {type: 'application/json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'cortes_manuais.json';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function resize() {
  const wrap = document.getElementById('net-wrap');
  const w = wrap.clientWidth;
  const h = Math.max(480, Math.min(720, window.innerHeight - 300));
  const dpr = window.devicePixelRatio || 1;
  canvas.width = w * dpr; canvas.height = h * dpr;
  canvas.style.height = h + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {w, h};
}

// --- quadtree para repulsao O(n log n) ---
function buildQuad(nodes, x0, y0, x1, y1) {
  const root = {x0, y0, x1, y1, mass: 0, cx: 0, cy: 0, node: null, kids: null};
  function insert(q, n, depth) {
    q.mass++; q.cx += n.x; q.cy += n.y;
    if (depth > 18) { return; }
    if (!q.kids && q.node === null) { q.node = n; return; }
    if (!q.kids) {
      const old = q.node; q.node = null;
      const mx = (q.x0+q.x1)/2, my = (q.y0+q.y1)/2;
      q.kids = [
        {x0:q.x0,y0:q.y0,x1:mx,y1:my,mass:0,cx:0,cy:0,node:null,kids:null},
        {x0:mx,y0:q.y0,x1:q.x1,y1:my,mass:0,cx:0,cy:0,node:null,kids:null},
        {x0:q.x0,y0:my,x1:mx,y1:q.y1,mass:0,cx:0,cy:0,node:null,kids:null},
        {x0:mx,y0:my,x1:q.x1,y1:q.y1,mass:0,cx:0,cy:0,node:null,kids:null}
      ];
      if (old) insert(q.kids[quadIdx(q, old)], old, depth+1);
    }
    insert(q.kids[quadIdx(q, n)], n, depth+1);
  }
  function quadIdx(q, n) {
    const mx = (q.x0+q.x1)/2, my = (q.y0+q.y1)/2;
    return (n.x >= mx ? 1 : 0) + (n.y >= my ? 2 : 0);
  }
  for (const n of nodes) insert(root, n, 0);
  return root;
}

function applyRepulsion(q, n, theta, strength) {
  if (!q || q.mass === 0) return;
  const cx = q.cx / q.mass, cy = q.cy / q.mass;
  let dx = cx - n.x, dy = cy - n.y;
  let d2 = dx*dx + dy*dy;
  if (d2 < 0.01) d2 = 0.01;
  const width = q.x1 - q.x0;
  if (!q.kids || (width * width / d2) < theta * theta) {
    if (q.node === n) return;
    const f = -strength * q.mass / d2;
    const d = Math.sqrt(d2);
    n.vx += (dx/d) * f; n.vy += (dy/d) * f;
    return;
  }
  for (const k of q.kids) applyRepulsion(k, n, theta, strength);
}

function startSim() {
  if (!cur || cur.error) return;
  const {w, h} = resize();

  // posicoes iniciais vindas do layout pre-calculado (converge bem mais
  // rapido do que comecar aleatorio)
  const N = cur.nodes.map(n => ({
    id: n.id, c: n.c, cp: n.cp, q: n.q, qc: n.qc, qcp: n.qcp, deg: n.deg,
    x: n.x * w, y: n.y * h, vx: 0, vy: 0, fixed: false
  }));
  const E = cur.edges.map(([a,b]) => [N[a], N[b]]);
  // indices crus, para o union-find do corte ao vivo (E guarda
  // referencias aos objetos, que a fisica precisa)
  const Eidx = cur.edges;

  const n_nodes = N.length;

  // A repulsao e somada sobre a massa do quadtree, entao a forca total num
  // no cresce com o NUMERO de nos. Manter `strength` quase constante fazia
  // uma rede de 3000 nos receber ~100x a forca de uma de 30: o no era
  // arremessado longe, puxado de volta, e ficava oscilando sem assentar.
  // Normalizando por n, a forca por no fica comparavel em qualquer tamanho.
  const repulsion = (260 + 900 / Math.sqrt(n_nodes)) * (60 / Math.max(60, n_nodes));
  const linkDist = Math.max(18, 60 / Math.sqrt(Math.max(1, n_nodes/100)));

  // Rede grande precisa de mais amortecimento e de esfriar mais rapido:
  // sao muitos corpos interagindo, e o layout util aparece bem antes da
  // convergencia perfeita. Rede pequena pode relaxar devagar sem tremer.
  const big = n_nodes > 800;
  const damping = big ? 0.55 : 0.72;      // fracao da velocidade mantida
  const decay = big ? 0.965 : 0.985;      // quanto alpha cai por quadro
  const maxStep = big ? 8 : 20;           // teto de deslocamento por quadro

  sim = {N, E, Eidx, activeEdges: E, alpha: 1, w, h,
         repulsion, linkDist, damping, decay, maxStep, running: true};
  setupCutoffBar();
  tick();
}

function tick() {
  if (!sim || !sim.running) return;
  const {N, E, w, h} = sim;
  const alpha = sim.alpha;

  // repulsao (Barnes-Hut)
  const quad = buildQuad(N, 0, 0, w, h);
  for (const n of N) applyRepulsion(quad, n, 0.9, sim.repulsion * alpha);

  // atracao pelas arestas (so as que passam o corte atual)
  for (const [a, b] of (sim.activeEdges || E)) {
    let dx = b.x - a.x, dy = b.y - a.y;
    let d = Math.sqrt(dx*dx + dy*dy) || 0.01;
    const f = (d - sim.linkDist) / d * 0.35 * alpha;
    const fx = dx * f, fy = dy * f;
    a.vx += fx; a.vy += fy;
    b.vx -= fx; b.vy -= fy;
  }

  // centralizacao suave
  const cx = w/2, cy = h/2;
  for (const n of N) {
    n.vx += (cx - n.x) * 0.006 * alpha;
    n.vy += (cy - n.y) * 0.006 * alpha;
  }

  // integracao
  const damp = sim.damping, maxStep = sim.maxStep;
  for (const n of N) {
    if (n.fixed) { n.vx = 0; n.vy = 0; continue; }
    n.vx *= damp; n.vy *= damp;
    // teto de deslocamento por quadro: e o que impede um unico empurrao
    // grande de arremessar o no para longe e disparar o vai-e-vem. Corta
    // a oscilacao sem mudar para onde o layout converge, so o caminho.
    const sp = Math.hypot(n.vx, n.vy);
    if (sp > maxStep) { const k = maxStep / sp; n.vx *= k; n.vy *= k; }
    n.x += n.vx; n.y += n.vy;
    n.x = Math.max(8, Math.min(w-8, n.x));
    n.y = Math.max(8, Math.min(h-8, n.y));
  }

  sim.alpha *= sim.decay;
  draw();

  if (sim.alpha > 0.005) {
    requestAnimationFrame(tick);
  } else {
    sim.running = false;
  }
}

function reheat(a) {
  if (!sim) return;
  sim.alpha = Math.max(sim.alpha, a || 0.35);
  if (!sim.running) { sim.running = true; tick(); }
}

function draw() {
  if (!sim) { ctx.clearRect(0,0,canvas.width,canvas.height); return; }
  const {N, E, w, h} = sim;
  ctx.clearRect(0, 0, w, h);
  ctx.save();
  ctx.translate(view.ox, view.oy);
  ctx.scale(view.scale, view.scale);

  // arestas mais visiveis: opacidade cai com a densidade, pra hairball nao
  // virar uma mancha solida, mas rede esparsa mostrar bem as conexoes
  const nActive = (sim.activeEdges || E).length;
  const edgeAlpha = Math.max(0.10, Math.min(0.55, 900 / Math.max(1, nActive)));
  ctx.strokeStyle = 'rgba(120,220,255,' + edgeAlpha + ')';
  ctx.lineWidth = 0.9 / view.scale;
  ctx.beginPath();
  const drawEdges = sim.activeEdges || E;
  for (const [a, b] of drawEdges) { ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); }
  ctx.stroke();

  // raio escala com o tamanho da rede -- fixo em 7/4.2/3px, uma rede de
  // 3000 nos virava uma bolota so, com os circulos se sobrepondo. O piso
  // (rMin) evita que fique pequeno demais para clicar/tocar em redes
  // enormes.
  const rScale = Math.max(0.28, Math.min(1, 26 / Math.sqrt(N.length)));
  const inGroup = n => groupMode === 'component'
    ? (n.qcpLive !== undefined ? n.qcpLive : n.qcp)
    : (n.qcLive !== undefined ? n.qcLive : n.qc);
  const order = N.map((n,i)=>i)
    .sort((a,b) => (N[a].q?2:inGroup(N[a])?1:0) - (N[b].q?2:inGroup(N[b])?1:0));
  for (const i of order) {
    const n = N[i];
    let r, fill, stroke = null;
    if (n.q)          { r = 7 * rScale; fill = '#f4fffc'; stroke = '#4fe8d6'; }
    else if (inGroup(n)) { r = 4.2 * rScale; fill = '#33d9c4'; }
    else              { r = 3 * rScale; fill = '#56688c'; }
    ctx.beginPath();
    if (n.q) { ctx.shadowColor = '#4fe8d6'; ctx.shadowBlur = 14 / view.scale; }
    ctx.arc(n.x, n.y, r / Math.sqrt(view.scale), 0, 6.2832);
    ctx.fillStyle = fill; ctx.fill();
    ctx.shadowBlur = 0;
    if (stroke) { ctx.lineWidth = 2/view.scale; ctx.strokeStyle = stroke; ctx.stroke(); }
  }
  ctx.restore();
}

// --- corte ao vivo: refaz os componentes conectados a cada ajuste ---
// Reproduz no navegador o que o criterio de closeness centrality otimiza
// (Hornung & Terrapon 2023): subfamilia = componente conectado. Union-find
// e barato o bastante pra rodar a cada movimento do slider.
// --- Louvain ao vivo ---
// Mesmo algoritmo que o backend usa (o atributo gravado no .graphml se
// chama "Louvain"), rodando aqui pra recalcular quando o corte muda --
// o backend so calcula uma vez, no corte automatico; isso significa que
// ao mexer no slider, a cor por Louvain ficava presa no corte antigo,
// so o modo "componente conectado" reagia. Debounced (nao roda a cada
// pixel arrastado) porque e mais pesado que o union-find dos componentes.
function louvainCompute(N, activeEdges, edgeWeights) {
  const n = N.length;
  let comm = new Int32Array(n);
  for (let i = 0; i < n; i++) comm[i] = i;

  // grau ponderado de cada no e peso total (2m), a partir das arestas ATIVAS
  const deg = new Float64Array(n);
  let m2 = 0;
  const adj = Array.from({length: n}, () => []);  // [[vizinho, peso], ...]
  for (let k = 0; k < activeEdges.length; k++) {
    const [a, b] = activeEdges[k];
    const w = edgeWeights ? edgeWeights[k] : 1;
    deg[a] += w; deg[b] += w; m2 += 2 * w;
    adj[a].push([b, w]); adj[b].push([a, w]);
  }
  if (m2 === 0) return comm;  // sem aresta, cada um na sua

  // fase 1: mover nos entre comunidades vizinhas pra maximizar ganho de
  // modularidade, ate nao ter mais melhora (ou limite de iteracoes)
  const commWeight = new Float64Array(n);  // soma dos graus dos nos de cada comunidade
  for (let i = 0; i < n; i++) commWeight[comm[i]] += deg[i];

  let improved = true, pass = 0;
  while (improved && pass < 12) {
    improved = false; pass++;
    for (let i = 0; i < n; i++) {
      const ci = comm[i];
      const neighComm = new Map();  // comunidade vizinha -> peso da aresta pra ela
      for (const [j, w] of adj[i]) {
        const cj = comm[j];
        neighComm.set(cj, (neighComm.get(cj) || 0) + w);
      }
      commWeight[ci] -= deg[i];

      let best = ci, bestGain = neighComm.has(ci)
        ? neighComm.get(ci) - deg[i] * commWeight[ci] / m2 : -Infinity;
      for (const [cj, wij] of neighComm) {
        if (cj === ci) continue;
        const gain = wij - deg[i] * commWeight[cj] / m2;
        if (gain > bestGain) { bestGain = gain; best = cj; }
      }
      commWeight[best] += deg[i];
      if (best !== ci) { comm[i] = best; improved = true; }
    }
  }

  // renumera 0..k-1 pra ficar limpo
  const relabel = new Map();
  const out = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    if (!relabel.has(comm[i])) relabel.set(comm[i], relabel.size);
    out[i] = relabel.get(comm[i]);
  }
  return out;
}

function recomputeComponents(cutoff) {
  if (!sim) return;
  const {N, E} = sim;
  const W = cur.edge_w;
  const parent = N.map((_, i) => i);
  function find(a) { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; }
  function union(a, b) { const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }

  const Ei = sim.Eidx;
  sim.activeEdges = [];
  sim.activeEdgeIdx = [];   // pares de INDICES (nao objetos) -- o Louvain precisa disso
  sim.activeEdgeW = [];     // pesos na mesma ordem/filtro
  for (let i = 0; i < Ei.length; i++) {
    if (!W || W[i] >= cutoff) {
      sim.activeEdges.push(E[i]);
      sim.activeEdgeIdx.push(Ei[i]);
      sim.activeEdgeW.push(W ? W[i] : 1);
      union(Ei[i][0], Ei[i][1]);
    }
  }
  for (let i = 0; i < N.length; i++) N[i].cpLive = find(i);

  // com varias referencias, o brilho de "query" fica em todas, mas o
  // GRUPO destacado ao vivo segue so a primeira (mesma logica do backend)
  const firstQ = N.find(n => n.q);
  const qRoot = firstQ ? find(N.indexOf(firstQ)) : null;
  for (const n of N) n.qcpLive = (qRoot !== null && n.cpLive === qRoot);

  const roots = new Set(N.map(n => n.cpLive));
  sim.nComponents = roots.size;
  sim.qCompSize = N.filter(n => n.qcpLive).length;
}

function setupCutoffBar() {
  const bar = document.getElementById('cutoff-bar');
  if (!cur || cur.error || !cur.has_weights || !cur.edge_w) {
    bar.innerHTML = '<span style="color:var(--faint)">Sem pesos de aresta neste grafo — ' +
      'corte interativo indisponível.</span>';
    return;
  }
  const W = cur.edge_w;
  // faixa FIXA (30-200), igual em todos os dominios -- nao adapta ao
  // min/max real de cada rede. Arestas fora dessa janela simplesmente
  // nunca somem/aparecem pelo slider (ficam sempre ativas se < 30,
  // sempre inativas se > 200); o aviso abaixo avisa quando isso acontece.
  const FIXED_LO = 30, FIXED_HI = 200;
  const lo = FIXED_LO, hi = FIXED_HI;
  const realLo = Math.floor(Math.min(...W)), realHi = Math.ceil(Math.max(...W));
  const start = (cur.auto_cutoff !== null && cur.auto_cutoff !== undefined)
    ? Math.max(lo, Math.min(hi, cur.auto_cutoff)) : lo;

  // faixa real sempre mostrada, nao so quando "estoura" a regua fixa --
  // o caso oposto (dado real nao alcanca 200, tipo maximo real = 100)
  // e igual confuso e antes nao disparava aviso nenhum: o slider deixava
  // arrastar ate 200 sem avisar que dali pra cima nao existe nada.
  let rangeWarn = '';
  if (cur.edges_from_matrix === false) {
    rangeWarn = '<span style="color:var(--amber)" title="as arestas vieram do grafo, ' +
      'que já foi gravado com o corte aplicado — baixar o controle não traz arestas novas">' +
      '⚠ sem matriz: só dá para SUBIR o corte</span>';
  } else {
    const foraDaRegua = realLo < FIXED_LO || realHi > FIXED_HI;
    const naoAlcancaRegua = realHi < FIXED_HI - 1 || realLo > FIXED_LO + 1;
    if (foraDaRegua) {
      rangeWarn = '<span style="color:var(--amber)" title="bitscore real desta rede: ' +
        realLo + '–' + realHi + '">⚠ dados fora de [' + FIXED_LO + ',' + FIXED_HI + ']</span>';
    } else if (naoAlcancaRegua) {
      rangeWarn = '<span style="color:var(--amber)" title="acima/abaixo disso não existe ' +
        'nenhum par nesta rede -- o slider deixa arrastar, mas o resultado fica vazio">' +
        '⚠ dados reais só entre ' + realLo + '–' + realHi + '</span>';
    } else {
      rangeWarn = '<span style="color:var(--faint)">dados: ' + realLo + '–' + realHi + '</span>';
    }
  }

  // sombreia no proprio slider a faixa fora do intervalo real -- pista
  // visual, sem precisar ler o texto do aviso
  const pctLo = ((Math.max(lo, realLo) - lo) / (hi - lo)) * 100;
  const pctHi = ((Math.min(hi, realHi) - lo) / (hi - lo)) * 100;
  const sliderBg = 'background:linear-gradient(90deg,' +
    '#1a2333 0%,#1a2333 ' + pctLo.toFixed(1) + '%,' +
    '#28405a ' + pctLo.toFixed(1) + '%,#28405a ' + pctHi.toFixed(1) + '%,' +
    '#1a2333 ' + pctHi.toFixed(1) + '%,#1a2333 100%);';

  bar.innerHTML =
    '<label>corte (bitscore)</label>' +
    '<input type="range" id="cutoff-slider" min="' + lo + '" max="' + hi + '" ' +
      'step="' + Math.max(0.1, ((hi-lo)/300).toFixed(2)) + '" value="' + start + '" ' +
      'style="' + sliderBg + '">' +
    '<span id="cutoff-val">' + start.toFixed(1) + '</span>' +
    '<span id="cutoff-stats"></span>' + rangeWarn +
    '<button class="reset-btn" id="cutoff-reset">voltar ao auto_cutoff</button>' +
    '<button class="use-btn" id="cutoff-use">usar este corte ✓</button>';

  const slider = document.getElementById('cutoff-slider');
  slider.oninput = () => applyCutoff(parseFloat(slider.value));
  document.getElementById('cutoff-reset').onclick = () => {
    const v = (cur.auto_cutoff !== null && cur.auto_cutoff !== undefined) ? cur.auto_cutoff : lo;
    slider.value = v;
    applyCutoff(v);
    delete manualCutoffs[cur.name];
    updateExportBtn();
  };
  document.getElementById('cutoff-use').onclick = () => {
    manualCutoffs[cur.name] = parseFloat(slider.value);
    updateExportBtn();
    markDomainOverridden(cur.name);
  };
  applyCutoff(start);
  if (cur.name in manualCutoffs) {
    slider.value = manualCutoffs[cur.name];
    applyCutoff(manualCutoffs[cur.name]);
  }
}

let louvainDebounce = null;

function applyCutoff(v) {
  document.getElementById('cutoff-val').textContent = v.toFixed(1);
  recomputeComponents(v);
  const auto = cur.auto_cutoff;
  const diff = (auto !== null && auto !== undefined && Math.abs(v - auto) > 0.05)
    ? '  (auto_cutoff: ' + auto.toFixed(1) + ')' : '  ← auto_cutoff';
  document.getElementById('cutoff-stats').innerHTML =
    sim.activeEdges.length.toLocaleString('pt-BR') + ' de ' +
    cur.edge_w.length.toLocaleString('pt-BR') + ' arestas · <b>' +
    sim.nComponents.toLocaleString('pt-BR') + '</b> componentes · grupo da query: <b>' +
    sim.qCompSize.toLocaleString('pt-BR') + '</b> nós' + diff;
  reheat(0.25);
  draw();
  updateMsaDimming();

  // Louvain e mais pesado que o union-find dos componentes -- recalcula
  // com atraso (debounce), so depois que o usuario para de mexer no
  // slider, pra nao travar o arrasto.
  clearTimeout(louvainDebounce);
  louvainDebounce = setTimeout(() => {
    const t0 = performance.now();
    const louvainResult = louvainCompute(sim.N, sim.activeEdgeIdx, sim.activeEdgeW);
    for (let i = 0; i < sim.N.length; i++) sim.N[i].louvainLive = louvainResult[i];
    const firstQIdx = sim.N.findIndex(n => n.q);
    if (firstQIdx >= 0) {
      const qComm = louvainResult[firstQIdx];
      for (let i = 0; i < sim.N.length; i++)
        sim.N[i].qcLive = (louvainResult[i] === qComm);
    }
    const dt = (performance.now() - t0).toFixed(0);
    const badge = document.getElementById('louvain-status');
    if (badge) badge.textContent = 'prévia Louvain (' + dt + 'ms) — clique "usar este corte" p/ resultado definitivo';
    draw();
  }, 250);
}

// --- interacao: arrastar no (ou comunidade inteira com Shift), pan, zoom ---
function toWorld(px, py) {
  return {x: (px - view.ox) / view.scale, y: (py - view.oy) / view.scale};
}

function nodeAt(px, py) {
  if (!sim) return null;
  const p = toWorld(px, py);
  const r2 = (12 / view.scale) ** 2;
  let best = null, bestD = r2;
  for (const n of sim.N) {
    const d = (n.x-p.x)**2 + (n.y-p.y)**2;
    if (d < bestD) { bestD = d; best = n; }
  }
  return best;
}

let dragNode = null, dragGroup = null, panStart = null;
let moveGroupMode = false;  // equivalente touch do Shift+arraste

function shouldMoveGroup(e) {
  return moveGroupMode || (e && e.shiftKey);
}

function dragTargets(n, e) {
  return shouldMoveGroup(e)
    ? sim.N.filter(x => groupMode === 'component'
        ? ((x.cpLive !== undefined ? x.cpLive : x.cp) === (n.cpLive !== undefined ? n.cpLive : n.cp))
        : x.c === n.c)
    : [n];
}

canvas.onmousedown = e => {
  const n = nodeAt(e.offsetX, e.offsetY);
  if (n) {
    dragNode = n;
    dragGroup = dragTargets(n, e);
    dragGroup.forEach(x => x.fixed = true);
    reheat(0.3);
  } else {
    panStart = {x: e.offsetX, y: e.offsetY, ox: view.ox, oy: view.oy};
  }
};

window.onmouseup = () => {
  if (dragGroup) dragGroup.forEach(x => x.fixed = false);
  dragNode = null; dragGroup = null; panStart = null;
};

canvas.onmousemove = e => {
  if (dragNode) {
    const p = toWorld(e.offsetX, e.offsetY);
    const dx = p.x - dragNode.x, dy = p.y - dragNode.y;
    dragGroup.forEach(x => { x.x += dx; x.y += dy; });
    reheat(0.3);
    draw();
    return;
  }
  if (panStart) {
    view.ox = panStart.ox + (e.offsetX - panStart.x);
    view.oy = panStart.oy + (e.offsetY - panStart.y);
    draw();
    return;
  }
  const n = nodeAt(e.offsetX, e.offsetY);
  if (n) {
    canvas.style.cursor = 'pointer';
    tooltip.style.display = 'block';
    tooltip.style.left = (e.offsetX + 12) + 'px';
    tooltip.style.top = (e.offsetY - 8) + 'px';
    tooltip.textContent = n.id + (n.q ? '  (query)' : '') +
      '  ·  Louvain ' + n.c + '  ·  comp. ' +
      (n.cpLive !== undefined ? n.cpLive : n.cp) + '  ·  grau ' + n.deg;
  } else {
    canvas.style.cursor = 'grab';
    tooltip.style.display = 'none';
  }
};

canvas.onmouseleave = () => { tooltip.style.display = 'none'; };

canvas.onwheel = e => {
  e.preventDefault();
  const f = e.deltaY < 0 ? 1.12 : 0.89;
  const p = toWorld(e.offsetX, e.offsetY);
  view.scale = Math.max(0.15, Math.min(15, view.scale * f));
  view.ox = e.offsetX - p.x * view.scale;
  view.oy = e.offsetY - p.y * view.scale;
  draw();
};

// --- toque (iPhone/Android): arrastar no com 1 dedo, deslocar o fundo com
// 1 dedo fora de um no, pinça de 2 dedos pra zoom. Nao existe Shift no
// toque, entao o botao "mover grupo" faz esse papel nas duas plataformas.
function touchPos(t) {
  const rect = canvas.getBoundingClientRect();
  return {x: t.clientX - rect.left, y: t.clientY - rect.top};
}

let pinchDist0 = null, pinchScale0 = null;

canvas.addEventListener('touchstart', e => {
  if (e.touches.length === 2) {
    const [a, b] = e.touches;
    pinchDist0 = Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY);
    pinchScale0 = view.scale;
    dragNode = null; dragGroup = null; panStart = null;
    return;
  }
  if (e.touches.length === 1) {
    const p = touchPos(e.touches[0]);
    const n = nodeAt(p.x, p.y);
    if (n) {
      dragNode = n;
      dragGroup = dragTargets(n);
      dragGroup.forEach(x => x.fixed = true);
      reheat(0.3);
    } else {
      panStart = {x: p.x, y: p.y, ox: view.ox, oy: view.oy};
    }
  }
}, {passive: true});

canvas.addEventListener('touchmove', e => {
  if (e.touches.length === 2 && pinchDist0) {
    e.preventDefault();
    const [a, b] = e.touches;
    const d = Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY);
    const rect = canvas.getBoundingClientRect();
    const midX = (a.clientX + b.clientX) / 2 - rect.left;
    const midY = (a.clientY + b.clientY) / 2 - rect.top;
    const p = toWorld(midX, midY);
    view.scale = Math.max(0.15, Math.min(15, pinchScale0 * (d / pinchDist0)));
    view.ox = midX - p.x * view.scale;
    view.oy = midY - p.y * view.scale;
    draw();
    return;
  }
  if (e.touches.length === 1) {
    e.preventDefault();  // trava o scroll da pagina enquanto mexe na rede
    const p = touchPos(e.touches[0]);
    if (dragNode) {
      const wp = toWorld(p.x, p.y);
      const dx = wp.x - dragNode.x, dy = wp.y - dragNode.y;
      dragGroup.forEach(x => { x.x += dx; x.y += dy; });
      reheat(0.3);
      draw();
      return;
    }
    if (panStart) {
      view.ox = panStart.ox + (p.x - panStart.x);
      view.oy = panStart.oy + (p.y - panStart.y);
      draw();
    }
  }
}, {passive: false});

canvas.addEventListener('touchend', e => {
  if (e.touches.length < 2) pinchDist0 = null;
  if (e.touches.length === 0) {
    if (dragGroup) dragGroup.forEach(x => x.fixed = false);
    dragNode = null; dragGroup = null; panStart = null;
  }
}, {passive: true});

window.onresize = () => { if (cur && !cur.error) startSim(); };
// ---------- curva do closeness_scan (aba Rede, embaixo da visualizacao) ----------
// ---------- aba Alinhamento ----------
// Esquema de cor por propriedade fisico-quimica (convencao tipo ClustalX,
// ajustada pra fundo escuro). Nao e o mesmo codigo do rdbs.sequence.view()
// -- nao tenho acesso ao fonte dele pra replicar exatamente -- mas segue
// o mesmo espirito de alinhamento colorido calibrado por residuo.
const AA_COLOR = {
  A:'#5b8fd4', V:'#5b8fd4', L:'#5b8fd4', I:'#5b8fd4', M:'#5b8fd4',
  F:'#7a6bd6', W:'#7a6bd6', Y:'#5bb8c4',
  K:'#e05f5f', R:'#e05f5f',
  D:'#c766d1', E:'#c766d1',
  N:'#4fbf7a', Q:'#4fbf7a', S:'#4fbf7a', T:'#4fbf7a',
  C:'#e0a15f', G:'#d4b45b', P:'#d4b45b',
  H:'#5bc4bb',
};

function renderMsa() {
  const el = document.getElementById('pane-msa');
  if (!cur || !cur.msa) { el.innerHTML = ''; return; }

  const recs = cur.msa.records;  // [[header, seq], ...]
  const alnLen = recs.length ? recs[0][1].length : 0;
  const sameLen = recs.every(([, s]) => s.length === alnLen);

  let rows = '';
  for (const [h, seq] of recs) {
    let spans = '';
    for (const c of seq) {
      const up = c.toUpperCase();
      const bg = AA_COLOR[up];
      spans += bg
        ? '<span style="background:' + bg + '22;color:' + bg + '">' + c + '</span>'
        : '<span style="color:var(--faint)">' + c + '</span>';
    }
    const isQuery = cur.query_nodes && cur.query_nodes.includes(h);
    rows += '<div class="msa-row' + (isQuery ? ' msa-query' : '') + '" data-id="' +
      h.replace(/"/g, '&quot;') + '">' +
      '<div class="msa-label">' + h + (isQuery ? ' (query)' : '') + '</div>' +
      '<div class="msa-seq">' + spans + '</div></div>';
  }

  const lenWarn = sameLen ? '' :
    '<div class="warn">As sequencias deste arquivo nao tem todas o mesmo comprimento -- ' +
    'pode nao ser um alinhamento de verdade, so um fasta comum.</div>';

  el.innerHTML =
    '<div class="cards">' +
    card('Arquivo', cur.msa.file) +
    card('Sequências', recs.length) +
    card('Comprimento', sameLen ? alnLen : '—') +
    '</div>' +
    lenWarn +
    '<p class="hint">Linhas esmaecidas = a sequência não está mais no grupo da query ' +
    'no corte atual do slider (aba Rede). Isso não realinha nada -- é só um filtro visual ' +
    'sobre o alinhamento já pronto.</p>' +
    '<div class="msa-wrap">' + rows + '</div>';

  updateMsaDimming();
}

// Chamado a cada mudanca de corte (barato -- so toggle de classe, sem
// reconstruir o HTML colorido inteiro). E a parte "pelo menos o
// necessario" do pedido de realinhar ao vivo: nao roda um alinhador de
// verdade (isso precisaria de um processo rodando, nao da num HTML
// estatico sem servidor) -- mas mostra, sem atraso, quem ainda pertence
// ao grupo da query no corte escolhido, e quem nao pertence mais.
function updateMsaDimming() {
  if (!cur || !cur.msa || !sim) return;
  const byId = new Map(sim.N.map(n => [n.id, n]));
  document.querySelectorAll('#pane-msa .msa-row').forEach(row => {
    const n = byId.get(row.dataset.id);
    const inGroup = n && (n.qcpLive !== undefined ? n.qcpLive : n.qcp);
    row.classList.toggle('out-of-group', n ? !inGroup : false);
  });
}

function renderScanSection() {
  const el = document.getElementById('scan-section');
  if (!cur || cur.error) { el.innerHTML = ''; return; }

  let scanWarn = '';
  if (cur.scan_truncated) {
    scanWarn = '<div class="warn">A curva ainda estava subindo quando parou' +
      (cur.scan_max ? ' (máximo alcançado: ' + cur.scan_max.toFixed(2) + ')' : '') +
      '. Isso normalmente <b>não</b> é falta de dado acima desse ponto — é a regra de ' +
      'parada por fragmentação (<code>max_singleton_fraction</code>, padrão 50%) ' +
      'interrompendo a varredura cedo. Aumentar <code>--cutoff-steps</code> sozinho não ' +
      'resolve, porque a parada é por fragmentação, não por número de passos. Se o slider ' +
      'abaixo alcança valores bem acima de onde esta curva parou, o pico real pode estar ' +
      'lá — vale testar direto nele antes de reprocessar.</div>';
  } else if (cur.scan_max !== null && cur.scan_max !== undefined && cur.scan_max < 0.7) {
    scanWarn = '<div class="warn">Closeness máximo de ' + cur.scan_max.toFixed(2) +
      ' — abaixo do que caracteriza uma separação limpa (referência: ~0,9–1,0). ' +
      'Os grupos podem não estar bem separados neste corte.</div>';
  }

  if (!cur.scan_png) {
    el.innerHTML = '<h3>Curva do closeness_scan</h3>' + scanWarn +
      '<p class="empty">Sem imagem da curva para este domínio ' +
      '(rode o estágio ssn para gerá-la).</p>';
    return;
  }

  el.innerHTML = '<h3>Curva do closeness_scan (escolha do corte)</h3>' + scanWarn +
    '<img class="scan-img" src="data:image/png;base64,' + cur.scan_png + '">';
}

// ---------- aba matriz ----------
function renderMatrix() {
  const el = document.getElementById('pane-matrix');
  const m = cur.matrix;
  if (!m || m.error) {
    el.innerHTML = '<p class="empty">Sem dados de matriz para este domínio' +
      (m && m.error ? ' (' + m.error + ')' : '') + '.</p>';
    return;
  }
  const maxC = Math.max(...m.hist_counts);
  const cutoff = cur.auto_cutoff;

  // histograma vertical (SVG) -- compacto, uma barra por faixa, com a
  // faixa do corte destacada em magenta. Substitui a antiga tabela de
  // barras horizontais, que ocupava uma linha inteira por faixa.
  const hW = 720, hH = 260, hPad = 4, barGap = 1.5;
  const nBins = m.hist_counts.length;
  const barW = (hW / nBins) - barGap;
  let cutoffBinIdx = -1;
  const svgBars = m.hist_counts.map((c, i) => {
    const lo = m.hist_edges[i], hi = m.hist_edges[i+1];
    const isCutoffBin = (cutoff !== null && cutoff !== undefined && cutoff >= lo && cutoff < hi);
    if (isCutoffBin) cutoffBinIdx = i;
    const barH = maxC ? (c / maxC) * (hH - hPad) : 0;
    const x = i * (hW / nBins);
    const y = hH - barH;
    const fill = isCutoffBin ? 'url(#gradCutoff)' : 'url(#gradBar)';
    return '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + barW.toFixed(1) +
      '" height="' + barH.toFixed(1) + '" fill="' + fill + '">' +
      '<title>' + lo + '–' + hi + ': ' + c.toLocaleString('pt-BR') + ' par(es)' +
      (isCutoffBin ? ' ◂ corte' : '') + '</title></rect>';
  }).join('');

  // marcador de linha vertical no centro exato do corte (mais preciso que
  // so destacar o bin inteiro, quando o corte cai no meio de uma faixa)
  let cutoffLine = '';
  if (cutoff !== null && cutoff !== undefined && m.hist_edges.length > 1) {
    const loEdge = m.hist_edges[0], hiEdge = m.hist_edges[m.hist_edges.length - 1];
    if (cutoff >= loEdge && cutoff <= hiEdge) {
      const xPos = ((cutoff - loEdge) / (hiEdge - loEdge)) * hW;
      cutoffLine = '<line x1="' + xPos.toFixed(1) + '" y1="0" x2="' + xPos.toFixed(1) +
        '" y2="' + hH + '" stroke="#ff5f9e" stroke-width="1.5" stroke-dasharray="3,2"/>';
    }
  }

  // grade horizontal com valores -- ajuda a ler a contagem real, nao so
  // comparar alturas de olho
  const nGrid = 4;
  let gridLines = '';
  for (let g = 1; g <= nGrid; g++) {
    const gy = hH - (g / nGrid) * (hH - hPad);
    const gVal = Math.round((g / nGrid) * maxC);
    gridLines +=
      '<line x1="0" y1="' + gy.toFixed(1) + '" x2="' + hW + '" y2="' + gy.toFixed(1) +
        '" stroke="#1f2c42" stroke-width="1"/>' +
      '<text x="4" y="' + (gy - 3).toFixed(1) + '" font-family="ui-monospace,monospace" ' +
        'font-size="9" fill="#445674">' + gVal.toLocaleString('pt-BR') + '</text>';
  }

  // mais marcacoes no eixo x (5 pontos, nao so as duas pontas)
  const nXTicks = 5;
  let xTicks = '';
  for (let t = 0; t <= nXTicks; t++) {
    const tx = (t / nXTicks) * hW;
    const tVal = m.hist_edges[0] + (t / nXTicks) * (m.hist_edges[m.hist_edges.length-1] - m.hist_edges[0]);
    const anchor = t === 0 ? 'start' : (t === nXTicks ? 'end' : 'middle');
    xTicks += '<text x="' + tx.toFixed(1) + '" y="' + (hH + 16) +
      '" font-family="ui-monospace,monospace" font-size="10" fill="#7488a8" ' +
      'text-anchor="' + anchor + '">' + tVal.toFixed(0) + '</text>';
  }

  const histogram =
    '<svg viewBox="0 0 ' + hW + ' ' + (hH + 22) + '" style="width:100%;height:auto;display:block">' +
      '<defs>' +
        '<linearGradient id="gradBar" x1="0" y1="0" x2="0" y2="1">' +
          '<stop offset="0%" stop-color="#4fe8d6"/><stop offset="100%" stop-color="#1f6b62"/>' +
        '</linearGradient>' +
        '<linearGradient id="gradCutoff" x1="0" y1="0" x2="0" y2="1">' +
          '<stop offset="0%" stop-color="#ff8fbd"/><stop offset="100%" stop-color="#b23d6d"/>' +
        '</linearGradient>' +
      '</defs>' +
      gridLines + svgBars + cutoffLine + xTicks +
    '</svg>';

  const histCaption =
    '<p class="hint">Cada barra conta PARES (arestas candidatas) da matriz bruta nessa faixa ' +
    'de bitscore — não indica quantos nós ficam conectados nem se a rede se separa bem. ' +
    'Para decidir o corte, use a curva de closeness na aba <b>Rede (SSN)</b>.</p>';

  const cutoffCard = cutoff !== null && cutoff !== undefined
    ? card('Corte escolhido (auto_cutoff)', cutoff.toFixed(1), 'warn')
    : card('Corte escolhido', '—');

  // segundo grafico: distribuicao de GRAU dos nos na rede ja filtrada
  // (pelo auto_cutoff). Eixo diferente do histograma de pares -- em vez
  // de "quantos pares tem esse score", responde "quantas conexoes cada
  // no tem", que e informacao sobre a TOPOLOGIA da rede resultante, nao
  // sobre a matriz bruta.
  const degrees = cur.nodes.map(n => n.deg).filter(d => d !== undefined);
  let degreeSection = '';
  if (degrees.length) {
    const dMax = Math.max(...degrees);
    const dBins = Math.min(20, dMax + 1);
    const dBinSize = Math.max(1, Math.ceil((dMax + 1) / dBins));
    const dCounts = new Array(Math.ceil((dMax + 1) / dBinSize)).fill(0);
    degrees.forEach(d => dCounts[Math.min(dCounts.length - 1, Math.floor(d / dBinSize))]++);
    const dMaxC = Math.max(...dCounts);

    const dW = 720, dH = 150;
    const dBarW = (dW / dCounts.length) - 1.5;
    const dBars = dCounts.map((c, i) => {
      const lo = i * dBinSize, hi = lo + dBinSize - 1;
      const barH = dMaxC ? (c / dMaxC) * (dH - 4) : 0;
      const x = i * (dW / dCounts.length);
      return '<rect x="' + x.toFixed(1) + '" y="' + (dH - barH).toFixed(1) + '" width="' +
        dBarW.toFixed(1) + '" height="' + barH.toFixed(1) + '" fill="url(#gradDeg)">' +
        '<title>grau ' + lo + (hi > lo ? '–' + hi : '') + ': ' + c + ' nó(s)</title></rect>';
    }).join('');

    const isolados = degrees.filter(d => d === 0).length;
    const medDeg = degrees.slice().sort((a,b)=>a-b)[Math.floor(degrees.length/2)];

    degreeSection =
      '<h3>Distribuição de grau dos nós (rede após o auto_cutoff)</h3>' +
      '<div class="cards">' +
      card('Grau mediano', medDeg) +
      card('Grau máximo', dMax) +
      card('Nós isolados', isolados) +
      '</div>' +
      '<div class="hist-wrap"><svg viewBox="0 0 ' + dW + ' ' + (dH+4) + '" ' +
        'style="width:100%;height:auto;display:block">' +
        '<defs><linearGradient id="gradDeg" x1="0" y1="0" x2="0" y2="1">' +
        '<stop offset="0%" stop-color="#33d9c4"/><stop offset="100%" stop-color="#175b52"/>' +
        '</linearGradient></defs>' + dBars + '</svg></div>' +
      '<p class="hint">Quantas conexões (arestas) cada nó tem, na rede já filtrada pelo corte ' +
      'automático — nós isolados (grau 0) sobrevivem no grafo mas sem nenhuma aresta acima do corte.</p>';
  }

  el.innerHTML =
    '<div class="cards">' +
    card('Pares na matriz', m.n_pairs.toLocaleString('pt-BR')) +
    card('Bitscore mín.', m.min.toFixed(1)) +
    card('Bitscore máx.', m.max.toFixed(1)) +
    card('Mediana', m.median.toFixed(1)) +
    cutoffCard +
    '</div>' +
    '<h3>Distribuição de bitscore (matriz bruta, antes do corte)</h3>' +
    '<div class="hist-wrap">' + histogram + '</div>' +
    histCaption +
    degreeSection;
}

// ---------- aba comunidades ----------
function renderStats() {
  const el = document.getElementById('pane-stats');
  if (!cur.comm_sizes || !cur.comm_sizes.length) {
    el.innerHTML = '<p class="empty">Sem informação de comunidade neste grafo.</p>';
    return;
  }
  const total = cur.n_nodes_total;
  const maxS = cur.comm_sizes[0][1];
  let rows = '';
  cur.comm_sizes.forEach(([cid, size]) => {
    const isQ = cur.query_comm !== null && cid === cur.query_comm;
    rows += '<tr class="' + (isQ ? 'is-query' : '') + '">' +
      '<td>' + cid + (isQ ? '  ← query' : '') + '</td>' +
      '<td style="width:60%"><div class="bar' + (isQ ? ' q' : '') +
        '" style="width:' + (size/maxS*100) + '%"></div></td>' +
      '<td class="num">' + size.toLocaleString('pt-BR') + '</td>' +
      '<td class="num">' + (size/total*100).toFixed(1) + '%</td></tr>';
  });
  el.innerHTML =
    '<div class="cards">' +
    card('Nós no grafo', total.toLocaleString('pt-BR')) +
    card('Comunidades', cur.n_communities) +
    card('Comunidade da query', cur.query_comm !== null ? cur.query_comm : '—', true) +
    card('Maior comunidade', maxS.toLocaleString('pt-BR') + ' nós') +
    '</div>' +
    '<h3>Tamanho por comunidade (15 maiores)</h3>' +
    '<table><thead><tr><th>comunidade</th><th></th><th class="num">nós</th>' +
    '<th class="num">%</th></tr></thead><tbody>' + rows + '</tbody></table>';
}

function card(label, value, variant) {
  const cls = variant === 'warn' ? ' warn' : (variant ? ' accent' : '');
  return '<div class="card"><div class="label">' + label + '</div>' +
         '<div class="value' + cls + '">' + value + '</div></div>';
}

document.getElementById('group-mode-btn').onclick = () => {
  moveGroupMode = !moveGroupMode;
  document.getElementById('group-mode-btn').textContent =
    'modo: ' + (moveGroupMode ? 'grupo inteiro' : 'nó único');
  document.getElementById('group-mode-btn').style.color =
    moveGroupMode ? 'var(--magenta)' : 'var(--cyan)';
};

document.getElementById('mode-btn').onclick = () => {
  groupMode = groupMode === 'louvain' ? 'component' : 'louvain';
  document.getElementById('mode-btn').textContent =
    'agrupar por: ' + (groupMode === 'component' ? 'componente conectado' : 'Louvain (ao vivo)');
  updateLegendCounts();
  draw();
};

function updateLegendCounts() {
  const el = document.getElementById('legend-counts');
  if (!cur || cur.error) { el.textContent = ''; return; }
  const liveSize = (sim && sim.qCompSize !== undefined) ? sim.qCompSize : cur.query_component_size;
  const liveN = (sim && sim.nComponents !== undefined) ? sim.nComponents : cur.n_components;
  el.textContent = groupMode === 'component'
    ? (liveSize || 0) + ' nós no componente da query (de ' + (liveN || '?') + ' componentes)'
    : (cur.query_comm_size || 0) + ' nós na comunidade da query (de ' +
      (cur.n_communities || '?') + ' comunidades)';
}

// inicia no primeiro dominio valido
// restaura marcacoes de corte manual salvas de uma sessao anterior
Object.keys(manualCutoffs).forEach(markDomainOverridden);
updateExportBtn();

const first = DATA.findIndex(d => !d.error);
if (first >= 0) select(first);
</script>
</body>
</html>
"""




# =============================================================================
# dashboard -- painel HTML interativo, autocontido
# =============================================================================


def apply_manual_cutoffs(overrides, matrix_dir, ssn_dir, out_dir=None,
                         cutoff_column='bitscore', add_community=True):
    """
    Apply cutoffs chosen by hand in the dashboard slider back onto the
    pipeline's SSN outputs.

    The panel is a static file and cannot write to disk, so it exports the
    chosen cutoffs as JSON; this rebuilds each named network at that
    threshold instead of the one auto_cutoff picked.

    Parameters
    ----------
    overrides : dict or str
        {unit: cutoff_value}, or a path to the JSON exported by the
        dashboard.
    matrix_dir : str
        Directory of `<unit>_hits.tsv` all-vs-all tables (same one the
        `ssn` stage read from).
    ssn_dir : str
        Directory holding the existing `<unit>.graphml` files to update.
    out_dir : str or None
        Where to write the rebuilt files. None overwrites in place.
    cutoff_column : str, default 'bitscore'
        Must match whatever column the original run scored on.
    add_community : bool, default True
        Recompute Louvain communities for the rebuilt network too.

    Returns
    -------
    dict
        {unit: {'nodes': int, 'edges': int, 'components': int}}, one entry
        per unit successfully rebuilt. Units that failed are logged as
        warnings and omitted, not raised, so one bad entry in a batch
        doesn't stop the rest.
    """
    if isinstance(overrides, str):
        with open(overrides) as fh:
            overrides = json.load(fh)

    out_dir = out_dir or ssn_dir
    os.makedirs(out_dir, exist_ok=True)
    results = {}

    for unit, cutoff in overrides.items():
        matrix_path = os.path.join(matrix_dir, f"{unit}_hits.tsv")
        if not os.path.exists(matrix_path):
            logger.warning(f"apply_manual_cutoffs: {unit}: no matrix at {matrix_path}, skipping")
            continue
        try:
            hits = pd.read_csv(matrix_path, sep="\t", low_memory=False)
            G = build_ssn(hits, criteria=f"{cutoff_column} > {cutoff}",
                          add_community=add_community, progress=False)
            G.graph["auto_cutoff"] = f"{cutoff_column} > {cutoff}"
            G.graph["cutoff_source"] = "manual"  # distingue de escolha automatica

            # keep the scan curve: it shows where a manual cutoff sits
            # relative to the peaks
            old_gml = os.path.join(ssn_dir, f"{unit}.graphml")
            if os.path.exists(old_gml):
                try:
                    old_G = nx.read_graphml(old_gml)
                    if "closeness_scan" in old_G.graph:
                        G.graph["closeness_scan"] = old_G.graph["closeness_scan"]
                except Exception:
                    pass

            out_path = os.path.join(out_dir, f"{unit}.graphml")
            export_ssn(G, out_path)

            n_comp = nx.number_connected_components(G) if G.number_of_nodes() else 0
            results[unit] = {"nodes": G.number_of_nodes(),
                             "edges": G.number_of_edges(), "components": n_comp}
            logger.info(f"apply_manual_cutoffs: {unit}: cutoff={cutoff} -> "
                       f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
                       f"{n_comp} components")
        except Exception as e:
            logger.warning(f"apply_manual_cutoffs: {unit}: failed ({e})")

    return results


def ssn_dashboard(ssn_dir, out, matrix_dir=None, aliases=None,
                  max_nodes=None, layout_iter=50, progress=True,
                  cutoff_overrides=None, msa_dir=None):
    """
    Build one self-contained HTML panel from a directory of SSN graphs.

    Data, layout and code all travel inside the file, so it opens by
    double-click with no server and no network access.

    Per network: the graph under a live force simulation (drag a node, or
    Shift-drag its whole group), a cutoff slider that recomputes connected
    components on the fly, the bitscore distribution with the cutoff
    marked, and the community table. The closeness scan plot is embedded
    when a PNG of the same basename sits next to the .graphml; `msa_dir`
    adds an alignment tab.

    Grouping defaults to CONNECTED COMPONENT rather than Louvain: that is
    what the closeness-centrality criterion optimises (Hornung & Terrapon
    2023), since the scan looks for the cutoff at which pieces come apart.
    Modularity partitions a still-connected graph anyway. Both views are
    available; the button switches live.

    The layout is computed in the browser rather than baked in, which is
    what makes nodes draggable. Repulsion uses a Barnes-Hut quadtree, so
    a few thousand nodes stay interactive.

    Parameters
    ----------
    ssn_dir : str
        Directory of `*.graphml` files, one per unit. Optional
        `<same-basename>.png` files are embedded as the scan plot.
    out : str
        Path of the HTML file to write.
    matrix_dir : str or None
        Directory of `<unit>_hits.tsv` all-vs-all tables, used for the
        bitscore distribution panel and for the live cutoff slider. Also
        what `cutoff_overrides` is evaluated against -- without it, an
        override is stored but the slider still starts from the graph's
        own baked-in cutoff, since there is no raw table to recompute
        components from.
    aliases : dict or None
        Optional {unit: [id, ...]} giving alternative accessions for each
        unit's reference sequence, for when it is stored in the network
        under a different id than the unit name.
    max_nodes : int or None, default None
        With None the whole network is drawn. Set an integer only if a
        very large network makes the simulation sluggish; the reference
        node's own group is always kept whole.
    layout_iter : int, default 50
        Iterations for the seed layout that the browser starts from.
    cutoff_overrides : dict or str, or None
        {unit: bitscore}, or a path to the JSON the dashboard's "baixar
        cortes manuais" button downloads. Where a unit has an override,
        the panel opens with the slider already at that value instead of
        the graph's own auto_cutoff -- no need to have rebuilt the
        `.graphml` first (`apply_manual_cutoffs`) just to preview this.
        Rebuilding is still what makes the choice permanent in the
        pipeline's own outputs (community assignment, `select`, ...);
        this only changes what the panel shows on open.
    msa_dir : str or None
        Directory to scan for `*.msa` (or `*.aln`/`*.fasta` with aligned,
        equal-length sequences) files. A file is attached to unit X's
        panel when X's reference sequence is one of its records -- files
        without it are skipped, since the point of the tab is to inspect
        the query's own alignment, not every alignment that happens to
        exist. Rendered as a coloured alignment view.

    Returns
    -------
    str
        The path written.
    """
    import glob as _glob
    graphmls = sorted(_glob.glob(os.path.join(ssn_dir, "*.graphml")))
    if not graphmls:
        raise FileNotFoundError(f"no .graphml files in {ssn_dir}")

    if isinstance(cutoff_overrides, str):
        with open(cutoff_overrides) as fh:
            cutoff_overrides = json.load(fh)
    cutoff_overrides = cutoff_overrides or {}

    msa_files = _dash_scan_msa_files(msa_dir) if msa_dir else []

    data = []
    for gml in graphmls:
        name = os.path.splitext(os.path.basename(gml))[0]
        mpath = os.path.join(matrix_dir, f"{name}_hits.tsv") if matrix_dir else None
        d = _dash_build_unit(gml, mpath, aliases or {}, max_nodes, layout_iter)

        if name in cutoff_overrides and not d.get("error"):
            d["auto_cutoff"] = float(cutoff_overrides[name])
            d["cutoff_is_manual"] = True

        if msa_files and not d.get("error") and d.get("query_nodes"):
            qset = set(d["query_nodes"])
            for fp, records in msa_files:
                headers = {h for h, _ in records}
                if qset & headers:
                    d["msa"] = {"file": os.path.basename(fp), "records": records}
                    break

        if progress:
            if d.get("error"):
                logger.warning(f"ssn_dashboard: {name}: {d['error']}")
            else:
                logger.info(f"ssn_dashboard: {name}: {d['n_nodes_total']} nodes, "
                            f"{d.get('n_components', '?')} components")
        data.append(d)

    html = (_DASH_TEMPLATE
            .replace("__DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__NDOM__", str(len(data))))
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info(f"ssn_dashboard: wrote {out} ({os.path.getsize(out)/1e6:.1f} MB)")
    return out

ENGINES = ('iqtree2', 'fasttree')


def _tree_command(engine, binary, aln, prefix, model, bb, alrt, moltype, gamma, cpu, seed, redo, extra_args):
    """
    Build the shell command for one tree-building engine. Both engines
    are made to write the final Newick tree to `{prefix}.treefile`, so
    the rest of build_tree never needs an engine-specific code path
    for reading the result back.
    """
    treefile = f'{prefix}.treefile'

    if engine == 'iqtree2':
        cmd = f'{binary} -s {aln} -m {model} -nt {cpu} -pre {prefix}'
        if bb:
            cmd += f' -bb {bb}'
        if alrt:
            cmd += f' -alrt {alrt}'
        if seed is not None:
            cmd += f' -seed {seed}'
        if redo:
            cmd += ' -redo'
        if extra_args:
            cmd += f' {extra_args}'

    elif engine == 'fasttree':
        flags = '-nt -gtr' if moltype == 'nucl' else ''
        if gamma:
            flags += ' -gamma'
        if extra_args:
            flags += f' {extra_args}'
        # FastTree has no -threads flag -- multithreading (with the MP
        # build) is controlled by the OMP_NUM_THREADS env var instead.
        env_prefix = f'OMP_NUM_THREADS={cpu} ' if cpu else ''
        cmd = f'{env_prefix}{binary} {flags} {aln} > {treefile} 2> {prefix}.log'

    else:
        raise ValueError(f"engine must be one of {ENGINES}, got {engine!r}")

    return cmd, treefile


def build_tree(seqobj, engine='fasttree', binary=None, model='MFP', bb=1000, alrt=1000,
               moltype='prot', gamma=True, cpu=8, seed=None, redo=True, extra_args='',
               outdir=None, path=None, progress=True):
    """
    Build a tree, called directly as a subprocess -- same pattern as
    all_vs_all_search (tempfile.TemporaryDirectory + Popen). Just runs
    the command with whatever parameters you set and returns the
    Newick text; no tree-parsing library involved at all.

    Parameters
    ----------
    seqobj : rotifer.devel.beta.sequence.sequence
        Must already be an ALIGNMENT (gaps kept, all sequences the same
        length) -- call `.align()` first if it isn't one yet.
    engine : str, default 'fasttree'
        'fasttree' or 'iqtree2'. These are NOT statistically
        equivalent, read on before treating their numbers as
        comparable:

        - 'fasttree' (default): approximate ML, much faster (built for
          alignments too large for iqtree2 to be practical). Its only
          support metric is a single Shimodaira-Hasegawa-like LOCAL
          support test, computed automatically -- there is no
          replicate count to set, so `bb`/`alrt` are simply ignored
          for this engine.
        - 'iqtree2': full maximum likelihood, ModelFinder (-m MFP by
          default), ultrafast bootstrap (-bb) and SH-aLRT (-alrt) as
          independently tunable replicate counts. `bb`/`alrt` already
          default to 1000/1000 -- just switch engine='iqtree2' and
          they're ready to go, no need to set them yourself. Slower,
          statistically the more rigorous of the two.

          Do not read a FastTree support value and an IQ-TREE UFBoot
          value as if they meant the same thing.
    binary : str or None, default None
        Executable name/path. Defaults to 'iqtree2' or 'fasttree'
        depending on `engine`. Use 'iqtree' or 'FastTreeMP' (the
        OpenMP-threaded build) if that's what's on your PATH instead.
    model : str, default 'MFP'
        iqtree2 only (-m). 'MFP' (ModelFinder Plus) selects the best
        substitution model automatically as part of the same run.
    bb : int or None, default 1000
        iqtree2 only. Ultrafast bootstrap replicates (-bb). IQ-TREE
        requires >= 1000 if set. None disables UFBoot.
    alrt : int or None, default 1000
        iqtree2 only. SH-aLRT replicates (-alrt). None disables it.
    moltype : str, default 'prot'
        fasttree only. 'prot' (default, no flag needed) or 'nucl'
        (adds -nt -gtr). iqtree2 infers this from the alignment itself.
    gamma : bool, default True
        fasttree only. Adds -gamma (rescales branch lengths under a
        discrete gamma model) -- FastTree's own recommendation for
        branch lengths comparable to other ML tools.
    cpu : int, default 8
        Threads. Mapped to -nt for iqtree2, to the OMP_NUM_THREADS
        environment variable for fasttree (only takes effect with an
        OpenMP-enabled build, e.g. binary='FastTreeMP').
    seed : int or None, default None
        iqtree2 only (-seed), for reproducible runs.
    redo : bool, default True
        iqtree2 only. Passes -redo, so re-running in the same
        directory doesn't error out on pre-existing output files.
    extra_args : str, default ''
        Any extra flags for the underlying engine, appended verbatim.
    outdir : str or None, default None
        If given, copies ALL of the engine's output there (for
        iqtree2: .treefile, .iqtree report, .log, .contree, .bionj...;
        for fasttree: .treefile and .log) -- otherwise everything but
        what `path` points to is lost when the temp dir is cleaned up.
    path : str or None, default None
        If given, copies just the Newick tree there.
    progress : bool, default True
        Show a tqdm progress bar with the current pipeline stage.

    Returns
    -------
    str
        The raw Newick string. Open it in TreeViewer, feed it to ete3
        yourself (`ete3.Tree(newick, format=1)`), or whatever else you
        use -- this function doesn't assume which one.

    Examples
    --------
    >>> aln = rdbs.sequence('alinhamento.msa').align()
    >>> newick_fast = rdao.build_tree(aln, cpu=8, outdir='resultados_fasttree')
    >>> newick = rdao.build_tree(aln, engine='iqtree2', cpu=8, outdir='resultados_iqtree')
    """
    if engine not in ENGINES:
        raise ValueError(f"engine must be one of {ENGINES}, got {engine!r}")
    if binary is None:
        binary = engine

    pbar = tqdm(total=3, desc='build_tree', unit='step', disable=not progress)

    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)

        pbar.set_description('writing alignment')
        aln_path = f'{tmp}/aln.fasta'
        seqobj.to_file(aln_path)
        pbar.update(1)

        pbar.set_description(f'running {engine}')
        prefix = f'{tmp}/tree_run'
        cmd, treefile = _tree_command(engine, binary, aln_path, prefix, model, bb, alrt,
                                       moltype, gamma, cpu, seed, redo, extra_args)
        proc = Popen(cmd, stdout=PIPE, stderr=STDOUT, shell=True)
        log, _ = proc.communicate()
        pbar.update(1)

        pbar.set_description('reading results')
        if not os.path.exists(treefile) or not os.path.getsize(treefile):
            os.chdir(cwd)
            pbar.close()
            logger.error(f'{binary} did not produce a tree -- log follows')
            raise RuntimeError(log.decode(errors='replace'))

        with open(treefile) as fh:
            newick = fh.read().strip()
        pbar.update(1)

        if outdir:
            import shutil
            os.makedirs(outdir, exist_ok=True)
            for f in os.listdir(tmp):
                if f.startswith('tree_run'):
                    shutil.copy(f'{tmp}/{f}', outdir)
            logger.info(f'all {engine} output copied to {outdir}')
        if path:
            import shutil
            shutil.copy(treefile, path)
            logger.info(f'tree saved to {path}')
    os.chdir(cwd)
    pbar.close()
    return newick


def split_support_labels(tree):
    """
    Optional, standalone helper -- entirely independent of build_tree.

    If you built the tree with both `bb` and `alrt` set, IQ-TREE writes
    a combined "SH-aLRT/UFBoot" label on each internal node (e.g.
    "95.5/98"). ete3 has no built-in way to split that into two
    values, so this does it for you, in place, on a tree YOU already
    loaded with ete3 -- this module never imports ete3 itself.

    Parameters
    ----------
    tree : ete3.Tree

    Returns
    -------
    ete3.Tree
        The same object, modified in place (also returned for
        convenience). Internal nodes gain `node.sh_alrt` and
        `node.ufboot` (floats); `node.name` is cleared on those nodes.

    Examples
    --------
    >>> import ete3
    >>> from rotifer.devel.alpha import mvroliveira as rdao
    >>> newick = rdao.build_tree(aln, cpu=8)
    >>> tree = ete3.Tree(newick, format=1)
    >>> rdao.split_support_labels(tree)
    >>> [n for n in tree.traverse() if not n.is_leaf() and n.ufboot >= 95]
    """
    for node in tree.traverse():
        if not node.is_leaf() and node.name and '/' in node.name:
            sh, uf = node.name.split('/')
            node.add_feature('sh_alrt', float(sh))
            node.add_feature('ufboot', float(uf))
            node.name = ''
    return tree


def to_tree(self, path=None, **kwargs):
    """
    Build a maximum-likelihood tree straight from this sequence object
    (IQ-TREE2 or FastTree, see build_tree's `engine`). Attached to
    every rotifer sequence object as `.to_tree()` -- just import this
    module once (same monkey-patch mechanism as `.to_ssn()`).

    `self` must already be an alignment (call `.align()` first if not).

    Parameters
    ----------
    path : str or None, default None
        If given, saves the tree (Newick) there.
    **kwargs
        Any argument accepted by build_tree (engine, binary, model, bb,
        alrt, moltype, gamma, cpu, seed, redo, extra_args, outdir,
        progress).

    Returns
    -------
    str
        The Newick tree string. Use split_support_labels() yourself,
        with your own ete3 import, if you want sh_alrt/ufboot split
        (see that function's docstring; iqtree2 only).

    Examples
    --------
    >>> seqs = rdbs.sequence('alinhamento.msa').align()
    >>> newick_fast = seqs.to_tree(cpu=8, outdir='resultados_fasttree')
    >>> newick = seqs.to_tree(engine='iqtree2', outdir='resultados_iqtree')
    """
    return build_tree(self, path=path, **kwargs)


def to_ssn(self, path=None, fmt='graphml', save_hits=True, **kwargs):
    """
    Run all-vs-all BLAST on this sequence object and build the ssn.

    Attached to every rotifer sequence object as `.to_ssn()` --
    just import this module once.

    Parameters
    ----------
    path : str or None, default None
        If given, saves the ssn there (see export_ssn).
        If None, nothing is written to disk -- only the graph is returned.
    fmt : str, default 'graphml'
        'graphml' or 'tsv', passed to export_ssn.
    save_hits : bool, default True
        If True and `path` is given, also saves the raw all-vs-all hit
        table right next to it -- same name, `_hits.tsv` suffix (e.g.
        path='LIC_11568_ssn.graphml' also writes
        'LIC_11568_ssn_hits.tsv'). Lets you rebuild the SSN later
        with a different `criteria`, without re-running the search:

            import pandas as pd
            hits = pd.read_csv('LIC_11568_ssn_hits.tsv', sep='\\t')
            G2 = rdao.build_ssn(hits, criteria='bitscore > 100')
    **kwargs
        Any argument accepted by sequence_ssn (criteria, cpu,
        keep_best, add_community, progress, ...).

    Returns
    -------
    networkx.Graph
        The raw hit table is also kept as `seq._ssn_hits` for the
        rest of the current session, regardless of `save_hits`/`path`.

    Examples
    --------
    >>> seqs = rdbs.sequence("teste.fasta")
    >>> G = seqs.to_ssn('rede.graphml')   # writes rede.graphml + rede_hits.tsv
    >>> rdao.to_cytoscape(G)
    """
    G, hits = sequence_ssn(self, **kwargs)
    self._ssn_hits = hits
    if path:
        export_ssn(G, path, fmt=fmt)
        if save_hits:
            hits_path = os.path.splitext(path)[0] + '_hits.tsv'
            hits.to_csv(hits_path, sep='\t', index=False)
            logger.info(f'raw hit table saved to {hits_path}')
    return G


def build_ssn_and_tree(seqobj, cpu=8, moltype='prot', progress=True,
                       ssn_kwargs=None, tree_kwargs=None):
    """
    Run sequence_ssn + build_tree in one call -- the similarity
    ssn AND the tree, without calling each function separately.

    `cpu`, `moltype` and `progress` are shared defaults applied to
    both steps. Anything specific to just one step (e.g. `criteria`
    for the ssn, `engine`/`bb`/`alrt` for the tree) goes in
    `ssn_kwargs`/`tree_kwargs` instead -- any argument
    sequence_ssn or build_tree accepts is valid there, and
    overrides the shared defaults if repeated.

    Parameters
    ----------
    seqobj : rotifer.devel.beta.sequence.sequence
        Must already be an alignment (call `.align()` first) -- both
        steps need one (all_vs_all_search strips gaps internally for
        its own purposes; build_tree needs them kept).
    cpu : int, default 8
    moltype : str, default 'prot'
    progress : bool, default True
    ssn_kwargs : dict or None, default None
        Extra/override arguments for sequence_ssn (tool, program,
        search_evalue, criteria, keep_best, add_community, ...).
    tree_kwargs : dict or None, default None
        Extra/override arguments for build_tree (engine, model, bb,
        alrt, gamma, outdir, path, ...).

    Returns
    -------
    dict with keys:
        'graph' : networkx.Graph
        'hits'  : pandas.DataFrame  (raw all-vs-all hit table)
        'tree'  : str               (Newick)

    Examples
    --------
    >>> aln = rdbs.sequence('alinhamento.msa').align()
    >>> result = rdao.build_ssn_and_tree(aln, cpu=8)
    >>> result['graph'], result['tree']

    >>> # overriding just one step:
    >>> result = rdao.build_ssn_and_tree(
    ...     aln, cpu=8,
    ...     ssn_kwargs={'tool': 'diamond', 'add_community': True},
    ...     tree_kwargs={'engine': 'iqtree2', 'outdir': 'resultados_iqtree'},
    ... )
    """
    ssn_kwargs = dict(ssn_kwargs or {})
    tree_kwargs = dict(tree_kwargs or {})
    ssn_kwargs.setdefault('cpu', cpu)
    ssn_kwargs.setdefault('moltype', moltype)
    ssn_kwargs.setdefault('progress', progress)
    tree_kwargs.setdefault('cpu', cpu)
    tree_kwargs.setdefault('moltype', moltype)
    tree_kwargs.setdefault('progress', progress)

    G, hits = sequence_ssn(seqobj, **ssn_kwargs)
    newick = build_tree(seqobj, **tree_kwargs)

    return {'graph': G, 'hits': hits, 'tree': newick}


def build_ssn_and_tree_method(self, **kwargs):
    """
    Run build_ssn_and_tree (ssn + tree) straight from this
    sequence object. Attached to every rotifer sequence object as
    `.build_ssn_and_tree()` -- just import this module once (same
    monkey-patch mechanism as `.to_ssn()`/`.to_tree()`).

    `self` must already be an alignment (call `.align()` first if not).

    Parameters
    ----------
    **kwargs
        Any argument accepted by build_ssn_and_tree (cpu, moltype,
        progress, ssn_kwargs, tree_kwargs).

    Returns
    -------
    dict with keys 'graph', 'hits', 'tree' -- see build_ssn_and_tree.

    Examples
    --------
    >>> seqs = rdbs.sequence('alinhamento.msa').align()
    >>> result = seqs.build_ssn_and_tree(cpu=8)
    >>> rdao.to_cytoscape(result['graph'])
    """
    return build_ssn_and_tree(self, **kwargs)


# Monkey-patch: attach to_ssn/to_tree/build_ssn_and_tree to the sequence class
# so any rdbs.sequence(...) object gets them just by importing this module.
rdbs.sequence.to_ssn = to_ssn
rdbs.sequence.to_tree = to_tree
rdbs.sequence.build_ssn_and_tree = build_ssn_and_tree_method
