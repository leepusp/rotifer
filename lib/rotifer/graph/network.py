"""
Build, view and analyse networks of protein domains found together.

Two domains are related when they occur in the same conserved gene
neighbourhood (a "block"): either fused in one protein, or in different
genes of the block. :func:`get_relations` turns a table of domain hits
into such relations (:func:`adjacency_edges` does the same for adjacent
features only), :func:`network_annotation` makes a styled networkx
graph from them, :func:`calculate_positions` lays the graph out for
:func:`rotifer.view.network.network_dash`, the interactive viewer, and
:func:`get_network_community`, :func:`community_stability`,
:func:`clique_jaccard_community` and :func:`subgraph_stats` find and
test groups of domains.

A typical session::

    edges = get_relations(domdf, gene_pos="start", feature_pos="dom_start")
    G = network_annotation(edges, node_color=colours)
    layouts = calculate_positions(G, iterations=50)
    positions = {name: {n: [float(c) for c in p] for n, p in layout.items()}
                 for name, layout in layouts.items()}
    schemes = get_network_community(G)
    # in rotifer.view.network:
    network_dash(G, positions, group_colours, node_groups, schemes)

The interactive viewer, :func:`rotifer.view.network.network_dash`, is in
:mod:`rotifer.view.network`, so this module does not need Dash. Nothing
here contacts a remote service.
"""
import networkx as nx
import numpy as np
import pandas as pd
import warnings

def network_annotation(networkdf,
                     node_color=None,
                     node_shape='ellipse',
                     second_shape=None,
                     second_shape_list=None,
                     node_outline=None,
                     node_outline_color=None,
                     node_outline_list=None,
                     edge_color=True,
                     unidirected_graph=True,
                     split_types=False):
    """
    Build a styled networkx graph from domain relations.

    Each domain becomes a node and each relation an edge. Nodes get a
    colour, shape, outline and degree, and edges get a line colour, style
    and opacity, as attributes that :func:`~rotifer.view.network.network_dash` and other drawing
    code read.

    Parameters
    ----------
    networkdf : pandas.DataFrame
        Relations, as returned by :func:`get_relations`: columns ``source``,
        ``target``, ``edge_type``, ``edge_count`` and ``edge_total``. When
        ``edge_total`` is missing it is estimated by summing ``edge_count``
        over the pair's edge types, which counts a block once for each type
        it shows, and a warning is raised.
    node_color : dict, optional
        Colour of each node, as a CSS colour (for example ``"#4C9BE8"``).
        Nodes not in it are drawn grey (``"#D3D3D3"``). By default every
        node is grey.
    node_shape : str, default 'ellipse'
        Cytoscape shape for the nodes, for example ``'ellipse'``,
        ``'rectangle'`` or ``'diamond'``.
    second_shape : dict or str, optional
        Shape for the nodes in ``second_shape_list``: a dict from node to
        shape, or one shape for all of them. By default they are hexagons.
    second_shape_list : list, optional
        Nodes drawn with ``second_shape`` instead of ``node_shape``.
    node_outline : bool, optional
        If true, outline the nodes in ``node_outline_list`` with
        ``node_outline_color``. Otherwise each outline has the node's own
        colour, so it does not show.
    node_outline_color : dict or str, optional
        Outline colour for the nodes in ``node_outline_list``: a dict from
        node to colour, or one colour for all of them. By default black.
    node_outline_list : list, optional
        Nodes that get an outline when ``node_outline`` is true.
    edge_color : bool, default True
        Style the edges: an edge between two nodes of the same colour is
        drawn in that colour, solid, with opacity 0.7; any other edge is
        grey, dotted, with opacity 0.2. If False, edges get no style
        attributes.
    unidirected_graph : bool, default True
        If true, build an undirected graph (see ``split_types``). If false,
        build a ``networkx.MultiDiGraph`` with one edge per row of
        ``networkdf``.
    split_types : bool, default False
        Undirected graphs only. If true, build a ``networkx.MultiGraph``
        with one edge per relation type of each pair, keyed by
        ``edge_type``. If false, build a ``networkx.Graph`` with one edge
        per pair, carrying the pair's most frequent type (ties broken
        alphabetically).

    Returns
    -------
    networkx.Graph, networkx.MultiGraph or networkx.MultiDiGraph
        The graph. Node attributes: ``node_color``, ``degree`` (number of
        neighbours), ``node_shape`` and ``node_outline``. Edge attributes:
        ``edge_type``, ``edge_count`` and ``edge_total`` from
        ``networkdf``; ``type_from`` in undirected graphs (see Notes); and
        ``line_style``, ``edge_color`` and ``opacity`` when ``edge_color``
        is true.

    Warns
    -----
    UserWarning
        If ``networkdf`` has no ``edge_total`` column, or if some nodes are
        missing from ``node_color``, ``second_shape`` or
        ``node_outline_color``. The warning names the nodes.

    See Also
    --------
    get_relations : Builds ``networkdf``.
    rotifer.view.network.network_dash : Interactive viewer for the returned graph.

    Notes
    -----
    In an undirected graph every relation is stored once. :func:`get_relations`
    lists each relation from both ends (A to B ``Cterminal`` and B to A
    ``Nterminal`` come from the same blocks), so the two rows are merged:
    each row is turned so that its source is the alphabetically first node,
    flipping the type when needed. Because an undirected edge does not keep
    a direction, ``type_from`` names the node that ``edge_type`` is read
    from. For example, ``type_from='ABC_tran'`` with ``edge_type='Downstream'``
    means the other domain's gene lies downstream of ABC_tran's gene.
    ``edge_total`` is the pair's number of blocks either way.
    """

    import networkx as nx
    import pandas as pd

    second_shape = {} if second_shape is None else second_shape
    second_shape_list = [] if second_shape_list is None else second_shape_list
    node_outline_color = {} if node_outline_color is None else node_outline_color
    node_outline_list = [] if node_outline_list is None else node_outline_list

    networkdf = networkdf.copy()
    if 'edge_total' not in networkdf.columns:
        warnings.warn("networkdf has no edge_total column; summing edge_count over edge types, "
                      "which counts a block once per type it occurs in")
        networkdf['edge_total'] = networkdf.groupby(['source', 'target'])['edge_count'].transform('sum')

    edge_attr = ['edge_type', 'edge_count', 'edge_total']
    if unidirected_graph:
        # One row per undirected relation: orient every row source <= target, flipping its type
        flip = {'Cterminal': 'Nterminal', 'Nterminal': 'Cterminal', 'Downstream': 'Upstream', 'Upstream': 'Downstream'}
        swap = networkdf.source > networkdf.target
        networkdf.loc[swap, ['source', 'target']] = networkdf.loc[swap, ['target', 'source']].values
        networkdf.loc[swap, 'edge_type'] = networkdf.loc[swap, 'edge_type'].replace(flip)
        # A self relation is its own mirror
        self_rel = networkdf.source == networkdf.target
        networkdf.loc[self_rel, 'edge_type'] = networkdf.loc[self_rel, 'edge_type'].replace({'Nterminal': 'Cterminal', 'Upstream': 'Downstream'})
        networkdf = networkdf.groupby(['source', 'target', 'edge_type'], as_index=False)[['edge_count', 'edge_total']].max()
        # An undirected edge does not keep its source: record which end edge_type refers to
        networkdf['type_from'] = networkdf.source
        edge_attr = edge_attr + ['type_from']
        if split_types:
            G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr, edge_key='edge_type', create_using=nx.MultiGraph())
        else:
            networkdf = networkdf.sort_values(['edge_count', 'edge_type'], ascending=[False, True]).drop_duplicates(['source', 'target'])
            G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr)
    else:
        G = nx.from_pandas_edgelist(networkdf, edge_attr=edge_attr, create_using=nx.MultiDiGraph())


    default_color = {item: "#D3D3D3" for item in G.nodes()}
    if node_color:
        missing = [node for node in G.nodes() if node not in node_color]
        if missing:
            warnings.warn(f"nodes missing from node_color, drawn grey: {sorted(missing)}")
        node_color = {**default_color, **node_color}
    else:
        node_color = default_color
    nx.set_node_attributes(G, node_color, name="node_color") 

    #Getting the size of nodes
    # In the future we can add more ways of getting the node size
    node_degree = dict(_simple_graph(G).degree())
    nx.set_node_attributes(G, node_degree, name="degree") 
    

    #Getting the shape of nodes
    #Creating a dictionary to shape nodes
    if second_shape == {}:
        second_shape={item: "hexagon" for item in second_shape_list}
    elif isinstance(second_shape, str):
        second_shape={item: second_shape for item in second_shape_list}
    else:
        pass

    missing = [node for node in G.nodes() if node in second_shape_list and node not in second_shape]
    if missing:
        warnings.warn(f"nodes missing from second_shape, drawn as {node_shape}: {sorted(missing)}")
    node_shape_dict  = {node: second_shape.get(node, node_shape) if node in second_shape_list else node_shape for node in G.nodes()}
    nx.set_node_attributes(G, node_shape_dict, name="node_shape")


    # Getting the outline of nodes:
    if node_outline_color == {}:
        node_outline_color={item: "#000000" for item in node_outline_list}
    elif isinstance(node_outline_color, str):
        node_outline_color={item: node_outline_color for item in node_outline_list}
    else:
        pass
    if node_outline:
        missing = [node for node in G.nodes() if node in node_outline_list and node not in node_outline_color]
        if missing:
            warnings.warn(f"nodes missing from node_outline_color, outlined in their node colour: {sorted(missing)}")
        node_outline_dict = {node: node_outline_color.get(node, node_color[node]) if node in node_outline_list else node_color[node] for node in G.nodes()}
    else:
        # The viewer always draws a border from node_outline: blend it into the node
        node_outline_dict = node_color
    nx.set_node_attributes(G, node_outline_dict, name="node_outline") 
    
    if edge_color:
        edge_colors=[]
        edge_alpha=[]
        edge_styles=[]
    #Getting the color of edges():
        if G.is_multigraph():
            edges_obj = G.edges(keys=True)
            for source,target, key in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
                   
            dict_edge_colors = { (u, v, key): color for (u, v, key), color in zip(edges_obj, edge_colors)}
            dict_edge_style = { (u, v, key): style for (u, v, key), style in zip(edges_obj, edge_styles)}
            dict_edge_alpha = { (u, v, key): alpha for (u, v, key), alpha in zip(edges_obj, edge_alpha)}

        else:    
            edges_obj = G.edges() 
            for source,target in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
            dict_edge_colors = {k: v for k, v in zip(G.edges(), edge_colors)}
            dict_edge_style = {k: v for k, v in zip(G.edges(), edge_styles)}
            dict_edge_alpha = {k: v for k, v in zip(G.edges(), edge_alpha)}

        nx.set_edge_attributes(G, dict_edge_style, 'line_style')
        nx.set_edge_attributes(G, dict_edge_colors, 'edge_color')
        nx.set_edge_attributes(G, dict_edge_alpha, 'opacity')        

    return G

def _simple_graph(G):
    # Parallel edges of an undirected MultiGraph (split_types) share one pair: collapse them
    # so degree, layouts and communities match the one-edge-per-pair graph
    """
    Collapse the parallel edges of an undirected multigraph.

    Parameters
    ----------
    G : networkx.Graph
        Any graph.

    Returns
    -------
    networkx.Graph
        ``G`` itself, unless it is an undirected multigraph (for example
        from ``network_annotation(..., split_types=True)``). Then a
        ``networkx.Graph`` with one edge per pair, so that degrees, layouts
        and communities match the graph with one edge per pair. The kept
        edge has the attributes of the pair's last parallel edge; in graphs
        from :func:`network_annotation` all of them carry the same
        ``edge_total``.
    """
    return nx.Graph(G) if G.is_multigraph() and not G.is_directed() else G

def calculate_positions(G, iterations):
    # Recalculate positions as needed
    """
    Compute four layouts of a graph for the viewer.

    :func:`~rotifer.view.network.network_dash` offers them as starting points: you can switch
    between them, move nodes by hand and save the result as a new layout.

    Parameters
    ----------
    G : networkx.Graph
        Graph to lay out. Parallel edges are collapsed first, so a pair
        with several relation types does not pull its nodes closer.
    iterations : int
        Iterations of the spring layout. More iterations let the layout
        settle further and take longer.

    Returns
    -------
    dict
        Layout name to ``{node: numpy.ndarray([x, y])}``. The names are
        ``'Spring Layout'`` (force-directed [1]_, seed 1, so it is
        reproducible), ``'Kamada-Kawai'`` (spring energy over graph
        distances [2]_), ``'Kamada-Kawai from Spring layout'`` (the same,
        starting from the spring layout) and ``'Circular'``.

    See Also
    --------
    rotifer.view.network.network_dash : Shows these layouts; the dict needs converting to
        lists first (see the example in that function).

    Notes
    -----
    The Kamada-Kawai layout uses the distances between all pairs of nodes,
    so its time and memory grow with the square of the number of nodes.

    References
    ----------
    .. [1] Fruchterman, T. M. J. & Reingold, E. M. "Graph drawing by
       force-directed placement." Softw. Pract. Exp. 21, 1129-1164 (1991).
       https://doi.org/10.1002/spe.4380211102
    .. [2] Kamada, T. & Kawai, S. "An algorithm for drawing general
       undirected graphs." Inf. Process. Lett. 31, 7-15 (1989).
       https://doi.org/10.1016/0020-0190(89)90102-6
    """
    G = _simple_graph(G)
    pos1 = nx.spring_layout(G, iterations=iterations, seed=1)
    pos2 = nx.kamada_kawai_layout(G, pos=pos1)
    pos3 = nx.kamada_kawai_layout(G)
    pos4 = nx.circular_layout(G)
    return {"Spring Layout": pos1,
            "Kamada-Kawai from Spring layout": pos2,
            "Kamada-Kawai": pos3,
            "Circular": pos4}

 

def _cpm_resolution(cpm_resolution, resolution_parameter):
    # resolution_parameter is the old name of cpm_resolution; it wins when given
    """
    Return the CPM resolution, accepting its deprecated name.

    Parameters
    ----------
    cpm_resolution : float
        Value given under the current name.
    resolution_parameter : float or None
        Value given under the old name, or None.

    Returns
    -------
    float
        ``resolution_parameter`` when it is given, otherwise
        ``cpm_resolution``.

    Warns
    -----
    DeprecationWarning
        When ``resolution_parameter`` is given. The warning points at the
        caller of the public function.
    """
    if resolution_parameter is None:
        return cpm_resolution
    warnings.warn("resolution_parameter is deprecated, use cpm_resolution "
                  "(it only sets the CPM resolution of the Leidein_2 partition)",
                  DeprecationWarning, stacklevel=3)
    return resolution_parameter

def _undirected_max(G, weight):
    # A pair's directed or parallel edges become one undirected edge with the largest weight
    # (net_functions graphs give every edge of a pair the same edge_total; summing would double it)
    """
    Turn a directed graph into an undirected one, keeping each pair's heaviest edge.

    Parameters
    ----------
    G : networkx.Graph
        Directed graph or multigraph.
    weight : str or None
        Edge attribute compared between the edges of a pair; edges without
        it count as 1. With None, the pair's first edge is kept.

    Returns
    -------
    networkx.Graph
        One edge per pair of nodes, with all the attributes of the kept
        edge. Node attributes are copied.

    Notes
    -----
    Graphs from :func:`network_annotation` give every edge of a pair the
    same ``edge_total``, so taking the largest weight keeps the pair's
    number of blocks, where a sum would count them twice.
    """
    H = nx.Graph()
    H.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if not H.has_edge(u, v) or (weight is not None and d.get(weight, 1) > H[u][v].get(weight, 1)):
            H.add_edge(u, v)
            H[u][v].clear()
            H[u][v].update(d)
    return H

def _community_scheme(nodes, communities, prefix):
    # [node -> label, label -> colour]: labels numbered by size, largest first (ties by smallest node name);
    # pastel has 10 colours and cycles, so more communities get evenly spaced hues
    """
    Turn a partition into a colour scheme for the viewer.

    Parameters
    ----------
    nodes : iterable
        Nodes in the order the returned mapping should follow. Every node
        must belong to one of ``communities``.
    communities : iterable of set
        The communities, as sets of nodes.
    prefix : str
        Text put before each community number, for example ``'Louvain_'``.

    Returns
    -------
    list
        ``[node_to_label, label_to_colour]``. Communities are numbered by
        size, 0 being the largest; ties go to the community whose smallest
        node name comes first. Colours come from the seaborn ``pastel``
        palette for up to 10 communities and from evenly spaced ``husl``
        hues beyond that, so no two communities share a colour.
    """
    import seaborn as sns
    ordered = sorted(communities, key=lambda c: (-len(c), min(str(n) for n in c)))
    label = {n: f"{prefix}{i}" for i, c in enumerate(ordered) for n in c}
    palette = sns.color_palette('pastel' if len(ordered) <= 10 else 'husl', len(ordered)).as_hex()
    return [{n: label[n] for n in nodes}, {f"{prefix}{i}": colour for i, colour in enumerate(palette)}]

def get_network_community(G,
                 community_to_color='leidein', weight='edge_total', cpm_resolution=0.05, seed=0,
                 *, resolution_parameter=None):
    """
    Find communities with Louvain, Leiden and Leiden-CPM, as viewer colour schemes.

    Partitions ``G`` three ways: Louvain [1]_ and Leiden [2]_, both
    maximising standard modularity [3]_ (resolution 1), and Leiden with the
    Constant Potts Model (CPM) [4]_ at resolution ``cpm_resolution``. Each
    partition is returned as a colour scheme that :func:`~rotifer.view.network.network_dash` can
    show.

    Parameters
    ----------
    G : networkx.Graph
        Graph to partition. A directed graph is partitioned as an
        undirected copy in which each pair keeps its heaviest edge. The
        parallel edges of a multigraph are collapsed to one edge per pair.
    community_to_color : str, default 'leidein'
        Unused. Kept so that positional calls still work; any other value
        raises a DeprecationWarning.
    weight : str or None, default 'edge_total'
        Edge attribute used as the weight, or None for an unweighted run.
        Edges without the attribute count as 1.
    cpm_resolution : float, default 0.05
        Resolution of ``Leidein_2``. A group of nodes becomes a community
        only if its edges weigh more than this per pair of nodes, on
        average, so higher values give more and smaller communities. In
        the unit of ``weight``: with ``edge_total``, a number of blocks.
        Louvain and ``Leidein`` do not use it.
    seed : int or None, default 0
        Random seed for all three algorithms, so a rerun gives the same
        partitions. None gives unseeded runs.
    resolution_parameter : float, optional
        Deprecated name of ``cpm_resolution``. When given, it is used, and
        a DeprecationWarning is raised.

    Returns
    -------
    dict
        ``{"Louvain" | "Leidein" | "Leidein_2": [node_to_label,
        label_to_colour]}``. Labels are ``Louvain_<i>``, ``Leiden_<i>`` and
        ``LeidenCPM_<i>``, numbered by community size (0 is the largest,
        ties broken by the smallest node name). Colours are pastel for up
        to 10 communities and evenly spaced hues beyond that. An empty
        graph gives empty dicts.

    Warns
    -----
    UserWarning
        If ``G`` is empty or directed, or if some or all edges lack
        ``weight``.
    DeprecationWarning
        If ``community_to_color`` or ``resolution_parameter`` is passed.

    See Also
    --------
    community_stability : How much these partitions depend on ``seed``.
    clique_jaccard_community : Communities of overlapping maximal cliques.
    rotifer.view.network.network_dash : Shows the returned colour schemes.

    Notes
    -----
    A fixed seed makes a partition reproducible, yet another seed can
    still give another partition. Modularity has many structurally
    different partitions with nearly the same score [5]_, so this happens
    on graphs with weak community structure.
    :func:`community_stability` measures how much that happens.

    References
    ----------
    .. [1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E.
       "Fast unfolding of communities in large networks." J. Stat. Mech.
       P10008 (2008). https://doi.org/10.1088/1742-5468/2008/10/P10008
    .. [2] Traag, V. A., Waltman, L. & van Eck, N. J. "From Louvain to
       Leiden: guaranteeing well-connected communities." Sci. Rep. 9, 5233
       (2019). https://doi.org/10.1038/s41598-019-41695-z
    .. [3] Newman, M. E. J. & Girvan, M. "Finding and evaluating community
       structure in networks." Phys. Rev. E 69, 026113 (2004).
       https://doi.org/10.1103/PhysRevE.69.026113
    .. [4] Traag, V. A., Van Dooren, P. & Nesterov, Y. "Narrow scope for
       resolution-limit-free community detection." Phys. Rev. E 84, 016114
       (2011). https://doi.org/10.1103/PhysRevE.84.016114
    .. [5] Good, B. H., de Montjoye, Y.-A. & Clauset, A. "Performance of
       modularity maximization in practical contexts." Phys. Rev. E 81,
       046106 (2010). https://doi.org/10.1103/PhysRevE.81.046106
    """
    import leidenalg as la
    import igraph as ig

    cpm_resolution = _cpm_resolution(cpm_resolution, resolution_parameter)
    if community_to_color != 'leidein':
        warnings.warn("community_to_color is unused and deprecated; all three partitions are always returned",
                      DeprecationWarning, stacklevel=2)
    if G.number_of_nodes() == 0:
        warnings.warn("G has no nodes; returning empty community schemes.")
        return {name: [{}, {}] for name in ["Louvain", "Leidein", "Leidein_2"]}
    if G.is_directed():
        warnings.warn("G is directed; detecting communities on an undirected copy (a pair's edges keep the largest weight).")
        G = _undirected_max(G, weight)
    else:
        G = _simple_graph(G)

    #Edge weights: weight=None means unweighted; edges lacking the attribute count as 1
    g = ig.Graph.from_networkx(G)
    leiden_weight = weight
    if weight is not None:
        n_missing = sum(1 for u, v, d in G.edges(data=True) if weight not in d)
        if G.number_of_edges() == 0:
            weight = leiden_weight = None
        elif n_missing == G.number_of_edges():
            warnings.warn(f"No edge has the '{weight}' attribute; running community detection unweighted.")
            weight = leiden_weight = None
        elif n_missing:
            warnings.warn(f"{n_missing} of {G.number_of_edges()} edges lack the '{weight}' attribute; using 1 for them.")
            leiden_weight = [1 if w is None else w for w in g.es[weight]]

    #Partitions as sets of node names, never by position
    names = g.vs['_nx_name']
    def leiden(partition_type, **kwargs):
        part = la.find_partition(g, partition_type, weights=leiden_weight, n_iterations=-1, seed=seed, **kwargs)
        return [{names[i] for i in c} for c in part]

    louvain = nx.community.louvain_communities(G, weight=weight, resolution=1, seed=seed)
    return {"Louvain": _community_scheme(G, louvain, 'Louvain_'),
            "Leidein": _community_scheme(G, leiden(la.ModularityVertexPartition), 'Leiden_'),
            #Leiden with CPM, the only one using the resolution:
            "Leidein_2": _community_scheme(G, leiden(la.CPMVertexPartition, resolution_parameter=cpm_resolution), 'LeidenCPM_')}

def community_stability(G, weight='edge_total', cpm_resolution=0.05, n_seeds=50, seed=0, warn_below=0.8):
    """
    Measure how much the partitions of get_network_community depend on the seed.

    Runs :func:`get_network_community` with the seeds ``seed``,
    ``seed + 1``, ..., ``seed + n_seeds - 1``. The run with ``seed`` is the
    reference: it gives the partition you would report. Modularity has many
    structurally different partitions with nearly the same score [1]_, and
    Louvain [2]_ and Leiden [3]_ settle on one of them depending on the
    seed. This function reports how much the result changes with the seed,
    for each partition and for each node.

    Parameters
    ----------
    G : networkx.Graph
        Graph to partition, as in :func:`get_network_community`.
    weight : str or None, default 'edge_total'
        Edge attribute used as the weight, or None for an unweighted run.
    cpm_resolution : float, default 0.05
        Resolution of the ``Leidein_2`` partition.
    n_seeds : int, default 50
        Number of runs, the reference included.
    seed : int, default 0
        Seed of the reference run. The other runs use the seeds after it.
    warn_below : float or None, default 0.8
        Warn for each partition whose median ARI is below this value. None
        turns the warning off.

    Returns
    -------
    summary : pandas.DataFrame
        One row per partition (``Louvain``, ``Leidein``, ``Leidein_2``),
        with columns ``communities`` (communities in the reference),
        ``distinct_partitions`` (different partitions among the runs), and
        ``ari_median`` and ``ari_min``: the median and smallest adjusted
        Rand index [4]_ over all pairs of runs. An ARI of 1 means identical
        partitions; about 0 means no more agreement than chance.
    nodes : pandas.DataFrame
        One row per node with two columns per partition: ``<name>``, the
        node's community in the reference, and ``<name>_stability``, the
        mean Jaccard index [5]_ between the node's community in each other
        run and its reference community. A stability of 1 means the node
        always had the same community; splits and merges both lower it.

    Raises
    ------
    ValueError
        If ``n_seeds`` is less than 2.

    Warns
    -----
    UserWarning
        For each partition whose median ARI is below ``warn_below``.
        Warnings from the repeated :func:`get_network_community` calls,
        such as a directed graph or missing weights, are shown once each.

    See Also
    --------
    get_network_community : The partitions being measured.

    Notes
    -----
    The node stability follows the co-assignment idea of consensus
    clustering [6]_, which counts how often nodes end up together across
    runs. Here it is measured against the reference, which avoids building
    a node by node matrix and costs time proportional to the number of
    nodes per run.

    This measures stability across seeds only. It does not tell how robust
    the communities are to sampling of the data, nor whether the graph has
    significant community structure. The runtime is about ``n_seeds`` calls
    of :func:`get_network_community`.

    References
    ----------
    .. [1] Good, B. H., de Montjoye, Y.-A. & Clauset, A. "Performance of
       modularity maximization in practical contexts." Phys. Rev. E 81,
       046106 (2010). https://doi.org/10.1103/PhysRevE.81.046106
    .. [2] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E.
       "Fast unfolding of communities in large networks." J. Stat. Mech.
       P10008 (2008). https://doi.org/10.1088/1742-5468/2008/10/P10008
    .. [3] Traag, V. A., Waltman, L. & van Eck, N. J. "From Louvain to
       Leiden: guaranteeing well-connected communities." Sci. Rep. 9, 5233
       (2019). https://doi.org/10.1038/s41598-019-41695-z
    .. [4] Hubert, L. & Arabie, P. "Comparing partitions." J. Classif. 2,
       193-218 (1985). https://doi.org/10.1007/BF01908075
    .. [5] Jaccard, P. "The distribution of the flora in the alpine zone."
       New Phytol. 11, 37-50 (1912).
       https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
    .. [6] Lancichinetti, A. & Fortunato, S. "Consensus clustering in
       complex networks." Sci. Rep. 2, 336 (2012).
       https://doi.org/10.1038/srep00336
    """
    import igraph as ig
    from collections import Counter

    if n_seeds < 2:
        raise ValueError(f"n_seeds must be >= 2 to compare runs, got {n_seeds}")
    # The repeated calls repeat their warnings (directed graph, missing weights): show each once
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        runs = [get_network_community(G, weight=weight, cpm_resolution=cpm_resolution, seed=s)
                for s in range(seed, seed + n_seeds)]
    for message, category in dict.fromkeys((str(w.message), w.category) for w in caught):
        warnings.warn(message, category, stacklevel=2)

    names = ["Louvain", "Leidein", "Leidein_2"]
    nodes = list(runs[0]["Louvain"][0])
    summary = pd.DataFrame(index=pd.Index(names, name='partition'),
                           columns=['communities', 'distinct_partitions', 'ari_median', 'ari_min'], dtype=float)
    node_table = pd.DataFrame(index=pd.Index(nodes, name='node'))
    if not nodes:
        summary[['communities', 'distinct_partitions']] = [0, 1]
        for name in names:
            node_table[name], node_table[f"{name}_stability"] = pd.Series(dtype=object), pd.Series(dtype=float)
        return summary, node_table

    for name in names:
        # Codes by first appearance in node order: relabelled copies of a partition get the same codes
        memberships = [pd.factorize(pd.Series([run[name][0][n] for n in nodes]))[0] for run in runs]
        counts = Counter(m.tobytes() for m in memberships)
        partitions = {m.tobytes(): m for m in memberships}

        # ARI over all pairs of runs, computed once per pair of distinct partitions
        aris, pairs = [], []
        keys = list(counts)
        for i, a in enumerate(keys):
            if counts[a] > 1:
                aris.append(1.0)
                pairs.append(counts[a] * (counts[a] - 1) // 2)
            for b in keys[i + 1:]:
                aris.append(ig.compare_communities(partitions[a].tolist(), partitions[b].tolist(), method='adjusted_rand'))
                pairs.append(counts[a] * counts[b])
        all_aris = np.repeat(aris, pairs)

        # Node stability: Jaccard of its community in each run with its reference community
        ref = memberships[0]
        ref_size = np.bincount(ref)[ref]
        stability = np.zeros(len(nodes))
        for m in memberships[1:]:
            pair = ref * (m.max() + 1) + m
            shared = np.bincount(pair)[pair]
            stability += shared / (ref_size + np.bincount(m)[m] - shared)

        summary.loc[name] = [ref.max() + 1, len(counts), np.median(all_aris), all_aris.min()]
        node_table[name] = [runs[0][name][0][n] for n in nodes]
        node_table[f"{name}_stability"] = stability / (n_seeds - 1)
        if warn_below is not None and np.median(all_aris) < warn_below:
            warnings.warn(f"{name}: median ARI {np.median(all_aris):.2f} across {n_seeds} seeds is below {warn_below}; "
                          f"its communities depend on the seed (the nodes table shows which nodes are unstable).",
                          stacklevel=2)
    summary[['communities', 'distinct_partitions']] = summary[['communities', 'distinct_partitions']].astype(int)
    return summary, node_table

def _clique_pairs(cliques, nodes, min_jaccard_index, warn_checked, warn_kept, block=1000):
    # Pairs of cliques sharing at least one node, from a sparse clique x node matrix times its transpose,
    # one block of cliques at a time so memory stays bounded; each pair once (i < j), no clique with itself.
    # Returns positions i, j (sorted) and their Jaccard index.
    """
    Find the pairs of cliques that share nodes, with their Jaccard index.

    The shared-node count of every pair comes from a sparse clique by node
    matrix multiplied by its transpose, one block of cliques at a time, so
    memory stays bounded. Pairs that share no node are never formed.

    Parameters
    ----------
    cliques : list of list
        Nodes of each clique.
    nodes : list
        All nodes of the graph.
    min_jaccard_index : float
        Pairs below this Jaccard index are dropped.
    warn_checked : int or None
        Warn when the estimated number of pairs to check is above this.
        None turns the warning off.
    warn_kept : int or None
        Warn when more pairs than this are kept. None turns the warning
        off.
    block : int, default 1000
        Number of cliques compared at a time.

    Returns
    -------
    i, j : numpy.ndarray
        Positions in ``cliques`` of each kept pair, with ``i < j``, sorted
        by ``i`` then ``j``. A clique is never paired with itself.
    jaccard : numpy.ndarray
        Jaccard index of each kept pair.

    Warns
    -----
    UserWarning
        After the first block, if the overlap rate seen there, extrapolated
        to all pairs, gives more than ``warn_checked`` pairs to check; and
        at the end, if more than ``warn_kept`` pairs are kept.
    """
    from scipy import sparse
    node_id = {n: k for k, n in enumerate(nodes)}
    size = np.array([len(c) for c in cliques], dtype=np.int64)
    n = len(cliques)
    rows = np.repeat(np.arange(n), size)
    cols = np.fromiter((node_id[v] for c in cliques for v in c), dtype=np.int64, count=size.sum())
    M = sparse.csr_matrix((np.ones(len(cols), dtype=np.float32), (rows, cols)), shape=(n, len(node_id)))
    MT = M.T.tocsr()
    kept_i, kept_j, kept_jaccard = [], [], []
    for start in range(0, n, block):
        part = (M[start:start + block] @ MT).tocoo()
        i = part.row.astype(np.int64) + start
        j = part.col.astype(np.int64)
        upper = j > i
        i, j, shared = i[upper], j[upper], part.data[upper].astype(np.float64)
        if start == 0 and warn_checked is not None:
            # Overlap rate in the first block, extrapolated to all pairs
            block_pairs = sum(n - 1 - r for r in range(min(block, n)))
            estimate = len(i) / block_pairs * n * (n - 1) / 2 if block_pairs else 0
            if estimate > warn_checked:
                warnings.warn(f"about {estimate:,.0f} pairs of the {n:,} cliques share nodes and must be checked "
                              f"(more than warn_checked={warn_checked:,}); this can take minutes. "
                              "Raise min_clique_size or min_jaccard_index to reduce it.", stacklevel=3)
        jaccard = shared / (size[i] + size[j] - shared)
        keep = jaccard >= min_jaccard_index
        order = np.lexsort((j[keep], i[keep]))
        kept_i.append(i[keep][order]); kept_j.append(j[keep][order]); kept_jaccard.append(jaccard[keep][order])
    i, j, jaccard = (np.concatenate(x) if x else np.array([], dtype=np.int64) for x in (kept_i, kept_j, kept_jaccard))
    if warn_kept is not None and len(i) > warn_kept:
        warnings.warn(f"{len(i):,} pairs of cliques pass min_jaccard_index (more than warn_kept={warn_kept:,}); "
                      "the result table will be very large. Raise min_clique_size or min_jaccard_index to reduce it.",
                      stacklevel=3)
    return i, j, jaccard

def clique_jaccard_community(G,
                             min_jaccard_index=0.4,
                             cpm_resolution=0.5,
                             min_clique_size=10,
                             seed=0,
                             *, resolution_parameter=None,
                             warn_checked=100_000_000,
                             warn_kept=5_000_000):
    """
    Find communities of overlapping maximal cliques.

    Finds the maximal cliques of ``G`` with at least ``min_clique_size``
    nodes [1]_ [2]_, connects each pair of cliques whose Jaccard index [3]_
    (nodes in both over nodes in either) is at least ``min_jaccard_index``,
    and partitions that graph of cliques with
    :func:`get_network_community`, using the Jaccard index as the edge
    weight. The idea resembles clique percolation [4]_, which joins
    k-cliques sharing k - 1 nodes; here the cliques are maximal and are
    joined by their Jaccard index.

    Parameters
    ----------
    G : networkx.Graph
        Graph whose maximal cliques are compared.
    min_jaccard_index : float, default 0.4
        Smallest Jaccard index for two cliques to be connected. Must be
        greater than 0.
    cpm_resolution : float, default 0.5
        Resolution of the ``Leidein_2`` partition (the
        ``*_leidein_resolution`` columns), passed to
        :func:`get_network_community`.
    min_clique_size : int, default 10
        Smallest clique compared, in nodes.
    seed : int or None, default 0
        Random seed passed to :func:`get_network_community`.
    resolution_parameter : float, optional
        Deprecated name of ``cpm_resolution``.
    warn_checked : int or None, default 100_000_000
        Warn when the estimated number of clique pairs that share a node,
        all of which must be checked, is above this. The estimate is made
        after the first 1,000 cliques. None turns the warning off.
    warn_kept : int or None, default 5_000_000
        Warn when more pairs than this pass ``min_jaccard_index``, which
        means the returned table will be very large. None turns the warning
        off.

    Returns
    -------
    pandas.DataFrame
        One row per connected pair of cliques, with ``source < target``.
        Columns:

        - ``source``, ``target``: clique numbers, their position in the
          output of ``networkx.find_cliques``.
        - ``cliques_x``, ``clique_size_x``, ``cliques_y``,
          ``clique_size_y``: the nodes and size of each clique.
        - ``t``: always 1.
        - ``uni``, ``inter``: the sets of nodes in either clique and in
          both.
        - ``Jaccard_index``.
        - ``source_leidein_resolution``, ``target_leidein_resolution``,
          ``source_leidein``, ``target_leidein``, ``source_louvain``,
          ``target_louvain``: the community of each clique in the
          ``Leidein_2``, ``Leidein`` and ``Louvain`` partitions.

        Empty, with the same columns, when no pair qualifies.

    Raises
    ------
    ValueError
        If ``min_jaccard_index`` is not greater than 0.

    Warns
    -----
    UserWarning
        If the estimated pairs to check exceed ``warn_checked``, or the
        kept pairs exceed ``warn_kept``.

    See Also
    --------
    get_network_community : The partitions computed on the graph of cliques.

    Notes
    -----
    Pairs come from the number of nodes each pair of cliques shares,
    computed as a sparse clique by node matrix times its transpose, 1,000
    cliques at a time. Pairs that share no node (Jaccard index 0) are never
    formed, which is why ``min_jaccard_index`` must be positive, and a
    clique that overlaps no other clique above the threshold is absent from
    the result.

    The cost grows with the number of pairs that share a node, and memory
    is bounded by the block of cliques. When nearly every pair overlaps,
    as with many maximal cliques on few nodes, the work grows with the
    square of the number of cliques again; raising ``min_clique_size`` or
    ``min_jaccard_index`` reduces it.

    References
    ----------
    .. [1] Bron, C. & Kerbosch, J. "Algorithm 457: finding all cliques of
       an undirected graph." Commun. ACM 16, 575-577 (1973).
       https://doi.org/10.1145/362342.362367
    .. [2] Tomita, E., Tanaka, A. & Takahashi, H. "The worst-case time
       complexity for generating all maximal cliques and computational
       experiments." Theor. Comput. Sci. 363, 28-42 (2006).
       https://doi.org/10.1016/j.tcs.2006.06.015
    .. [3] Jaccard, P. "The distribution of the flora in the alpine zone."
       New Phytol. 11, 37-50 (1912).
       https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
    .. [4] Palla, G., Derenyi, I., Farkas, I. & Vicsek, T. "Uncovering the
       overlapping community structure of complex networks in nature and
       society." Nature 435, 814-818 (2005).
       https://doi.org/10.1038/nature03607
    """
    cpm_resolution = _cpm_resolution(cpm_resolution, resolution_parameter)
    if not min_jaccard_index > 0:
        raise ValueError(f"min_jaccard_index must be > 0 (pairs sharing no node are never connected), got {min_jaccard_index}")
    cliques = list(nx.find_cliques(G))
    cliques = pd.Series(cliques, dtype=object).rename('cliques').to_frame()
    cliques['clique_size'] = cliques.cliques.apply(lambda x: len(x))
    cliques = cliques.query('clique_size >= @min_clique_size')
    cliques['t']  = 1
    cliques = cliques.reset_index().rename({'index':'clique_n'}, axis=1)
    i, j, jaccard = _clique_pairs(cliques.cliques.tolist(), list(G), min_jaccard_index, warn_checked, warn_kept)
    # Same columns, order and index as the former all-pairs merge of cliques with itself
    left = cliques.iloc[i].reset_index(drop=True).add_suffix('_x').rename({'t_x': 't'}, axis=1)
    right = cliques.iloc[j].reset_index(drop=True).add_suffix('_y').drop(columns='t_y')
    cm = pd.concat([left, right], axis=1)
    cm.index = i * len(cliques) + j
    cm['uni'] = [set(a) | set(b) for a, b in zip(cm.cliques_x, cm.cliques_y)]
    cm['inter'] = [set(a) & set(b) for a, b in zip(cm.cliques_x, cm.cliques_y)]
    cm['Jaccard_index'] = jaccard.astype(float)
    cm.rename({"clique_n_x" :'source', 'clique_n_y' :'target'}, axis=1, inplace=True)
    community_columns = ['source_leidein_resolution', 'target_leidein_resolution', 'source_leidein',
                         'target_leidein', 'source_louvain', 'target_louvain']
    if cm.empty:
        return cm.assign(**{c: pd.Series(dtype=object) for c in community_columns})
    G = nx.from_pandas_edgelist(cm, edge_attr='Jaccard_index')
    z = get_network_community(G, weight='Jaccard_index', cpm_resolution=cpm_resolution, seed=seed)
    cm['source_leidein_resolution'] = cm.source.replace(z['Leidein_2'][0])
    cm['target_leidein_resolution'] = cm.target.replace(z['Leidein_2'][0])
    cm['source_leidein'] = cm.source.replace(z['Leidein'][0])
    cm['target_leidein'] = cm.target.replace(z['Leidein'][0])
    cm['source_louvain'] = cm.source.replace(z['Louvain'][0])
    cm['target_louvain'] = cm.target.replace(z['Louvain'][0])
    return cm

def subgraph_stats(net,
                   subnet,
                   theme_dict,
                   theme_to_stats,
                   alternative = "two-sided",
                   return_df = False):
    """
    Test whether a theme is enriched in a subnetwork.

    Compares how often nodes of one theme, such as a functional category,
    occur in a subnetwork and in the rest of the network, with Fisher's
    exact test [1]_ on the 2x2 table::

                           theme   not theme
        subnetwork         a       b
        net - subnet       c       d

    "Not theme" means annotated with another theme. Nodes without an
    annotation are left out of the table (see Notes).

    Parameters
    ----------
    net : networkx.Graph
        The whole network. Its nodes outside ``subnet`` form the second
        row of the table.
    subnet : networkx.Graph
        The subnetwork tested, usually a subgraph of ``net`` such as a
        community.
    theme_dict : dict or pandas.Series
        Theme of each node. Nodes missing from it, or mapped to a missing
        value (``None``, NaN), count as unannotated.
    theme_to_stats : object
        The theme tested, compared with ``==`` to the values of
        ``theme_dict``.
    alternative : {'two-sided', 'greater', 'less'}, default 'two-sided'
        Alternative hypothesis, passed to ``scipy.stats.fisher_exact``.
        ``'greater'`` tests for enrichment in the subnetwork, ``'less'``
        for depletion.
    return_df : bool, default False
        Also return the table of counts.

    Returns
    -------
    result : SignificanceResult
        The result of ``scipy.stats.fisher_exact``, with ``statistic`` (the
        sample odds ratio, ``a*d / (b*c)``) and ``pvalue``. Returned alone
        when ``return_df`` is false.
    table : pandas.DataFrame
        Only when ``return_df`` is true, returned as ``(result, table)``.
        Rows ``'Subnetwork'`` and ``'net - subnet'``; columns ``'Theme'``
        and ``'Not theme'`` (the counts tested) and ``'Unannotated'`` (the
        nodes left out of the test).

    Warns
    -----
    UserWarning
        When nodes are left out because they have no annotation, with
        their number in each row.

    See Also
    --------
    get_network_community : Partitions whose communities can be tested.
    community_stability : Whether those communities depend on the seed.

    Notes
    -----
    Leaving unannotated nodes out tests enrichment among annotated nodes.
    That is unbiased when annotation is missing at random. When it is not,
    for example when well-studied themes are better annotated, the result
    can still be biased, and the warning is the cue to check. Fewer nodes
    also mean less power.

    A node is counted once per row, however often it is listed. Nodes of
    ``subnet`` that are not in ``net`` are still counted in the subnetwork
    row.

    References
    ----------
    .. [1] Fisher, R. A. "On the interpretation of chi-squared from
       contingency tables, and the calculation of P." J. R. Stat. Soc. 85,
       87-94 (1922). https://doi.org/10.2307/2340521
    """
    from  scipy.stats import fisher_exact
    import pandas as pd

    def count(nodes):
        # Nodes missing from theme_dict, or mapped to a missing value, are left out of the test
        annotated = [n for n in nodes if n in theme_dict and not (np.ndim(theme_dict[n]) == 0 and pd.isna(theme_dict[n]))]
        theme = sum(1 for n in annotated if theme_dict[n] == theme_to_stats)
        return theme, len(annotated) - theme, len(nodes) - len(annotated)

    sub_nodes = set(subnet.nodes)
    rest_nodes = set(net.nodes) - sub_nodes
    sub_theme, sub_other, sub_missing = count(sub_nodes)
    net_theme, net_other, net_missing = count(rest_nodes)
    if sub_missing or net_missing:
        warnings.warn(f"{sub_missing} of {len(sub_nodes)} subnetwork nodes and {net_missing} of {len(rest_nodes)} "
                      f"other nodes are missing from theme_dict and were left out of the test", stacklevel=2)
    result = fisher_exact([[sub_theme, sub_other], [net_theme, net_other]], alternative=alternative)
    if return_df:
        table = pd.DataFrame([[sub_theme, sub_other, sub_missing], [net_theme, net_other, net_missing]],
                             columns=['Theme', 'Not theme', 'Unannotated'], index=['Subnetwork', 'net - subnet'])
        return result, table
    return result

def count_series(
        series,
        normalize=False,
        cut_off=False,
        count='domain',
        as_list=False
        ):
    """
    Count the domains or architectures in a series and list them with their counts.

    Each value is an architecture: domain names joined by ``+``, such as
    ``"SP+TM"``. By default every architecture is split into its domains
    and the domains are counted; with ``count='architecture'`` whole
    architectures are counted instead. :func:`get_relations` uses it to
    list the partners of each domain in its YAML output.

    Parameters
    ----------
    series : pandas.Series of str
        Architectures, one per row. Missing values are skipped.
    normalize : bool, default False
        Report each item's frequency as a percentage instead of its count.
    cut_off : float or False, default False
        Keep only items whose count is at least this value. With
        ``normalize=True`` the comparison uses the frequency as a fraction
        between 0 and 1, so pass 0.1 for 10 %. False keeps everything.
    count : {'domain', 'architecture'}, default 'domain'
        What to count: single domains or whole architectures.
    as_list : bool, default False
        Return a list of strings instead of one comma-separated string.

    Returns
    -------
    str or list of str
        Items formatted as ``"name(count)"``, or ``"name(12.34%)"`` with
        ``normalize=True``, most frequent first. A string joins them with
        ``", "`` (empty when there is nothing to count); a list has one item
        per entry, or ``[""]`` when there is nothing to count.

    Notes
    -----
    Moved unchanged from ``rotifer.devel.alpha.gian_func`` (issue #150),
    which still imports it from here. Its old docstring said ``cut_off`` is
    a percentage when ``normalize=True``; the code compares it with the
    fraction, as described above.

    Examples
    --------
    >>> import pandas as pd
    >>> count_series(pd.Series(["SP+TM", "TM", None]))
    'TM(2), SP(1)'
    >>> count_series(pd.Series(["SP+TM", "TM", "TM"]), count="architecture", as_list=True)
    ['TM(2)', 'SP+TM(1)']
    """
    flattened_list = []
    if count == 'architecture':
        s = series.dropna().value_counts(normalize=normalize)
    else:
        s = series.str.split('+').explode().dropna().value_counts(
                normalize=normalize
                )

    if cut_off:
        #cut_off = cut_off/100
        s = s.where(lambda x: x >= cut_off).dropna()
    for y, z in s.items():
        if normalize:
            flattened_list.append(f'{y}({100 * z:.2f}%)')
        else:
            flattened_list.append(f'{y}({z})')
    
    if as_list:
        if len(flattened_list) ==0:
            return list([''])
        else:
            return flattened_list

    return ', '.join(flattened_list)

def get_relations(domdf,
                  dom='full',
                  min_connections = 5,
                  self_relation = False,
                  to_yaml=False,
                  iteration=1,
                  filter_rename_yaml  = False,
                  *,
                  gene_pos,
                  feature_pos):
    """
    List the relations between domains that share gene neighbourhoods.

    Two domains are related when they occur in the same block (a conserved
    gene neighbourhood). Domains of the same protein are related by their
    order along the protein (a fusion); domains of different genes by the
    order of the genes in the block. Each relation is counted once per
    block, and only pairs seen in at least ``min_connections`` blocks are
    kept.

    Parameters
    ----------
    domdf : pandas.DataFrame
        One row per domain hit, with columns ``block_id`` (the
        neighbourhood), ``pid`` (the protein), ``dom`` (the domain name),
        and the columns named by ``gene_pos`` and ``feature_pos``. It is
        not modified.
    dom : str or list of str, default 'full'
        ``'full'`` uses every block. Otherwise, one or more seed domains:
        only blocks that contain them are used (see ``iteration``).
    min_connections : int, default 5
        Smallest number of blocks in which a pair of domains must occur,
        in any relation, for its relations to be kept.
    self_relation : bool, default False
        Keep relations of a domain with itself, for example two copies in
        one protein or in neighbouring genes.
    to_yaml : bool, default False
        Return a dict ready to be written as YAML instead of a DataFrame.
    iteration : int, default 1
        With seed domains, how far to expand from them. With 1, the blocks
        used are those containing a seed domain. With 2, the blocks
        containing any domain found in those blocks, and so on. Must be at
        least 1.
    filter_rename_yaml : str or False, default False
        Path to a YAML file that renames and selects domains. Each key is
        a domain name with two fields: ``Display_name`` (the new name; an
        empty string keeps the key) and ``Include`` (1 to keep the domain).
        Domains are renamed first; then only domains whose display name
        belongs to an included entry are kept, so domains missing from the
        file are dropped.
    gene_pos : str
        Keyword only. Column that orders the genes of a block along the
        genome, smaller meaning upstream, for example the start
        coordinate.
    feature_pos : str
        Keyword only. Column that orders the domains of a protein from the
        N- to the C-terminus, smaller meaning more N-terminal, for example
        the domain start in amino-acid coordinates.

    Returns
    -------
    pandas.DataFrame or dict
        With ``to_yaml`` false, a DataFrame with one row per (``source``,
        ``target``, ``edge_type``) and columns:

        - ``source``, ``target``: the two domains.
        - ``edge_type``: where ``target`` lies relative to ``source``.
          ``'Cterminal'`` or ``'Nterminal'`` in the same protein,
          ``'Downstream'`` or ``'Upstream'`` in another gene of the block.
        - ``edge_count``: number of blocks with this relation.
        - ``edge_total``: number of blocks with the pair in any relation,
          the value ``min_connections`` is compared to.

        Every relation appears from both ends: A to B ``Cterminal`` comes
        with B to A ``Nterminal``, from the same blocks. Row order carries
        no meaning.

        With ``to_yaml`` true, a dict ``{source: {'total': n, 'Cterminal':
        [...], 'Downstream': [...], 'Nterminal': [...], 'Upstream':
        [...]}}``, where ``total`` is the number of blocks in which
        ``source`` has a kept relation, and each list holds
        ``"target(count)"`` strings, or ``['']`` when there is none. An
        empty dict when no pair passes ``min_connections``.

    Raises
    ------
    ValueError
        If seed domains are given and ``iteration`` is less than 1.

    Warns
    -----
    UserWarning
        If some seed domains do not occur in ``domdf``.

    See Also
    --------
    network_annotation : Builds a graph from the returned DataFrame.
    adjacency_edges : Edges between adjacent features only, from a table
        already in genome order.

    Notes
    -----
    Relation types come from the ``gene_pos`` and ``feature_pos``
    values, so row order does not matter. Two domains with the same
    ``feature_pos`` in one protein, or two genes with the same
    ``gene_pos``, get no relation.

    ``gene_pos`` is used as it is, whatever the strand. For the output of
    ``gian_func.split_domain``, whose ``domp`` follows genome order (C to N
    on the minus strand), use ``gene_pos='blockp'`` and a ``feature_pos``
    column built as::

        domdf['dom_aa_order'] = np.where(domdf.strand == -1, -domdf.domp, domdf.domp)

    To keep only the relations among the seed domains themselves, filter
    the result: ``edges[edges.source.isin(dom) & edges.target.isin(dom)]``.
    """
    domdf = domdf.copy()
    domdf.dom = domdf.dom.str.strip()
    if filter_rename_yaml:
        import yaml
        with open(filter_rename_yaml, "r") as file:
            data = yaml.safe_load(file)
            data = pd.DataFrame(data).T
            data.loc[data.Display_name =="", 'Display_name'] = data.loc[data.Display_name ==""].index.tolist()
            rename_dict = data.Display_name.to_dict()
            domdf.dom = domdf.dom.replace(rename_dict)
            to_display = data.query('Include ==1').Display_name.unique().tolist()
            domdf = domdf.query('dom in @to_display')

    if isinstance(dom, str):
        dom = [dom]
    if dom[0] =='full':
        block_id = domdf.block_id.tolist()
    else:
        if iteration < 1:
            raise ValueError(f"iteration must be >= 1, got {iteration}")
        missing = set(dom) - set(domdf.dom)
        if missing:
            warnings.warn(f"seed domains not found in domdf: {sorted(missing)}")
        for x in range(iteration):
            block_id = domdf.query('dom in @dom').block_id.tolist()
            dom = domdf.query('block_id in @block_id').dom.unique().tolist()
    to_net = domdf.query('block_id in @block_id')[['pid', 'block_id', 'dom', gene_pos, feature_pos]]
    to_net.columns = ['pid', 'block_id', 'dom', 'gpos', 'fpos']
    to_net = to_net.merge(to_net, on ='block_id')
    # Edge type from positions, not row order: genes by gene_pos, domains of the same protein by feature_pos
    same_pid = to_net.pid_x == to_net.pid_y
    to_net['edge_type'] = np.select(
            [~same_pid & (to_net.gpos_x < to_net.gpos_y),
             ~same_pid & (to_net.gpos_x > to_net.gpos_y),
             same_pid & (to_net.fpos_x < to_net.fpos_y),
             same_pid & (to_net.fpos_x > to_net.fpos_y)],
            ['Downstream', 'Upstream', 'Cterminal', 'Nterminal'],
            default='')

    to_net= to_net.query('edge_type != ""')
    to_net = to_net[['block_id','dom_x','dom_y','edge_type']]
    to_net.columns = ['block_id','source','target','edge_type']
    # Count each relation once per block, however many times it occurs there
    to_net = to_net.drop_duplicates()
    if not self_relation:
        to_net = to_net[to_net.source != to_net.target]
    if to_yaml:
        to_net['edge_count'] = to_net.groupby(['source', 'target']).block_id.transform('nunique')
        to_net = to_net.query('edge_count >= @ min_connections').drop(['edge_count'], axis=1)
        if to_net.empty:
            return {}
        yamldf = pd.pivot_table(to_net,
                            index='source',
                            columns='edge_type',
                            values = ['target'],
                            aggfunc={
                                'target':lambda x : count_series(x, as_list=True)}
                            )
        # Edge types absent from the filtered relations get no pivot column
        yamldf = yamldf.loc[:,'target'].reindex(columns=['Cterminal', 'Downstream', 'Nterminal', 'Upstream'], fill_value='')
        yamldf = to_net.groupby('source').block_id.nunique().rename('total').to_frame().join(yamldf).sort_values('total')
        yamldf.fillna('', inplace=True)
        yamldf.Nterminal = yamldf.Nterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Cterminal = yamldf.Cterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Downstream = yamldf.Downstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Upstream = yamldf.Upstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf['L'] = yamldf[['Cterminal', 'Downstream', 'Nterminal', 'Upstream']].apply(
                lambda row: [item for sublist in row for item in sublist if item != ''],
                axis=1
                )
        yamldf.L = yamldf.L.apply(len)
        yamldf = yamldf.query('L >0')
        yamldf = yamldf.drop(['L'], axis=1)

        return(yamldf.to_dict(orient='index'))
    else:
        to_net['edge_total'] = to_net.groupby(['source', 'target']).block_id.transform('nunique')
        to_net = to_net.groupby(['source', 'target', 'edge_type']).agg(
                {'block_id' :'count',
                 'edge_total' : 'first'}).reset_index().sort_values('block_id')
        to_net = to_net.rename({'block_id':'edge_count'},axis=1)
        to_net = to_net.query('edge_total >= @min_connections')
        return(to_net)


def adjacency_edges(df,
                    target='pfam',
                    *,
                    ftype='CDS',
                    interaction=True,
                    ignore=None,
                    separator='+',
                    separator_regex=False,
                    replace=None,
                    replace_regex=False,
                    orient=False,
                    stats=None):
    """
    List the edges between features that are adjacent in a block.

    Each feature is labelled by one of its annotation columns, and each
    label is joined to the next one in the same block. A protein with
    several domains, written ``'A+B'``, contributes one label per domain,
    so its domains are joined to each other (a fusion) and its last domain
    to the first label of the next protein (a neighbour). Identical edges
    are then counted.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per feature, already in genome order within each block (the
        function does not sort). Needs the columns ``block_id``, ``pid``
        and ``type``, and the column or columns named by ``target``. It is
        not modified.
    target : str or list of str, default 'pfam'
        Column holding the label of each feature. With a list, the label is
        the first value that is not missing, going through the columns in
        the order given. A feature with no label in any of them is labelled
        ``'?'``.
    ftype : str or list of str, default 'CDS'
        Values of the ``type`` column to keep. Other rows are dropped
        before the edges are made.
    interaction : bool, default True
        Add an ``interaction`` column that tells fusions from neighbours,
        and count each kind separately. With False, both are counted
        together.
    ignore : iterable of str, optional
        Labels to remove, after ``replace`` is applied and before the
        edges are made. The labels around a removed one become adjacent. A
        single string is one label, not a list of characters.
    separator : str or None, default '+'
        Splits each label into several, one per domain, in order. None
        keeps each label whole.
    separator_regex : bool, default False
        Read ``separator`` as a regular expression.
    replace : mapping, optional
        Renames labels, after splitting. Keys are the old labels and values
        the new ones. None renames nothing.
    replace_regex : bool, default False
        Read the keys of ``replace`` as regular expressions.
    orient : bool, default False
        Make neighbour edges undirected by writing each pair with the
        smaller label first, so that A to B and B to A become one edge.
        Fusion edges are not reordered, since the order of domains in a
        protein is meaningful.
    stats : dict, optional
        Named aggregations passed to ``groupby(...).agg``, as
        ``{name: (column, function)}``, over the columns ``block_id``
        and ``pid`` (the protein of ``source``), and the ``source``,
        ``target`` and ``interaction`` columns.
        Defaults to ``{'weight': ('block_id', 'count'), 'blocks':
        ('block_id', 'nunique')}``.

    Returns
    -------
    pandas.DataFrame
        One row per distinct edge, sorted by ``source``, ``target`` and
        ``interaction``, with the columns:

        - ``source``, ``target``: the two labels, ``target`` being the one
          that follows ``source`` in the block.
        - ``interaction``: ``'fusion'`` when both labels come from the same
          protein, ``'neighbor'`` otherwise. Only with ``interaction``.
        - one column per entry of ``stats``; by default ``weight``, the
          number of rows with this edge, and ``blocks``, the number of
          blocks with it.

        The frame is empty, with these columns, when no edge is found.

    Raises
    ------
    KeyError
        If ``block_id``, ``pid``, ``type`` or a ``target`` column is
        missing from ``df``.

    See Also
    --------
    get_relations : Relates every pair of domains in a block, typed by
        relative position.
    network_annotation : Builds a graph from relations.
    rotifer.genome.data.NeighborhoodDF.to_network : Applies this function to a
        neighborhood table, ordering each block along the query.

    Notes
    -----
    Unlike :func:`get_relations`, which relates every pair of domains found
    in a block and types them as ``Cterminal``, ``Nterminal``,
    ``Upstream`` or ``Downstream``, this function relates only features
    that are next to each other, and types them as ``fusion`` or
    ``neighbor``. Its edges go from a label to the next one, so each pair
    is listed once, not from both ends.

    ``weight`` counts rows of the edge table, so a block in which the edge
    occurs twice adds 2; ``blocks`` counts each block once. ``orient``
    changes neighbour edges only, so it can merge two neighbour edges but
    never a fusion with a neighbour. Labels that are integers stay
    integers in ``source`` and ``target``.

    Two labels are a fusion when the two rows they come from have the
    same ``pid``. Labels of one protein are always consecutive, because
    the split keeps them in place, unless ``ignore`` removes the ones
    between them.

    The edges follow the order of the rows, so sort ``df`` into the order
    in which you want the block read (for example by strand-aware genome
    position) before calling this function.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     'block_id': ['b1', 'b1', 'b1', 'b2', 'b2'],
    ...     'pid': ['p1', 'p2', 'p3', 'p4', 'p5'],
    ...     'type': ['CDS'] * 5,
    ...     'pfam': ['A', 'B+C', 'D', 'A', 'B+C']})
    >>> adjacency_edges(df)
      source target interaction  weight  blocks
    0      A      B    neighbor       2       2
    1      B      C      fusion       2       2
    2      C      D    neighbor       1       1
    """
    cols = [target] if isinstance(target, str) else list(target)
    missing = [c for c in ['block_id', 'pid', 'type'] + cols if c not in df.columns]
    if missing:
        raise KeyError(f"columns missing from df: {missing}")
    if isinstance(ftype, str):
        ftype = [ftype]
    if stats is None:
        stats = {'weight': ('block_id', 'count'), 'blocks': ('block_id', 'nunique')}

    w = df.loc[df['type'].isin(ftype), ['block_id', 'pid'] + cols].reset_index(drop=True)
    label = w[cols[0]]
    for col in cols[1:]:
        label = label.where(label.notna(), w[col])
    w = w[['block_id', 'pid']].assign(source=label.fillna('?'))

    if separator is not None:
        w['source'] = w['source'].astype(str).str.split(separator, regex=separator_regex)
        w = w.explode('source').reset_index(drop=True)
    if replace is not None:
        w['source'] = w['source'].replace(replace, regex=replace_regex)
    if ignore is not None:
        ignored = [ignore] if isinstance(ignore, str) else list(ignore)
        w = w[~w['source'].isin(ignored)].reset_index(drop=True)

    # The target is the next row's label. Built from a numpy array so the
    # label dtype is kept (shift() would upcast integer labels to float).
    # The last row's target is a placeholder: it never has a same-block
    # successor, so it is dropped below.
    src = w['source'].to_numpy()
    nxt = w[['block_id', 'pid']].shift(-1)
    w = w.assign(target=np.concatenate([src[1:], src[-1:]]))
    same_block = (w['block_id'] == nxt['block_id']).to_numpy()
    fusion = same_block & (w['pid'] == nxt['pid']).to_numpy()
    w = w.loc[same_block, ['block_id', 'pid', 'source', 'target']].assign(
        interaction=np.where(fusion[same_block], 'fusion', 'neighbor'))

    if orient:
        swap = ((w['interaction'] == 'neighbor') & (w['source'] > w['target'])).to_numpy()
        w.loc[swap, ['source', 'target']] = w.loc[swap, ['target', 'source']].to_numpy()

    keys = ['source', 'target'] + (['interaction'] if interaction else [])
    edges = w.groupby(keys).agg(**stats).reset_index()
    return edges.sort_values(keys).reset_index(drop=True)
