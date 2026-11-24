# Network functions

Domains that occur in the same conserved gene neighbourhood (a block)
become nodes in a network; shared blocks become edges. The analysis
functions live in `rotifer.graph.network` and need the ``graph`` extra
(``pip install rotifer[graph]``). The interactive viewer,
`rotifer.view.network.network_dash`, needs the ``view`` extra
(``pip install rotifer[view]``). Nothing contacts a remote service.

## Quickstart

Every snippet below runs offline from the root of a repository checkout
(one reads `test/data/graph/toy_domains.tsv`). The input
table has columns ``block_id``, ``pid``, ``start``, ``end``, ``strand``,
``dom``, ``dom_start``.

::::{tab-set}

:::{tab-item} Relations

```python
import pandas as pd
from rotifer.graph import network as gn

# In-memory table like the two_blocks fixture
rows = [
    ('b1', 'p1', 100, 200, 1, 'A', 1),
    ('b1', 'p1', 100, 200, 1, 'B', 50),
    ('b1', 'p2', 300, 400, 1, 'C', 1),
    ('b2', 'q2', 300, 400, 1, 'C', 1),
    ('b2', 'q1', 100, 200, 1, 'B', 50),
    ('b2', 'q1', 100, 200, 1, 'A', 1),
]
df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'end', 'strand', 'dom', 'dom_start'])

edges = gn.get_relations(df, min_connections=1, gene_pos="start", feature_pos="dom_start")
print(edges[['source', 'target', 'edge_type', 'edge_count', 'edge_total']].to_string(index=False))
#   source target  edge_type  edge_count  edge_total
#        A      B  Cterminal           2           2
#        A      C Downstream           2           2
#        B      A  Nterminal           2           2
#        B      C Downstream           2           2
#        C      A   Upstream           2           2
#        C      B   Upstream           2           2
```
:::

:::{tab-item} Graph and layouts

```python
import pandas as pd
from rotifer.graph import network as gn

rows = [
    ('b1', 'p1', 100, 200, 1, 'A', 1),
    ('b1', 'p1', 100, 200, 1, 'B', 50),
    ('b1', 'p2', 300, 400, 1, 'C', 1),
    ('b2', 'q2', 300, 400, 1, 'C', 1),
    ('b2', 'q1', 100, 200, 1, 'B', 50),
    ('b2', 'q1', 100, 200, 1, 'A', 1),
]
df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'end', 'strand', 'dom', 'dom_start'])

edges = gn.get_relations(df, min_connections=1, gene_pos="start", feature_pos="dom_start")
G = gn.network_annotation(edges)
layouts = gn.calculate_positions(G, iterations=10)
print(list(layouts))
# ['Spring Layout', 'Kamada-Kawai from Spring layout', 'Kamada-Kawai', 'Circular']
print(sorted(layouts['Spring Layout']))
# ['A', 'B', 'C']
```
:::

:::{tab-item} Communities

```python
import pandas as pd
from rotifer.graph import network as gn

rows = [
    ('b1', 'p1', 100, 200, 1, 'A', 1),
    ('b1', 'p1', 100, 200, 1, 'B', 50),
    ('b1', 'p2', 300, 400, 1, 'C', 1),
    ('b2', 'q2', 300, 400, 1, 'C', 1),
    ('b2', 'q1', 100, 200, 1, 'B', 50),
    ('b2', 'q1', 100, 200, 1, 'A', 1),
]
df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'end', 'strand', 'dom', 'dom_start'])

edges = gn.get_relations(df, min_connections=1, gene_pos="start", feature_pos="dom_start")
G = gn.network_annotation(edges)
schemes = gn.get_network_community(G, cpm_resolution=0.5)
for name, (node_to_label, label_to_colour) in schemes.items():
    print(f"{name}: {node_to_label}")
# Louvain: {'A': 'Louvain_0', 'B': 'Louvain_0', 'C': 'Louvain_0'}
# Leidein: {'A': 'Leiden_0', 'B': 'Leiden_0', 'C': 'Leiden_0'}
# Leidein_2: {'A': 'LeidenCPM_0', 'B': 'LeidenCPM_0', 'C': 'LeidenCPM_0'}
```
:::

:::{tab-item} Domain counts

```python
import pandas as pd
from rotifer.graph import network as gn

s = pd.Series(["SP+TM", "TM", "SP+TM", None])
print("Domain counts:", gn.count_series(s))
print("Architecture counts:", gn.count_series(s, count="architecture", as_list=True))
print("Normalized, cutoff 0.5:", gn.count_series(s, normalize=True, cut_off=0.5))
# Domain counts: TM(3), SP(2)
# Architecture counts: ['SP+TM(2)', 'TM(1)']
# Normalized, cutoff 0.5: TM(60.00%)
```
:::
::::

## Functions

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Function
  - Module
  - Input
  - Output
  - Needs
* - {py:func}`rotifer.graph.network.get_relations`
  - rotifer.graph.network
  - DataFrame with block_id, pid, dom, plus gene_pos and feature_pos columns
  - DataFrame: source, target, edge_type, edge_count, edge_total
  - graph extra
* - {py:func}`rotifer.graph.network.adjacency_edges`
  - rotifer.graph.network
  - DataFrame with block_id, pid, type and a label column (pfam by default), in genome order
  - DataFrame: source, target, interaction (fusion or neighbor), weight, blocks
  - pandas only
* - {py:func}`rotifer.graph.network.network_annotation`
  - rotifer.graph.network
  - Relations DataFrame from get_relations
  - networkx.Graph (or MultiGraph/MultiDiGraph) with styled nodes and edges
  - graph extra
* - {py:func}`rotifer.graph.network.calculate_positions`
  - rotifer.graph.network
  - networkx.Graph, iterations (int)
  - Dict of layout name to {node: array([x, y])}
  - graph extra
* - {py:func}`rotifer.graph.network.get_network_community`
  - rotifer.graph.network
  - networkx.Graph, weight, cpm_resolution, seed
  - Dict: {"Louvain"|"Leidein"|"Leidein_2": [node_to_label, label_to_colour]}
  - graph extra (igraph, leidenalg, seaborn)
* - {py:func}`rotifer.graph.network.community_stability`
  - rotifer.graph.network
  - networkx.Graph, weight, cpm_resolution, n_seeds, seed
  - (summary DataFrame, nodes DataFrame)
  - graph extra
* - {py:func}`rotifer.graph.network.clique_jaccard_community`
  - rotifer.graph.network
  - networkx.Graph, min_jaccard_index, cpm_resolution, min_clique_size, seed
  - DataFrame of connected clique pairs with communities
  - graph extra
* - {py:func}`rotifer.graph.network.subgraph_stats`
  - rotifer.graph.network
  - net, subnet (networkx.Graph), theme_dict, theme_to_stats
  - SignificanceResult (or (result, table) with return_df=True)
  - graph extra (scipy)
* - {py:func}`rotifer.graph.network.count_series`
  - rotifer.graph.network
  - pandas.Series of architectures (domain names joined by +)
  - String or list of "name(count)" or "name(%)"
  - pandas only
* - {py:func}`rotifer.view.network.convert_to_cytoscape`
  - rotifer.view.network
  - Graph, positions dict, layout_key, scaling_factor, x_scale, y_scale
  - List of Cytoscape node and edge elements
  - view extra (dash, dash-cytoscape)
* - {py:func}`rotifer.view.network.network_dash`
  - rotifer.view.network
  - Graph, pos_dict, xx (group->colour), function (node->group), color_schemes, and keyword arguments for Dash.run such as port
  - None, once the local Dash server stops
  - view extra
```

## Common patterns

Check community stability before reporting a partition:

```python
import pandas as pd
from rotifer.graph import network as gn

rows = [
    ('b1', 'p1', 100, 200, 1, 'A', 1),
    ('b1', 'p1', 100, 200, 1, 'B', 50),
    ('b1', 'p2', 300, 400, 1, 'C', 1),
    ('b2', 'q2', 300, 400, 1, 'C', 1),
    ('b2', 'q1', 100, 200, 1, 'B', 50),
    ('b2', 'q1', 100, 200, 1, 'A', 1),
]
df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'end', 'strand', 'dom', 'dom_start'])

edges = gn.get_relations(df, min_connections=1, gene_pos="start", feature_pos="dom_start")
G = gn.network_annotation(edges)
summary, nodes = gn.community_stability(G, weight='edge_total', n_seeds=5)
print(summary)
#            communities  distinct_partitions  ari_median  ari_min
# partition
# Louvain              1                    1         1.0      1.0
# Leidein              1                    1         1.0      1.0
# Leidein_2            1                    1         1.0      1.0
```

Test whether a theme is enriched in a community with `subgraph_stats`:

```python
import pandas as pd
from rotifer.graph import network as gn

domdf = pd.read_csv("test/data/graph/toy_domains.tsv", sep="\t")
edges = gn.get_relations(domdf, min_connections=2, gene_pos="start", feature_pos="dom_start")
G = gn.network_annotation(edges)
labels = gn.get_network_community(G)['Louvain'][0]
subnet = G.subgraph([n for n, lab in labels.items() if lab == 'Louvain_0'])
theme = {
    **dict.fromkeys(["SBP_bac_5", "BPD_transp_1", "ABC_tran", "oligo_HPY"], "transport"),
    **dict.fromkeys(["Response_reg", "Trans_reg_C", "HAMP", "HisKA", "HATPase_c"], "signalling"),
    **dict.fromkeys(["PhdYeFM_antitox", "ParE_toxin", "MethyltransfD12", "RE_AlwI"], "defence"),
}
result, table = gn.subgraph_stats(G, subnet, theme, "transport", return_df=True)
print(f"odds ratio {float(result.statistic):.2f}, p {float(result.pvalue):.3f}")
print(table)
# odds ratio 0.00, p 0.105
#               Theme  Not theme  Unannotated
# Subnetwork        0          5            0
# net - subnet      4          4            0
```

Per-cluster domain composition with `groupby(...).agg(count_series)`:

```python
import pandas as pd
from rotifer.graph import network as gn

rows = [
    ('b1', 'p1', 100, 200, 1, 'A', 1),
    ('b1', 'p1', 100, 200, 1, 'B', 50),
    ('b1', 'p2', 300, 400, 1, 'C', 1),
    ('b2', 'q2', 300, 400, 1, 'C', 1),
    ('b2', 'q1', 100, 200, 1, 'B', 50),
    ('b2', 'q1', 100, 200, 1, 'A', 1),
]
df = pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'end', 'strand', 'dom', 'dom_start'])

edges = gn.get_relations(df, min_connections=1, gene_pos="start", feature_pos="dom_start")
G = gn.network_annotation(edges)
schemes = gn.get_network_community(G)
labels = schemes['Louvain'][0]
df_nodes = pd.DataFrame({'node': list(labels.keys()), 'community': list(labels.values())})
df_nodes['architecture'] = df_nodes['node']  # each node is one domain here
comp = df_nodes.groupby('community')['architecture'].agg(gn.count_series)
print(comp)
# community
# Louvain_0    A(1), B(1), C(1)
# Name: architecture, dtype: str
```

`get_relations` relates every pair of domains in a block. To join only
features that sit next to each other, use `adjacency_edges`. It splits
multidomain labels on `+`, so domains of one protein become fusion edges and
the last domain of a protein joins the first label of the next one:

```python
import pandas as pd
from rotifer.graph import network as gn

df = pd.DataFrame({
    "block_id": ["b1", "b1", "b1", "b2", "b2"],
    "pid": ["p1", "p2", "p3", "p4", "p5"],
    "type": ["CDS"] * 5,
    "pfam": ["A", "B+C", "D", "A", "B+C"],
})
print(gn.adjacency_edges(df))
#   source target interaction  weight  blocks
# 0      A      B    neighbor       2       2
# 1      B      C      fusion       2       2
# 2      C      D    neighbor       1       1
```

The function does not sort, so put the rows in genome order first.

## Migrating

The old module `rotifer.devel.alpha.net_functions` is deprecated and will
be removed on 2026-11-30 (issue #153). Importing it emits a
`FutureWarning`:

```
FutureWarning: rotifer.devel.alpha.net_functions is deprecated and will be removed on 2026-11-30 (issue #153); import from rotifer.graph.network (analysis) and rotifer.view.network (network_dash) instead
```

Replace imports as follows:

- Analysis functions: `from rotifer.graph.network import get_relations, network_annotation, calculate_positions, get_network_community, community_stability, clique_jaccard_community, subgraph_stats, count_series`
- Viewer: `from rotifer.view.network import convert_to_cytoscape, network_dash`

The parameter `resolution_parameter` in `get_network_community` and
`clique_jaccard_community` is a deprecated alias of `cpm_resolution`.
Use `cpm_resolution` instead; the old name still works but raises a
`DeprecationWarning`.

## Viewer

The interactive viewer `rotifer.view.network.network_dash` opens a Dash
app with a Cytoscape canvas. It takes a graph from
`network_annotation`, layouts from `calculate_positions` (converted to
lists), a group-to-colour map `xx`, a node-to-group map `function`, and
optional `color_schemes` (for example the output of
`get_network_community`). In the viewer you can switch layouts, filter
by minimum connections, recolour nodes and edges, add or remove edges,
save positions as new layouts, and export the picture. Outside Jupyter
the call blocks until you stop the server.

A working demo is in `test/data/graph/run_dash_demo.py`:

```bash
python test/data/graph/run_dash_demo.py
```

Then browse to `http://127.0.0.1:8050`. Use `--split` for one edge per
edge type.

The call itself (not run here because it starts a server):

```python
# Starts a local web server; not run here.
from rotifer.view.network import network_dash

network_dash(G, pos_dict, xx, function, color_schemes)
```