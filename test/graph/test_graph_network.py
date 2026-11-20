#!/usr/bin/env python3
"""Tests for rotifer.graph.network (issue #121).

What has to hold, and why:

  * get_relations reads the edge type from positions, not row order: genes
    by gene_pos, domains of one protein by feature_pos. Each relation counts
    once per block, edge_count is the number of blocks with that relation
    and edge_total the number of blocks with the pair in any type, which is
    what min_connections filters on. The input table is left unmodified
  * the YAML export returns {} when nothing passes the filter and honours
    self_relation
  * network_annotation keeps get_relations' block counts (it used to sum
    edge_count over types, counting a block once per type) and, in
    undirected mode, merges A->B Cterminal with its mirror B->A Nterminal.
    split_types=True keeps one edge per type; type_from names the node the
    type is read from
  * get_network_community is reproducible for a given seed, numbers
    communities by size and gives every community its own colour (the
    pastel palette has 10 colours and used to repeat); directed graphs and
    empty graphs no longer crash
  * cpm_resolution only changes the CPM partition; resolution_parameter is
    a deprecated alias for it
  * community_stability reports agreement across seeds against the seed's
    own partition
  * clique_jaccard_community pairs each clique once and never with itself
    (self-loops split a chain of overlapping cliques into one community per
    clique), and rejects min_jaccard_index <= 0
  * subgraph_stats leaves nodes missing from theme_dict out of the Fisher
    test, says so, and with return_df=True returns the result with a
    correctly labelled table

The data are tiny synthetic tables plus the toy table in
test/data/graph/toy_domains.tsv; nothing touches the network.

Run under pytest: PYTHONPATH=lib pytest test/graph/test_graph_network.py
"""

import os
import warnings

import pytest

nx = pytest.importorskip('networkx')
pytest.importorskip('scipy')
pytest.importorskip('igraph')
pytest.importorskip('leidenalg')
pytest.importorskip('seaborn')

import pandas as pd
from scipy.stats import fisher_exact

from rotifer.graph import network as gn

TOY = os.path.join(os.path.dirname(__file__), '..', 'data', 'graph', 'toy_domains.tsv')


def relations(df, **kwargs):
    return gn.get_relations(df, gene_pos='start', feature_pos='dom_start', **kwargs)


@pytest.fixture
def toy():
    return pd.read_csv(TOY, sep='\t')


@pytest.fixture
def two_blocks():
    # Block 1: A-B fused (A N-terminal), C in the next gene; block 2: the same, rows shuffled
    rows = [
        ('b1', 'p1', 100, 'A', 1), ('b1', 'p1', 100, 'B', 50), ('b1', 'p2', 300, 'C', 1),
        ('b2', 'q2', 300, 'C', 1), ('b2', 'q1', 100, 'B', 50), ('b2', 'q1', 100, 'A', 1),
    ]
    return pd.DataFrame(rows, columns=['block_id', 'pid', 'start', 'dom', 'dom_start'])


# get_relations

def test_edge_types_come_from_positions_not_row_order(two_blocks):
    edges = relations(two_blocks, min_connections=1)
    types = {(r.source, r.target): r.edge_type for r in edges.itertuples()}
    assert types[('A', 'B')] == 'Cterminal'
    assert types[('B', 'A')] == 'Nterminal'
    assert types[('A', 'C')] == 'Downstream'
    assert types[('C', 'A')] == 'Upstream'


def test_counts_are_blocks(two_blocks):
    doubled = pd.concat([two_blocks, two_blocks])      # every hit twice
    edges = relations(doubled, min_connections=1).set_index(['source', 'target', 'edge_type'])
    assert edges.loc[('A', 'B', 'Cterminal'), 'edge_count'] == 2
    assert edges.loc[('A', 'B', 'Cterminal'), 'edge_total'] == 2


def test_edge_total_counts_blocks_across_types(toy):
    edges = relations(toy, min_connections=2)
    pair = edges[(edges.source == 'ABC_tran') & (edges.target == 'oligo_HPY')]
    assert pair.edge_type.nunique() == 3
    assert pair.edge_count.sum() == 30 and set(pair.edge_total) == {10}


def test_min_connections_filters_on_edge_total(toy):
    assert (relations(toy, min_connections=10).edge_total >= 10).all()


def test_input_is_not_modified(toy):
    before = toy.copy()
    relations(toy, min_connections=2)
    pd.testing.assert_frame_equal(toy, before)


def test_iteration_must_be_positive_and_missing_seeds_warn(toy):
    with pytest.raises(ValueError):
        relations(toy, dom=['HisKA'], iteration=0)
    with pytest.warns(UserWarning, match='NotADomain'):
        relations(toy, dom=['HisKA', 'NotADomain'])


def test_yaml_empty_and_self_relation(toy):
    assert relations(toy, to_yaml=True, min_connections=1000) == {}
    with_self = relations(toy, to_yaml=True, min_connections=2, self_relation=True)
    without = relations(toy, to_yaml=True, min_connections=2)
    assert 'BPD_transp_1' in with_self['BPD_transp_1']['Downstream'][0] or \
        any('BPD_transp_1' in item for item in with_self['BPD_transp_1']['Downstream'])
    assert not any('BPD_transp_1(' in item for item in without['BPD_transp_1']['Downstream'])


# network_annotation

def test_annotation_keeps_block_counts_and_merges_mirrors(toy):
    G = gn.network_annotation(relations(toy, min_connections=2))
    assert type(G) is nx.Graph
    assert G['ABC_tran']['oligo_HPY']['edge_total'] == 10
    assert G['ABC_tran']['oligo_HPY']['type_from'] == 'ABC_tran'


def test_split_types(toy):
    edges = relations(toy, min_connections=2)
    G = gn.network_annotation(edges, split_types=True)
    assert type(G) is nx.MultiGraph
    assert G.number_of_edges('ABC_tran', 'oligo_HPY') == 3
    assert all(key == data['edge_type'] for *_, key, data in G.edges(keys=True, data=True))
    simple = gn.network_annotation(edges)
    assert dict(G.degree()) != {} and nx.get_node_attributes(G, 'degree') == nx.get_node_attributes(simple, 'degree')


def test_mirror_type_flips_when_swapped():
    edges = pd.DataFrame([['C', 'B', 'Cterminal', 1, 1]],
                         columns=['source', 'target', 'edge_type', 'edge_count', 'edge_total'])
    G = gn.network_annotation(edges, split_types=True)
    assert [(d['type_from'], d['edge_type']) for *_, d in G.edges(data=True)] == [('B', 'Nterminal')]


def test_directed_mode_and_missing_edge_total(toy):
    edges = relations(toy, min_connections=2)
    assert type(gn.network_annotation(edges, unidirected_graph=False)) is nx.MultiDiGraph
    with pytest.warns(UserWarning, match='edge_total'):
        gn.network_annotation(edges.drop(columns='edge_total'))


# get_network_community

def test_communities_are_reproducible_and_numbered_by_size():
    G = nx.les_miserables_graph()
    first = gn.get_network_community(G, weight='weight')
    assert first == gn.get_network_community(G, weight='weight')
    for name, (labels, colours) in first.items():
        sizes = pd.Series(labels).value_counts()
        prefix = sizes.index[0].rsplit('_', 1)[0] + '_'
        assert sizes.index[0] == prefix + '0'
        assert set(labels.values()) == set(colours)


def test_every_community_has_its_own_colour():
    labels, colours = gn.get_network_community(nx.les_miserables_graph(), weight='weight', cpm_resolution=1)['Leidein_2']
    assert len(colours) > 10 and len(set(colours.values())) == len(colours)


def test_directed_and_empty_graphs(toy):
    G = gn.network_annotation(relations(toy, min_connections=2), unidirected_graph=False)
    with pytest.warns(UserWarning, match='directed'):
        result = gn.get_network_community(G)
    assert set(result['Louvain'][0]) == set(G)
    with pytest.warns(UserWarning, match='no nodes'):
        assert gn.get_network_community(nx.Graph()) == {k: [{}, {}] for k in ['Louvain', 'Leidein', 'Leidein_2']}


def test_cpm_resolution_only_changes_cpm_and_old_name_warns():
    G = nx.les_miserables_graph()
    low = gn.get_network_community(G, weight='weight', cpm_resolution=0.05)
    high = gn.get_network_community(G, weight='weight', cpm_resolution=5)
    assert low['Louvain'] == high['Louvain'] and low['Leidein'] == high['Leidein']
    assert len(set(high['Leidein_2'][0].values())) > len(set(low['Leidein_2'][0].values()))
    with pytest.warns(DeprecationWarning, match='cpm_resolution'):
        old = gn.get_network_community(G, weight='weight', resolution_parameter=5)
    assert old == high


# community_stability

def test_stability_on_clear_structure():
    G = nx.planted_partition_graph(4, 15, 0.9, 0.01, seed=2)
    with warnings.catch_warnings():
        warnings.simplefilter('error')                      # no low-agreement warning
        summary, nodes = gn.community_stability(G, weight=None, n_seeds=10)
    assert summary.loc['Louvain', 'ari_min'] == 1 and summary.loc['Louvain', 'distinct_partitions'] == 1
    assert (nodes['Louvain_stability'] == 1).all()
    assert nodes['Louvain'].to_dict() == gn.get_network_community(G, weight=None)['Louvain'][0]


def test_stability_needs_two_seeds():
    with pytest.raises(ValueError):
        gn.community_stability(nx.les_miserables_graph(), n_seeds=1)


# clique_jaccard_community

def chain_of_cliques(n=6, size=6, step=4):
    G = nx.Graph()
    for i in range(n):
        nodes = range(i * step, i * step + size)
        G.add_edges_from((a, b) for a in nodes for b in nodes if a < b)
    return G


def test_cliques_pair_once_never_with_themselves():
    table = gn.clique_jaccard_community(chain_of_cliques(), 0.1, 0.5, 6)
    assert (table.source < table.target).all()
    assert len(table) == len(set(zip(table.source, table.target))) == 5


def test_chain_of_cliques_is_one_community_per_overlap_group():
    table = gn.clique_jaccard_community(chain_of_cliques(), 0.1, 0.5, 6)
    assert pd.concat([table.source_louvain, table.target_louvain]).nunique() == 3


def test_clique_threshold_and_empty_input():
    with pytest.raises(ValueError):
        gn.clique_jaccard_community(chain_of_cliques(), 0, 0.5, 6)
    empty = gn.clique_jaccard_community(nx.path_graph(5), 0.2, 0.5, 6)
    assert empty.empty and 'Jaccard_index' in empty.columns


# subgraph_stats

@pytest.fixture
def themed():
    net = nx.path_graph(13)
    sub = net.subgraph(range(9))
    theme = {n: ('TCS' if n < 5 else 'other') for n in range(9)}
    theme.update({n: 'other' for n in range(9, 13)})
    return net, sub, theme


def test_full_annotation_is_plain_fisher(themed):
    net, sub, theme = themed
    assert tuple(gn.subgraph_stats(net, sub, theme, 'TCS')) == tuple(fisher_exact([[5, 4], [0, 4]]))


def test_unannotated_nodes_are_left_out(themed):
    net, sub, theme = themed
    for n in (3, 4):
        del theme[n]
    with pytest.warns(UserWarning, match='2 of 9 subnetwork nodes'):
        result = gn.subgraph_stats(net, sub, theme, 'TCS')
    assert tuple(result) == tuple(fisher_exact([[3, 4], [0, 4]]))


def test_return_df_returns_result_and_labelled_table(themed, capsys):
    net, sub, theme = themed
    result, table = gn.subgraph_stats(net, sub, theme, 'TCS', return_df=True)
    assert capsys.readouterr().out == ''
    assert list(table.columns) == ['Theme', 'Not theme', 'Unannotated']
    assert table.loc['Subnetwork'].tolist() == [5, 4, 0]
    assert tuple(result) == tuple(gn.subgraph_stats(net, sub, theme, 'TCS'))


# count_series (moved from rotifer.devel.alpha.gian_func)

def test_count_series():
    s = pd.Series(['SP+TM', 'TM', None])
    assert gn.count_series(s) == 'TM(2), SP(1)'
    assert gn.count_series(s, normalize=True, cut_off=0.5) == 'TM(66.67%)'
    assert gn.count_series(pd.Series([], dtype=object), as_list=True) == ['']
