#!/usr/bin/env python3
"""Tests for rotifer.view.network, the network_dash viewer (issue #121).

What has to hold, and why:

  * network_dash builds a Dash app with its callbacks registered on that
    app; it used to register some on Dash's global registry, so a second
    viewer in the same session broke the first
  * the redraw callback returns every edge of the graph, also for directed
    graphs (multigraph edges need their key; the filter used to fail with
    HTTP 500)
  * convert_to_cytoscape places each node from the chosen layout, scaled,
    with the layout's two coordinates swapped (x from the second)

The server is never started: Dash.run is replaced by a stub, and callbacks
are called through Flask's test client.

Run under pytest: PYTHONPATH=lib pytest test/view/test_view_network.py
"""

import os

import pytest

pytest.importorskip('dash')
pytest.importorskip('dash_cytoscape')
nx = pytest.importorskip('networkx')

import dash
import pandas as pd
from networkx.readwrite import json_graph

from rotifer.graph import network as gn
from rotifer.view import network as vn

TOY = os.path.join(os.path.dirname(__file__), '..', 'data', 'graph', 'toy_domains.tsv')
GROUP = {
    **dict.fromkeys(['SBP_bac_5', 'BPD_transp_1', 'ABC_tran', 'oligo_HPY'], 'transport'),
    **dict.fromkeys(['Response_reg', 'Trans_reg_C', 'HAMP', 'HisKA', 'HATPase_c'], 'signalling'),
    **dict.fromkeys(['PhdYeFM_antitox', 'ParE_toxin', 'MethyltransfD12', 'RE_AlwI'], 'defence'),
}
COLOUR = {'transport': '#4C9BE8', 'signalling': '#E8A84C', 'defence': '#D9534F'}


def viewer(monkeypatch, directed=False):
    edges = gn.get_relations(pd.read_csv(TOY, sep='\t'), min_connections=2, gene_pos='start', feature_pos='dom_start')
    G = gn.network_annotation(edges, node_color={n: COLOUR[g] for n, g in GROUP.items()}, unidirected_graph=not directed)
    pos = {name: {n: [float(c) for c in p] for n, p in layout.items()}
           for name, layout in gn.calculate_positions(G, iterations=5).items()}
    apps = []
    monkeypatch.setattr(dash.Dash, 'run', lambda self, *a, **k: apps.append(self))
    vn.network_dash(G, pos, COLOUR, GROUP, {'Theme': [GROUP, COLOUR]})
    return G, pos, apps[0]


def redraw(app, G, pos):
    key = next(k for k in app.callback_map if 'cytoscape.elements' in k and 'last-edge-click' in k)
    outputs = [{'id': p.rsplit('.', 1)[0], 'property': p.rsplit('.', 1)[1]} for p in key[2:-2].split('...')]
    def items(lst):
        return [{'id': i, 'property': p, 'value': v} for i, p, v in lst]
    inputs = [('cytoscape', 'tapEdgeData', None), ('cytoscape', 'tapNodeData', None),
              ('layout-selector', 'value', list(pos)[0]), ('group_selection', 'value', list(COLOUR)),
              ('connections-input', 'value', 1), ('Node_scale', 'value', 10), ('x_scale', 'value', 1),
              ('y_scale', 'value', 1), ('scale-size-input', 'value', 900), ('pos-dict-store', 'data', pos),
              ('edge-modification-trigger', 'data', 1)]
    state = [('cytoscape', 'elements', []), ('last-edge-click', 'data', {'id': None, 'time': 0}),
             ('edge-color-cycle', 'data', []), ('edge-style-overrides', 'data', {}),
             ('last-node-click', 'data', {'id': None, 'time': 0}), ('node-color-cycle', 'data', []),
             ('xx-store', 'data', COLOUR), ('function-store', 'data', GROUP),
             ('pos-dict-store', 'data', pos), ('graph-store', 'data', json_graph.node_link_data(G))]
    body = {'output': key, 'outputs': outputs, 'inputs': items(inputs), 'state': items(state),
            'changedPropIds': ['edge-modification-trigger.data']}
    response = app.server.test_client().post('/_dash-update-component', json=body)
    assert response.status_code == 200
    return response.get_json()['response']['cytoscape']['elements']


def test_callbacks_belong_to_each_app(monkeypatch):
    _, _, first = viewer(monkeypatch)
    _, _, second = viewer(monkeypatch)
    assert first.callback_map and set(first.callback_map) == set(second.callback_map)


@pytest.mark.parametrize('directed', [False, True])
def test_redraw_returns_every_edge(monkeypatch, directed):
    G, pos, app = viewer(monkeypatch, directed)
    elements = redraw(app, G, pos)
    assert sum('source' in e['data'] for e in elements) == G.number_of_edges()


def test_convert_to_cytoscape_positions():
    G = nx.path_graph(['a', 'b'])
    elements = vn.convert_to_cytoscape(G, {'L': {'a': [0.1, 0.2], 'b': [0.0, 0.0]}}, 'L', scaling_factor=100, x_scale=2)
    node = next(e for e in elements if e['data'].get('id') == 'a')
    assert node['position'] == {'x': int(0.2 * 100 * 2), 'y': int(0.1 * 100)}
