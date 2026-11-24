#!/usr/bin/env python3
"""Tests for the deprecated rotifer.devel.alpha.net_functions shim (issue #121).

What has to hold, and why:

  * importing the old path raises a FutureWarning (shown by default, unlike
    DeprecationWarning) that names the new modules and the removal date,
    and points at the importing line, not at Python's import machinery
  * the shim holds no code of its own: every name, including the private
    helpers old code may reach for, is the object from rotifer.graph.network
    or rotifer.view.network

The import runs in a fresh interpreter, because Python warns only on the
first import of a module in a session.

Run under pytest: PYTHONPATH=lib pytest test/devel/test_net_functions_shim.py
"""

import os
import subprocess
import sys

import pytest

pytest.importorskip('dash')
pytest.importorskip('dash_cytoscape')
pytest.importorskip('networkx')

LIB = os.path.join(os.path.dirname(__file__), '..', '..', 'lib')


def test_import_warns_at_the_importing_line(tmp_path):
    script = tmp_path / 'old_notebook.py'
    script.write_text('from rotifer.devel.alpha import net_functions\n')
    env = dict(os.environ, PYTHONPATH=os.path.abspath(LIB))
    result = subprocess.run([sys.executable, str(script)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert f'{script}:1: FutureWarning' in result.stderr
    for text in ('rotifer.graph.network', 'rotifer.view.network', '2026-11-30'):
        assert text in result.stderr


def test_shim_reexports_the_new_objects():
    with pytest.warns(FutureWarning):
        sys.modules.pop('rotifer.devel.alpha.net_functions', None)
        from rotifer.devel.alpha import net_functions as old
    from rotifer.graph import network as graph
    from rotifer.view import network as view
    for name in ('get_relations', 'network_annotation', 'calculate_positions', 'get_network_community',
                 'community_stability', 'clique_jaccard_community', 'subgraph_stats', 'count_series',
                 '_simple_graph', '_cpm_resolution', '_undirected_max', '_community_scheme', '_clique_pairs'):
        assert getattr(old, name) is getattr(graph, name)
    assert old.network_dash is view.network_dash
    assert old.convert_to_cytoscape is view.convert_to_cytoscape
