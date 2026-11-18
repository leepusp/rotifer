"""
Deprecated location of :mod:`rotifer.graph.network` and :mod:`rotifer.view.network`.

.. deprecated::
   ``rotifer.devel.alpha.net_functions`` moved to :mod:`rotifer.graph.network`
   (analysis) and :mod:`rotifer.view.network` (the ``network_dash`` viewer).
   This module only re-exports them and will be removed on 2026-11-30
   (issue #153).
"""
import warnings

# stacklevel=2 names the importing line: warnings skips the import machinery's frames
warnings.warn(
    "rotifer.devel.alpha.net_functions is deprecated and will be removed on 2026-11-30 "
    "(issue #153); import from rotifer.graph.network (analysis) and "
    "rotifer.view.network (network_dash) instead",
    FutureWarning,
    stacklevel=2,
)

from rotifer.graph.network import *  # noqa: E402,F401,F403
from rotifer.graph.network import (  # noqa: E402,F401
    _clique_pairs,
    _community_scheme,
    _cpm_resolution,
    _simple_graph,
    _undirected_max,
)
from rotifer.view.network import convert_to_cytoscape, network_dash  # noqa: E402,F401
