# Graphs

The `rotifer.graph` package builds and analyses networks of protein
domains that occur together in conserved gene neighbourhoods. Two
domains are related when they appear in the same block: fused in one
protein or in different genes of the block. The analysis functions are
in `rotifer.graph.network`; the interactive viewer is in
`rotifer.view.network`. Both need the ``graph`` extra
(``pip install rotifer[graph]``); the viewer also needs the ``view``
extra (``pip install rotifer[view]``). Nothing here contacts a remote
service.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:graph` rotifer.graph
:link: /api/generated/rotifer.graph
:link-type: doc

Top level package. Imports nothing; import `rotifer.graph.network`
directly.
:::

:::{grid-item-card} {iconify}`tabler:network` rotifer.graph.network
:link: /api/generated/rotifer.graph.network
:link-type: doc

Relations from domain tables, styled networkx graphs, layouts,
communities, enrichment tests and domain counts.
:::

:::{grid-item-card} {iconify}`tabler:eye` rotifer.view.network
:link: /api/generated/rotifer.view.network
:link-type: doc

Interactive Dash Cytoscape viewer for domain networks.
:::

:::{grid-item-card} {iconify}`tabler:map` Network functions
:link: network
:link-type: doc

Build a network from a domain table, find its communities, test
enrichment and open the viewer.
:::
::::

```{toctree}
:hidden:

network
```