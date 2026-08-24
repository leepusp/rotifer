# ROTIFER

ROTIFER (Rapid Open-source Tools and Infrastructure For data Exploration and
Research) is a multi-language collection of high-level libraries for building
data analysis pipelines in comparative genomics and computational analysis of
biological sequences, plus command line tools built on that framework.

This site documents the Python package that lives under `lib/rotifer` in the
[repository](https://github.com/leepbioinfo/rotifer).

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {iconify}`tabler:database` Data access
:link: db/index
:link-type: doc

Cursors for biological databases: NCBI, UniProt and local mirrors.
The reference example of how ROTIFER modules are documented.
:::

:::{grid-item-card} {iconify}`tabler:book-2` API reference
:link: api/index
:link-type: doc

The generated reference for every documented module in the package.
:::

:::{grid-item-card} {iconify}`tabler:tools` Contributing
:link: contributing/index
:link-type: doc

How to write docstrings, build this site, and keep the design consistent.
:::
::::

```{toctree}
:hidden:
:maxdepth: 2

db/index
api/index
contributing/index
```
