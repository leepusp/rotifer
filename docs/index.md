---
layout: rotifer-landing
---

:::{grid-item-card} {iconify}`tabler:terminal-2` Command line tools
:link: cli/index
:link-type: doc

Shell commands for HMMER, BLAST and Phobius reports, domain
architectures and database mirrors, each backed by a library function.
:::

:::{grid-item-card} {iconify}`tabler:book-2` API reference
:link: api/index
:link-type: doc

The generated reference for every documented module in the package.
:::

:::{grid-item-card} {iconify}`tabler:tools` Contributing
:link: contributing/index
:link-type: doc

How to write docstrings, build this site, keep the design
consistent, and find which file controls what.
:::

:::{grid-item-card} {iconify}`tabler:git-merge` Changelog
:link: changelog
:link-type: doc

What changed in each release, and what to do when upgrading.
:::

:::{grid-item-card} {iconify}`tabler:license` License
:link: licensing/index
:link-type: doc
# ROTIFER

<!-- The landing page. Its sections come from _templates/landing/ and its
     copy from _data/landing.en.yml; this file only selects the layout and
     holds the root of the site's navigation. -->

```{toctree}
:hidden:
:maxdepth: 2

Documentation <overview>
install
db/index
pandas/blocks
cli/index
api/index
contributing/index
changelog
licensing/index
```
