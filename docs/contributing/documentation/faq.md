# FAQ

Short answers to the question that comes up most often: *which file
do I edit to change this?* The site has no single layout file. Seven
layers each own one part of the result, and editing the wrong one
is the usual reason a change does nothing.

## Which file controls what

```{list-table}
:header-rows: 1

* - I want to change
  - Edit
* - The links in the header, the theme, the accent color, which
    source files are ignored
  - `docs/conf.py`
* - Which sections exist, their order, which page nests under which
  - The `toctree` in each `index.md` or `index.rst`
* - Which modules get an API page, and how they are grouped
  - `docs/api/index.rst`
* - The anatomy of a generated API page
  - `docs/_templates/autosummary/module.rst` and `class.rst`
* - How the sidebar nests and labels the API entries
  - `docs/_ext/sidebar_tree.py`
* - Either sidebar's collapse control
  - `docs/_templates/partials/globaltoc-above.html` (left),
    `docs/_templates/sidebars/localtoc.html` (right) and
    `docs/_static/sidebar-toggle.js` (both)
* - Where the built HTML lands
  - `docs/Makefile` locally, `.readthedocs.yaml` when hosted
```

## How is the navigation built?

Two independent mechanisms, and both have to agree.

The header links are a hand-written list in the `nav_links` key of
`html_theme_options` in `docs/conf.py`. They are plain strings, not
cross-references: Sphinx does not check them, so a typo there
produces a dead link that `make strict` will not catch.

Everything else, the sidebar, the previous and next buttons, the
page hierarchy, comes from the `toctree` directives. `docs/index.md`
holds the root one, which lists the top-level sections in the order
they appear. Each section index then holds its own `toctree`
listing its pages. A page reachable from the root toctree, through
however many levels, is part of the site; a page that is not is an
orphan.

Which sidebar widgets appear is set by `html_sidebars` in
`docs/conf.py`. This site deliberately keeps only the local table
of contents, because the other widgets the theme offers call out to
external services.

## Why is the API sidebar a tree when `api/index.rst` is a flat list?

Because `docs/_ext/sidebar_tree.py` rewrites it. autosummary emits
one flat `toctree` entry per documented object, which put
`rotifer.db.ncbi.entrez` beside `rotifer.db` rather than inside it
and repeated the parent name in every label. The extension regroups
the sidebar entries whose titles are dotted Python names into the
tree those names describe and shortens each label to its own last
component, so `rotifer.db.ncbi.entrez` reads as `entrez` inside
`ncbi` inside `db` inside `rotifer`.

Three consequences are worth knowing:

- Every node with children gets a collapse toggle, because the theme
  adds one to any sidebar item that has children. Deepening the tree
  is what produces the toggles.
- Siblings are sorted, packages first and then pages, each group
  alphabetical. The curated order of the `autosummary` blocks in
  `docs/api/index.rst` still decides the order *on the API reference
  page*, but it cannot survive in the sidebar: it is an order over a
  flat list, and the sidebar is a tree.
- A package listed only through its submodules, such as
  `rotifer.core.io`, has no page of its own, so its label only
  toggles the branch instead of linking anywhere. Adding the package
  to `docs/api/index.rst` turns it into a link.

The rewrite is HTML post-processing and touches no source file. Set
`NEST_API_SIDEBAR = False` in `docs/conf.py` to get the flat list
back.

## How do the collapsible sidebars work?

Both sidebars collapse to a 3rem rail, and both work the same way:
a class on `<html>`, a rail in CSS, and `.sy-main` widening to
reclaim the difference. Shibuya computes `.sy-main`'s width from the
sidebars rather than letting it flex, so every combination of the
two is spelled out in the stylesheet.

**The controls.** The left one lives in
`docs/_templates/partials/globaltoc-above.html`. The theme renders
the left sidebar inline in its layout with no sidebar template to
override, and this partial is the hook it does provide: it sits
inside the scroll container, directly above the toctree. It is a
bare `tabler:menu-2` hamburger with no label, because the navigation
has no heading of its own to collapse into.

The right one is the panel's own heading, in
`docs/_templates/sidebars/localtoc.html`, which overrides the
Shibuya template of the same name. The whole "On this page" row is
the button, with a `tabler:layout-sidebar-right-collapse` icon that
flips to the `-expand` variant when collapsed. Every icon is
vendored in `_static/vendor/iconify-preload.js`, like the rest of
the site.

**The state.** `docs/_static/sidebar-toggle.js` puts
`lside-collapsed` or `rside-collapsed` on `<html>` and remembers
each choice in `localStorage`, so both stay as you left them from
page to page. Two details there are load-bearing:

- It binds each click to the button itself rather than delegating
  from the document. Shibuya registers a `stopPropagation` click
  handler on the aside to drive its drawer, so no click inside a
  sidebar ever reaches the document and a delegated handler
  silently never fires.
- It is loaded from `<head>` through
  `docs/_templates/partials/extra-head.html`, not through
  `html_js_files`, because Shibuya renders those at the end of
  `<body>`, which would paint both sidebars expanded for one frame
  first.

**The appearance.** `docs/_static/theme.css` owns the rails, which
icon shows, and the animation. The width eases over
`--rot-transition-panel`, but only while a toggle is in flight: the
script adds a `sidebars-animating` class for the duration, because a
transition left on permanently would also animate every window
resize, which reads as lag rather than as motion.

**The breakpoints.** Each panel collapses only where the theme has
it in flow: 768px for the left, 1280px for the right. Below that the
panel is an overlay drawer with its own close button, so the
stylesheet makes the control inert and the script ignores it, and
the two mechanisms never fight over the same panel.

Right-panel entry labels are shortened by
`docs/_ext/sidebar_tree.py`, the same file that nests the left
sidebar. autodoc writes a method entry as `ClassName.method()`, one
unbreakable token wider than the panel, so the names used to be
clipped mid-word against its edge; the child now drops the prefix it
inherits from its parent and reads `method()`.

## How do I add a new page?

Write the file, then add its name, without the extension, to the
`toctree` of the index of the section it belongs to. Position in
that list is position in the sidebar; the list is not sorted.

If the page is a new top-level section, create a directory with an
`index.md` inside it, add that `index` to the root toctree in
`docs/index.md`, and add a matching entry to `nav_links` in
`docs/conf.py` if it deserves a header link. The landing page grid
cards in `docs/index.md` are a separate, hand-written copy of the
same links, so update them too or the two will disagree.

## Why does Sphinx say my page is not in any toctree?

Because nothing links to it. Every source file under `docs/` is
read, and a file no toctree references is an orphan. Under
`make strict` the warning is an error, so the build fails.

Three ways out: add it to a toctree, which is almost always what you
want; add it to `exclude_patterns` in `docs/conf.py` if it is a
working note that should not be published, which is how
`OPEN_QUESTIONS.md` is handled; or put `:orphan:` at the top of the
file if it is genuinely meant to be reachable only by direct link.

## How do I add a module to the API reference?

Add its dotted name to one of the `autosummary` blocks in
`docs/api/index.rst`, under the heading that fits it. The section
headings in that file are the grouping you see on the page, and the
order inside each block is the order rendered.

Everything under `docs/api/generated/` is written by autosummary at
build time from the templates in `docs/_templates/autosummary/`.
Never edit those files by hand and never commit them: your changes
will be silently overwritten on the next build. To change what a
generated page contains, for every module at once, edit the
template.

If adding a module breaks the build, do not fix the source file as
part of a docs change. See the troubleshooting table in
{doc}`/contributing/documentation/building-the-docs`.

## Where does the built HTML go?

Locally, into `docs/_build/html`, set by `BUILDDIR` in
`docs/Makefile`. The directory is build output and is not committed.

On Read the Docs the build is driven by `.readthedocs.yaml` at the
repository root, which points Sphinx at `docs/conf.py` and sets
`fail_on_warning: true`. That makes the hosted build exactly as
strict as `make strict`, so a clean local strict build means a clean
hosted one.

## Why are some pages Markdown and some reStructuredText?

Both work. Markdown pages are parsed by MyST, which supports every
Sphinx directive through the `{directive}` fence syntax, so nothing
is lost by writing prose in Markdown. reStructuredText is used where
a page is mostly directives, such as `docs/api/index.rst`, where the
reST form is shorter. Pick whichever makes the page you are writing
easier to read.

## Further reading

This site is a fairly ordinary Sphinx project, so the upstream
documentation answers anything not covered above.

- [Sphinx documentation](https://www.sphinx-doc.org/en/master/),
  the starting point for how a build works.
- [The toctree directive](https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-toctree),
  which defines the page hierarchy.
- [HTML configuration](https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output),
  including `html_theme_options` and `html_sidebars`.
- [autosummary](https://www.sphinx-doc.org/en/master/usage/extensions/autosummary.html)
  and [autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html),
  which generate the API pages.
- [MyST parser](https://myst-parser.readthedocs.io/en/latest/),
  for the Markdown syntax used here.
- [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html),
  for the docstring format.
- [Shibuya](https://shibuya.lepture.com/), the theme, for the
  options accepted by `html_theme_options`.
- [Read the Docs configuration](https://docs.readthedocs.io/en/stable/config-file/v2.html),
  for `.readthedocs.yaml`.
