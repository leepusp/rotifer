# Design system

## Where the tokens live

Every color, font stack, corner radius and transition used by this
site is defined in exactly one place: the two token blocks in
`docs/_static/css/tokens.css`. The first `:root` block holds the
light palette; the `html.dark` block directly beneath it holds the
dark palette (Shibuya switches modes by setting the `dark` class on
the root element). Every other stylesheet only references the
tokens; no other CSS file, template or directive may contain a
color literal.

The other files in `docs/_static/css/` each have one job, and
`html_css_files` in `docs/conf.py` loads them in this order:
`fonts.css` (the self-hosted typefaces), `tokens.css`,
`theme-map.css` (tokens mapped onto the Shibuya and sphinx-design
variables), `content.css` (prose, code, cards, tables),
`navigation.css` (header and sidebars) and `footer.css`. The landing
page adds `docs/_static/landing/landing.css`.

To change the accent, edit the `--rot-accent-*` tokens in both
blocks, and set the `ACCENT_COLOR` tunable in
`docs/conf.py` to the nearest radix family name so the surfaces
Shibuya owns stay in sync. Nothing else needs to change.

## The palette

One accent color, used for links, active navigation, focus rings,
card hover borders and the landing page's main button (the brighter
`--rot-accent-fill`, which always carries `--rot-on-accent` text). Nothing else gets a brand color; admonition
colors are semantic and stay functional. All pairs pass WCAG AA:

The values themselves live only in `tokens.css`; this page cites
tokens by name so the numbers cannot go stale. Contrast ratios from
the last audit, measured against `--rot-bg` in each mode:

```{list-table}
:header-rows: 1

* - Token
  - Light mode
  - Dark mode
* - `--rot-accent`
  - 5.05:1
  - 9.11:1
* - `--rot-text`
  - 17.42:1
  - 15.54:1
* - `--rot-text-muted`
  - 6.06:1
  - 7.34:1
```

Run a contrast check on both modes whenever a token changes; the
accent must keep at least 4.5:1 against `--rot-bg` and
`--rot-surface`.

## Radius, typography, motion

- One corner radius scale: `--rot-radius` (4px), applied to cards,
  dropdowns, images, tables and code blocks. Do not introduce a
  second radius.
- Two typefaces, each behind one variable: `--rot-font-sans`
  (Chivo) for body and headings, `--rot-font-mono` (Chivo Mono) for
  code. Both are self-hosted in `_static/fonts` with
  `font-display: swap`, with system fonts as the fallback. No Google
  Fonts.
- Prose is capped near 72 characters. Tables and code blocks may
  use the full content column but scroll inside their own
  containers; the page never scrolls sideways.
- Motion is limited to hover and focus transitions under 150ms
  (`--rot-transition`), with one exception: collapsing a sidebar is a
  movement rather than a state change, so it eases its width over
  240ms (`--rot-transition-panel`). The landing page's drifting
  sequences are the only continuous animation. All of it honours
  `prefers-reduced-motion`, which disables it.

## Which component for which job

```{list-table}
:header-rows: 1

* - Job
  - Component
* - A set of links to related pages, each with a one-line pitch
  - `grid` of `grid-item-card` (a card must link somewhere or
    group a real unit; no decorative cards)
* - The same task shown against different backends or languages
  - `tab-set` with one `tab-item` per variant
* - Advanced or rarely used options that would bloat the main flow
  - `dropdown`, closed by default
* - A dataset the reader will want to sort or search
  - `list-table` with `:class: sphinx-datatable`
* - A short warning or side note
  - A standard admonition (`note`, `warning`)
```

## Page skeleton for a new module section

Copy the structure of the {doc}`data access section </db/index>`,
the worked example:

1. An intent paragraph: what the subpackage is for and what design
   idea unifies it.
2. A card grid, one card per submodule, each with an
   `` {iconify}`tabler:...` `` icon, the resource name, and one
   line of description, linking to the generated API page.
3. A quickstart as tabs, one minimal working snippet per data
   source or task.
4. Optional: a sortable table for enumerable facts, dropdowns for
   tuning knobs.
5. Common patterns: three to five realistic end-to-end snippets.

The generated API pages come for free: add the module to
`docs/api/index.rst` and autosummary builds the reference page the
cards link to. Icons must be present in
`_static/vendor/iconify-preload.js` to render offline; extend that
file when you use a new icon (see the URL pattern in the file
header of `conf.py`'s iconify section).
