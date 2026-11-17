---
name: rotifer
description: Lagoa: one aqua accent on green-tinted neutrals, Chivo and Chivo Mono, and the Filtro rotifer mark.
colors:
  accent: "#0b7869"
  accent-hover: "#075b50"
  accent-fill: "#2fb3a1"
  accent-tint: "#ddf1ec"
  on-accent: "#101413"
  bg: "#f7f8f6"
  surface: "#eef1ee"
  border: "#dce1dd"
  text: "#101413"
  text-muted: "#59605c"
  code-bg: "#eef1ee"
  dark-accent: "#4cc9b5"
  dark-accent-hover: "#7ddccb"
  dark-accent-tint: "#173a34"
  dark-bg: "#111413"
  dark-surface: "#181d1b"
  dark-border: "#28302d"
  dark-text: "#e8ece9"
  dark-text-muted: "#9da5a0"
  dark-code-bg: "#161b19"
typography:
  display:
    fontFamily: "Chivo, system-ui, sans-serif"
    fontSize: "clamp(3rem, 9vw, 6rem)"
    fontWeight: 800
    lineHeight: 0.95
    letterSpacing: "-0.04em"
  section:
    fontFamily: "Chivo, system-ui, sans-serif"
    fontSize: "clamp(1.75rem, 3.2vw, 2.5rem)"
    fontWeight: 800
    lineHeight: 1.1
    letterSpacing: "-0.025em"
  body:
    fontFamily: "Chivo, system-ui, sans-serif"
    fontSize: "1rem"
    fontWeight: 400
    lineHeight: 1.6
  code:
    fontFamily: "Chivo Mono, ui-monospace, monospace"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.65
rounded:
  base: "4px"
components:
  action-primary:
    backgroundColor: "{colors.accent-fill}"
    textColor: "{colors.on-accent}"
    rounded: "{rounded.base}"
    padding: "12px 22px"
  code-sample:
    backgroundColor: "{colors.code-bg}"
    textColor: "{colors.text}"
    rounded: "{rounded.base}"
---

# Design System: rotifer

## Overview

rotifer's site is a technical tool's home: calm neutrals, one aqua accent,
real code and real data. The identity comes from three decisions taken on
2026-11-17: the Lagoa palette, the Chivo + Chivo Mono type pair and
the Filtro logo, a minimal rotifer. The landing page is the one place
with an authored moment: real protein sequences drifting slowly behind the
headline "Genomics made easy". The docs stay quiet and technical, and their
polish comes from type, spacing and clear navigation.

The product context (audience, positioning, constraints) is in PRODUCT.md,
next to this file.

Key characteristics:

- One accent, aqua, used for links, focus, the active page and the primary
  button. Nothing else carries a brand color.
- Neutrals tinted slightly toward green, to sit with the aqua.
- Mono type only for code, commands, identifiers and sequences.
- A single 4px radius; hairline borders; no shadows.
- No external requests: fonts, icons and scripts are all self-hosted.

## Colors

All values live in `docs/_static/css/tokens.css`, in a light block and an
`html.dark` block. Every other stylesheet uses the tokens.

- Accent (`--rot-accent`, #0b7869 / dark #4cc9b5): links, focus rings,
  the current page in the sidebar, card hover borders. It's the deep teal on
  light backgrounds because the brighter aqua fails as text on white.
- Accent fill (`--rot-accent-fill`, #2fb3a1): filled shapes only, such
  as the landing's Install button and the logo's wheels. Always carries ink
  text (`--rot-on-accent`, 7.14:1), never white.
- Accent tint (`--rot-accent-tint`): text selection and highlights.
- Fixed surfaces (`--rot-ink`, `--rot-paper`, `--rot-aqua-on-ink`): the
  same in both modes, for the landing hero, which is always ink.
- Neutrals: `--rot-bg`, `--rot-surface`, `--rot-border`, `--rot-text`,
  `--rot-text-muted`, `--rot-code-bg`.

Contrast against `--rot-bg` (WCAG AA is 4.5:1):

| Token | Light | Dark |
|---|---|---|
| `--rot-accent` | 5.05:1 | 9.11:1 |
| `--rot-text` | 17.42:1 | 15.54:1 |
| `--rot-text-muted` | 6.06:1 | 7.34:1 |

For emphasis, use weight or size first. Color is for links, state and
the primary action.

## Typography

- Chivo for everything readable: headings at weights 700 to 800 with tight
  tracking, body at 400. One variable file covers every weight.
- Chivo Mono for code blocks, inline code, command and module names,
  and the hero sequences. Never as a "technical" costume for prose.
- Prose is capped near 72 characters (65 to 75 on the landing).
- Both faces are self-hosted from `docs/_static/fonts/`, built by
  `docs/_tools/build_fonts.py`.

## Layout

- Docs: Shibuya's three-column layout: navigation tree, page, and "On
  this page". Both side panels collapse to rails
  (`_static/sidebar-toggle.js`).
- Landing: a single column up to 76rem with a fluid gutter. The hero comes
  first, centred on ink (`--rot-ink`) in both color modes, with paper
  text (`--rot-paper`), the reverse lockup, and "Genomics" in the bright
  aqua (`--rot-aqua-on-ink`, 9.13:1). Its primary button is that aqua with
  ink text. After the hero, then "What it does" as rows separated by hairlines (text
  on the left, a real usage sample on the right), then a "Get started" band
  on the surface color. On narrow screens each row stacks.

## Elevation and shape

No shadows. Objects are separated by 1px `--rot-border` hairlines and the
surface color. Every rounded corner uses `--rot-radius` (4px).

## Motion

- Hover and focus transitions under 150ms (`--rot-transition`).
- Sidebar collapse eases its width over 240ms (`--rot-transition-panel`).
- Sequence drift (landing only): 14 rows of real protein sequence, in
  mono at 14% paper over the ink, move horizontally at different slow speeds
  (140 to 230 seconds per loop), alternating direction. A mask fades them at the
  edges and quiets them behind the centred text. Pure CSS;
  `sequence-drift.js` only pauses them when the hero is off
  screen.
- `prefers-reduced-motion` turns all of it off, and the drift stays as a
  still texture.

## Components

- Header: the horizontal lockup (`rotifer-lockup-horizontal-*.svg`) at
  30px, then Install, Data access, API, Contributing, License.
- Primary action: aqua fill, ink text, 4px radius. Secondary actions
  are underlined accent links.
- Code sample: code-bg surface, hairline border, mono 0.875rem. Console
  prompts (`$ `) are dimmed and not selectable.
- Name chips: module and tool names in mono inside a hairline outline.
- Footer: "Made at" followed by the psel horizontal lockup (color or
  reverse per mode), then the Sphinx/Shibuya credit.

## Logo

The kit is generated by `brand/logo-kit-src/build.py`, which is kept out of
git. The site uses copies in `docs/_static/brand/`: the horizontal lockup
and mark in color and reverse, the favicons, and the psel lockup for the
credit. Only the wheels take color. Never retype "rotifer" in a font; use
the outlined files.

## Do's and don'ts

Do:

- Show real modules, commands and sequences, and record where data comes
  from (the hero names its PDB entries).
- Keep landing copy in `docs/_data/landing.<lang>.yml` and markup in
  `docs/_templates/landing/`.
- Theme browser surfaces: selection, focus rings, scrollbars.

Don't:

- Add a second accent, gradients, glass or shadows.
- Put white text on the aqua fill.
- Add eyebrows or kicker labels above headings, or decorative section
  numbers.
- Load anything from a CDN or font service.
