# Tables as text

`rotifer.pandas` turns DataFrames into text you can read in a terminal or edit
by hand. It renders frames as aligned columns, pads and colours strings, writes
groups of rows as block reports in plain text or YAML and reads the YAML back,
and has a few display helpers. Everything is a plain function; no accessor
is registered on pandas objects. Padding uses
pyarrow and the YAML reports use pyyaml.

Block reports have their own guide: {doc}`blocks`.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:align-left` Text
:link: /api/generated/rotifer.pandas.text
:link-type: doc

`to_string` renders a frame as aligned columns, `pad` pads strings to a common
width and `to_color` adds terminal colours to each character.
:::

:::{grid-item-card} {iconify}`tabler:layout-rows` Blocks
:link: /api/generated/rotifer.pandas.blocks
:link-type: doc

`to_blocks` writes one text block per group, each with a two-line header and
the group's rows.
:::

:::{grid-item-card} {iconify}`tabler:file-code` YAML blocks
:link: /api/generated/rotifer.pandas.yaml_blocks
:link-type: doc

Markdown tables (`md_table`, `read_md_table`) and YAML block reports that
round trip (`to_yaml`, `read_yaml`).
:::

:::{grid-item-card} {iconify}`tabler:adjustments` Display
:link: /api/generated/rotifer.pandas.display
:link-type: doc

`optimize_memory_usage` picks smaller dtypes; `print_everything` and
`show_display_options` set and list the pandas display options.
:::

:::{grid-item-card} {iconify}`tabler:history` Functions (compatibility)
:link: /api/generated/rotifer.pandas.functions
:link-type: doc

The old import path. It re-exports `pad`, `to_color`, `config` and the display
helpers. New code imports from `rotifer.pandas`.
:::

:::{grid-item-card} {iconify}`tabler:book` Block reports guide
:link: blocks
:link-type: doc

How to write block reports, edit them by hand and read the edits back.
:::
::::

## Quickstart

Every example runs offline.

::::{tab-set}

:::{tab-item} Aligned text

```python
import pandas as pd
import rotifer.pandas as rp

df = pd.DataFrame({"gene": ["a", "bb", "ccc"], "score": [1.5, None, 3.0]})
print(rp.to_string(df, float_format="%.2f", na_rep="NA", align="right"))
#   : gene : score
# 0 :    a :  1.50
# 1 :   bb :    NA
# 2 :  ccc :  3.00

print(rp.pad(df["gene"], side="right", fillchar=".").tolist())
# ['a..', 'bb.', 'ccc']
```
:::

:::{tab-item} Colour

`to_color` wraps each character in a 256-colour ANSI code. Without `colors` it
uses the amino acid colours. Printing the `repr` shows the codes:

```python
import pandas as pd
import rotifer.pandas as rp

s = pd.Series(["AC", "CA"])
colored = rp.to_color(s, colors={"A": "196", "C": "46"})
print(repr(colored.iloc[0]))
# '\x1b[38;5;196mA\x1b[m\x1b[38;5;46mC\x1b[m'

# one colour for every character, on the background
print(repr(rp.to_color(s, colors="21", style="bg").iloc[1]))
# '\x1b[48;5;21mC\x1b[0m\x1b[48;5;21mA\x1b[0m'
```
:::

:::{tab-item} Markdown and YAML

```python
import pandas as pd
import rotifer.pandas as rp

df = pd.DataFrame({"g": ["A", "A", "B"], "id": ["p1", "p2", "p3"], "kb": [10, 20, 30]})
print(rp.md_table(df))
#     | g | id | kb |
#     | - | -- | -- |
#     | A | p1 | 10 |
#     | A | p2 | 20 |
#     | B | p3 | 30 |

text = rp.to_yaml(df, "g", stats={"n": ("id", "nunique")})
print(text)
# - parameters:
#     groupby: [g]
#
# # -----
# - key: [A]
#   stats: |
#     | n |
#     | - |
#     | 2 |
#   rows: |
#     | id | kb |
#     | -- | -- |
#     | p1 | 10 |
#     | p2 | 20 |
# ...

print(rp.read_yaml(text))
#    g  id  kb  n
# 0  A  p1  10  2
# 1  A  p2  20  2
# 2  B  p3  30  1
```
:::

:::{tab-item} Smaller frames

```python
import pandas as pd
import rotifer.pandas as rp

df = pd.DataFrame({"n": [1, 2, 3, 4], "tag": ["a", "a", "a", "b"], "val": [1.0, 2.5, 3.0, 4.0]})
small = rp.optimize_memory_usage(df)
print({column: str(dtype) for column, dtype in small.dtypes.items()})
# {'n': 'int8', 'tag': 'category', 'val': 'float32'}
```
:::
::::

## Functions

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Function
  - Module
  - Takes
  - Returns
  - Use
* - {py:func}`~rotifer.pandas.text.to_string`
  - text
  - DataFrame
  - str
  - Render a frame as aligned columns, with your choice of separator, missing value text and alignment.
* - {py:func}`~rotifer.pandas.text.pad`
  - text
  - DataFrame or Series
  - Same type
  - Pad every value of a column to the column's width, optionally counting the header.
* - {py:func}`~rotifer.pandas.text.to_color`
  - text
  - DataFrame or Series
  - Same type
  - Colour each character, foreground or background, from a map, a YAML file or one colour code.
* - {py:func}`~rotifer.pandas.blocks.to_blocks`
  - blocks
  - DataFrame, grouping columns
  - str, or None with `buf`
  - Write one text block per group.
* - {py:func}`~rotifer.pandas.blocks.default_header`
  - blocks
  - key, grouping columns, stats
  - str
  - The header `to_blocks` writes unless you pass your own.
* - {py:func}`~rotifer.pandas.blocks.default_block_rows`
  - blocks
  - DataFrame, grouping columns
  - DataFrame
  - Deprecated. Use `to_blocks`.
* - {py:func}`~rotifer.pandas.yaml_blocks.md_table`
  - yaml_blocks
  - DataFrame
  - str
  - Write a padded Markdown table.
* - {py:func}`~rotifer.pandas.yaml_blocks.read_md_table`
  - yaml_blocks
  - str
  - DataFrame
  - Read a Markdown table back.
* - {py:func}`~rotifer.pandas.yaml_blocks.to_yaml`
  - yaml_blocks
  - DataFrame, grouping columns, extra tables
  - str, or None with `buf`
  - Write a YAML block report with Markdown tables inside.
* - {py:func}`~rotifer.pandas.yaml_blocks.read_yaml`
  - yaml_blocks
  - text, path or file
  - DataFrame, or a dict of tables with `merge=False`
  - Read a report written by `to_yaml`, including hand edits.
* - {py:func}`~rotifer.pandas.display.optimize_memory_usage`
  - display
  - DataFrame
  - New DataFrame
  - Downcast integers, use float32 when the values fit and category for repetitive text.
* - {py:func}`~rotifer.pandas.display.print_everything`
  - display
  - display limits
  - None
  - Raise the pandas display limits so whole frames print.
* - {py:func}`~rotifer.pandas.display.show_display_options`
  - display
  - nothing
  - DataFrame
  - List the current pandas display options.
```

## Common patterns

For a quick look at a frame in the terminal, `print(rp.to_string(df))` shows
every row with ` : ` between columns. `print_everything()` raises the pandas
display limits for the rest of the session, so plain `print(df)` shows every
row and column too.

To review or curate groups of rows by hand, write them with `to_yaml`, edit the
file and read it back with `read_yaml`. The {doc}`blocks` guide walks through
it.

To show sequences aligned and coloured, pass `padding` to `to_color`. It pads
the column first, so the coloured values line up:

```python
import pandas as pd
import rotifer.pandas as rp

seqs = pd.DataFrame({"id": ["p1", "p2"], "seq": ["MKV", "MK"]})
seqs["seq"] = rp.to_color(seqs["seq"], padding="right", fillchar="-")
print(repr(seqs["seq"].iloc[1]))
# '\x1b[38;5;033mM\x1b[m\x1b[38;5;124mK\x1b[m\x1b[39;5;000m-\x1b[m'
```

Before saving a large frame, or keeping it in memory, run
`optimize_memory_usage` on it. It returns a new frame and leaves the input
alone.

## Migrating

These names still work but warn with a `FutureWarning`. They will be removed
on 2026-11-30 (issue #160).

```{list-table}
:header-rows: 1

* - Old
  - New
* - `rotifer.pandas.to_network`
  - {py:func}`rotifer.graph.network.adjacency_edges`
* - `to_blocks(..., apply=...)`
  - The `header`, `columns` and `stats` arguments of `to_blocks`
* - `default_header(df, ...)`, with a DataFrame first
  - `default_header(key, groupby, stats)`
* - `default_block_rows`
  - `to_blocks`
* - `padding_df` in `rotifer.devel.alpha.gian_func`, `rotifer.devel.alpha.rodolfo` and `rotifer.devel.beta.blast`
  - `rotifer.pandas.pad(df, side="right", header=True)`
```

```{toctree}
:hidden:

blocks
```
