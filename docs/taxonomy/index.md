# Taxonomy

The `rotifer.taxonomy` package reduces NCBI classifications such as
`Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia`
to short lineages that keep only clades on a preferred-taxa list. By default,
the names are lowercased and joined by `>`, giving
`bacteria>pseudomonadota>gammaproteobacteria>enterobacterales`.
The list ships as `share/rotifer/taxonomy/taxonomy.txt` and is read by
`load_preferred_taxa()`. The NCBI, UniProt and ete3 taxonomy cursors in
`rotifer.db` use `lineage()` to fill their `lineage` columns. The full
reference is in the generated API pages linked below.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:binary-tree-2` rotifer.taxonomy.utils
:link: /api/generated/rotifer.taxonomy.utils
:link-type: doc

Reduce classifications with `lineage()` and read preferred clades with
`load_preferred_taxa()`.
:::

:::{grid-item-card} {iconify}`tabler:list-details` rotifer.taxonomy
:link: /api/generated/rotifer.taxonomy
:link-type: doc

Package overview. Import the functions from `rotifer.taxonomy.utils`; the
package itself exports nothing.
:::
::::

## Quickstart

Run these examples from the repository root. They use only local data.

::::{tab-set}

:::{tab-item} A list

Pass an iterable of classifications to get a list of lineages in the same
order. A missing classification stays missing:

```python
from rotifer.taxonomy.utils import lineage

classifications = [
    "Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia",
    "Archaea; Methanobacteriota; Methanobacteria",
    None,
]
print(lineage(classifications))
# ['bacteria>pseudomonadota>gammaproteobacteria>enterobacterales',
#  'archaea>methanobacteria', None]
```
:::

:::{tab-item} A table

Read the fixture included in the repository and add a lineage for each row:

```python
import pandas as pd
from rotifer.taxonomy.utils import lineage

path = "test/data/taxonomy/ncbi_classifications.tsv"
df = pd.read_csv(path, sep="\t", comment="#")
df["lineage"] = lineage(df)
print(df[["taxid", "lineage"]].head())
```

With the shipped list, the first printed row pairs taxid `2160` with
`archaea>methanobacteria>methanobacteriales`. DataFrame input uses the `classification`
column by default and returns a Series with the same index.
:::

:::{tab-item} Your own clades

Pass your own preferred taxa and output separator. With `case=None`, matching
is case-sensitive and the output keeps the input case:

```python
from rotifer.taxonomy.utils import lineage

classifications = [
    "Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia"
]
result = lineage(
    classifications,
    preferred_taxa=["Bacteria", "Enterobacteriaceae", "Escherichia"],
    outsep="|",
    case=None,
)
print(result)
# ['Bacteria|Enterobacteriaceae|Escherichia']
```
:::
::::

## Functions

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Name
  - Input
  - Returns
  - Notes
* - {py:func}`rotifer.taxonomy.utils.lineage`
  - A Series, a DataFrame with a `classification` column, or an iterable of
    strings. A single string raises TypeError.
  - A Series for Series or DataFrame input, otherwise a list. Missing values
    become NaN in a Series and None in a list.
  - `preferred_taxa=None` loads the default list on every call. An empty
    iterable keeps nothing. Cost scales with the number of distinct
    classifications.
* - {py:func}`rotifer.taxonomy.utils.load_preferred_taxa`
  - An optional file path or data-file name. The default is
    `':taxonomy.taxonomy.txt'`.
  - A set of lowercase clade names, or an empty set if no file is found.
  - Reads every matching user and shipped file. Comment lines and blank
    lines are skipped.
* - ``loadPreferredTaxa()``
  - No arguments.
  - The same set as `load_preferred_taxa()`.
  - Deprecated. Emits FutureWarning and will be removed on 2026-11-30.
```

## Common patterns

Load the preferred taxa once and reuse the set when making several calls:

```python
from rotifer.taxonomy.utils import lineage, load_preferred_taxa

taxa = load_preferred_taxa()
bacteria = lineage(
    ["Bacteria; Bacillota; Bacilli"],
    preferred_taxa=taxa,
)
archaea = lineage(
    ["Archaea; Methanobacteriota; Methanobacteria"],
    preferred_taxa=taxa,
)
print(bacteria)
# ['bacteria>bacillota>bacilli']
print(archaea)
# ['archaea>methanobacteria']
```

Extend the shipped list by putting one clade per line in
`taxonomy/taxonomy.txt` under `~/.rotifer/share` or `$ROTIFER_DATA`.
The default loader reads these copies as well as the shipped one, so a user
file adds names to the shipped list. Names are stripped of surrounding
whitespace and lowercased.

Missing classifications stay missing. Series input also keeps its name,
input order and duplicate index labels:

```python
import pandas as pd
from rotifer.taxonomy.utils import lineage

classifications = pd.Series(
    ["Bacteria; Bacillota; Bacilli", None],
    index=["a", "a"],
    name="classification",
)
result = lineage(classifications)
print(result.iloc[0])
# bacteria>bacillota>bacilli
print(result.isna().tolist())
# [False, True]
print(result.index.tolist())
# ['a', 'a']
```

For other iterables, a missing classification gives `None`. An empty
classification, or one with no preferred clades, gives an empty string, so
you can still tell it apart from a missing one.

## Migrating

Replace `loadPreferredTaxa()` with `load_preferred_taxa()`. The old name
emits FutureWarning and will be removed on 2026-11-30 (issue #160).
`bin/rlineage` was removed. Call `lineage()` from `rotifer.taxonomy.utils`
instead.
