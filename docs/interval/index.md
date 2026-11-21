# Intervals and architectures

The `rotifer.interval` package works with sequence regions in pandas
DataFrames. For each sequence, `select_regions` picks the best
non-overlapping set of regions. `domain2architecture` writes the kept
regions as a domain architecture: names joined by `+`, as in `A+B`,
with nested regions in brackets, as in `A[B]`. `architecture2table`
expands an architecture table back to one row per region. The
`domain2architecture` and `architecture2table` shell commands are
documented on the {doc}`command line tools page </cli/index>`.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:layout-columns` Region selection
:link: /api/generated/rotifer.interval.architecture
:link-type: doc

Select regions, load domain tables, and build or expand domain
architectures.
:::

:::{grid-item-card} {iconify}`tabler:arrows-split` Interval utilities
:link: /api/generated/rotifer.interval.utils
:link-type: doc

Find the gaps between regions. Also holds the two deprecated region
filters.
:::
::::

## Quickstart

Every example runs offline with pandas and NumPy.

::::{tab-set}

:::{tab-item} Select regions

A has the lowest E-value and is kept first. B overlaps A by 11
residues, more than 40% of B's 20, so the default drops B. C does not
overlap A. Raising the limit to 15 residues keeps all three:

```python
import pandas as pd
from rotifer.interval.architecture import select_regions

regions = pd.DataFrame({
    "sequence": ["p1", "p1", "p1"],
    "domain": ["A", "B", "C"],
    "start": [1, 50, 70],
    "end": [60, 69, 120],
    "evalue": [1e-10, 1e-5, 1e-3],
})

print(select_regions(regions)["domain"].tolist())
# ['A', 'C']
print(select_regions(regions, overlap=15)["domain"].tolist())
# ['A', 'B', 'C']
```
:::

:::{tab-item} Architectures

With a residue overlap limit, a region nested inside another can stay.
Here B lies inside A, and C follows both:

```python
import pandas as pd
from rotifer.interval.architecture import domain2architecture

regions = pd.DataFrame({
    "sequence": ["p1", "p1", "p1"],
    "domain": ["A", "B", "C"],
    "start": [1, 20, 70],
    "end": [60, 40, 120],
    "evalue": [1e-10, 1e-5, 1e-3],
})

architectures = domain2architecture(regions, overlap="30")
print(architectures[["ID", "architecture"]].to_dict("records"))
# [{'ID': 'p1', 'architecture': 'A[B]+C'}]
```

The result also has `regions` and `evalue` columns. Use `read_regions`
to load headerless domain tables from file names, lists, open files, or
compressed files. `architecture_table` accepts regions that are already
loaded and filtered.
:::

:::{tab-item} Gaps

`complement` merges overlapping regions before listing the gaps. The
`length` column gives each sequence length, so the trailing gap is
included:

```python
import pandas as pd
from rotifer.interval.utils import complement

regions = pd.DataFrame({
    "ID": ["s1", "s1", "s1"],
    "start": [10, 20, 50],
    "end": [30, 25, 80],
    "length": [100, 100, 100],
})

gaps = complement(regions)
print(gaps.to_dict("records"))
# [{'ID': 's1', 'start': 1, 'end': 9, 'length': 9},
#  {'ID': 's1', 'start': 31, 'end': 49, 'length': 19},
#  {'ID': 's1', 'start': 81, 'end': 100, 'length': 20}]
```
:::
::::

## Presets

`PRESETS` holds `select_regions` arguments for common table formats.
Each preset uses `overlap=0.4`. Pass one with
`select_regions(df, **PRESETS['hmmer'])`.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Preset
  - reference
  - domain
  - start / end
  - Ranking
  - Use for
* - `hmmer`
  - `sequence`
  - `model`
  - `estart` / `eend`
  - Lowest `evalue`, then longest region
  - `hmmer2table` or `rotifer.io.search` output
* - `hhsuite`
  - `query`
  - `hit`
  - `qstart` / `qend`
  - Largest `sum_probs`, `probability`, `score`, then longest region
  - HH-suite tables
* - `archtsv`
  - `ID`
  - `domain`
  - `start` / `end`
  - Lowest `evalue`, then longest region
  - Expanded architecture tables
```

## Options

The main options of `select_regions` are:

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Option
  - Default
  - Meaning
* - `reference`
  - `'sequence'`
  - Column or columns naming the sequence of each region.
* - `domain`
  - `None`
  - Region-name column used by `ignore_overlap`. Uses `domain` when present, otherwise `model`.
* - `start`, `end`
  - `'start'`, `'end'`
  - Columns with inclusive region coordinates.
* - `criteria`
  - `None`
  - Ranking columns in precedence order, mapped to `True` for smallest first or `False` for largest first. `'region_length'` means `end - start + 1`, unless `df` already has a `region_length` column. The default is lowest `evalue`, then longest, then smallest start.
* - `by_length`
  - `False`
  - Rank longest regions first. Cannot be combined with `criteria`.
* - `overlap`
  - `0.4`
  - A float is a fraction of the shorter region, an int is a number of residues, and a string is read as written. Limits are exclusive.
* - `nested`
  - `True`
  - Allow a region inside a kept region when the overlap limit allows it. `False` always drops contained regions.
* - `ignore_overlap`
  - `()`
  - Domain names that are always kept and never drop another region. `'all'` keeps every region.
```

After the ranking criteria, ties go to the smallest start and then input
order. The result keeps all columns and the original index, sorted by
reference, numeric start, and numeric end.

## Common patterns

Use a preset and override one argument without changing `PRESETS`. The
default fraction drops B here; allowing 10 shared residues keeps it:

```python
import pandas as pd
from rotifer.interval.architecture import PRESETS, select_regions

regions = pd.DataFrame({
    "sequence": ["p1", "p1"],
    "model": ["A", "B"],
    "estart": [1, 52],
    "eend": [60, 71],
    "evalue": [1e-10, 1e-5],
})

print(select_regions(regions, **PRESETS["hmmer"])["model"].tolist())
# ['A']
options = {**PRESETS["hmmer"], "overlap": 10}
print(select_regions(regions, **options)["model"].tolist())
# ['A', 'B']
```

Rank by region length instead of E-value:

```python
import pandas as pd
from rotifer.interval.architecture import select_regions

regions = pd.DataFrame({
    "sequence": ["p1", "p1"],
    "domain": ["A", "B"],
    "start": [1, 10],
    "end": [20, 100],
    "evalue": [1e-10, 1e-5],
})

print(select_regions(regions)["domain"].tolist())
# ['A']
print(select_regions(regions, by_length=True)["domain"].tolist())
# ['B']
```

Compare a fractional limit with a residue limit. A float, or a string
with a decimal point or a negative exponent such as `'1e-2'`, is a fraction; an integer or integer string is a
number of residues:

```python
import pandas as pd
from rotifer.interval.architecture import select_regions

regions = pd.DataFrame({
    "sequence": ["p1", "p1", "p1"],
    "domain": ["A", "B", "C"],
    "start": [1, 50, 70],
    "end": [60, 69, 120],
    "evalue": [1e-10, 1e-5, 1e-3],
})

print(select_regions(regions, overlap=0.4)["domain"].tolist())
# ['A', 'C']
print(select_regions(regions, overlap="15")["domain"].tolist())
# ['A', 'B', 'C']
```

List domains in `ignore_overlap` to always keep them. They never cause
another region to be dropped either:

```python
import pandas as pd
from rotifer.interval.architecture import select_regions

regions = pd.DataFrame({
    "sequence": ["p1", "p1", "p1"],
    "domain": ["A", "B", "C"],
    "start": [1, 50, 70],
    "end": [60, 69, 120],
    "evalue": [1e-10, 1e-5, 1e-3],
})

kept = select_regions(regions, ignore_overlap=["B"])
print(kept["domain"].tolist())
# ['A', 'B', 'C']
```

Build an architecture and expand it back to one row per region:

```python
import pandas as pd
from rotifer.interval.architecture import architecture2table, domain2architecture

regions = pd.DataFrame({
    "sequence": ["p1", "p1", "p1"],
    "domain": ["A", "B", "C"],
    "start": [1, 20, 70],
    "end": [60, 40, 120],
    "evalue": [1e-10, 1e-5, 1e-3],
})

architectures = domain2architecture(regions, overlap="30")
expanded = architecture2table(architectures)
print(architectures["architecture"].tolist())
# ['A[B]+C']
print(expanded[["ID", "domain", "start", "end"]].to_dict("records"))
# [{'ID': 'p1', 'domain': 'A', 'start': 1, 'end': 60},
#  {'ID': 'p1', 'domain': 'B', 'start': 20, 'end': 40},
#  {'ID': 'p1', 'domain': 'C', 'start': 70, 'end': 120}]
```

The expanded table uses `TABLE_COLUMNS`: `ID`, `domain`, `start`, `end`,
and `evalue`.

## Migrating

`rotifer.interval.utils.filter_nonoverlapping_regions` and
`rotifer.interval.utils.architecture` are deprecated. Both raise
`FutureWarning` and will be removed on 2026-11-30 (issue #173).

```{list-table}
:header-rows: 1

* - Deprecated function
  - Replacement
* - `rotifer.interval.utils.filter_nonoverlapping_regions`
  - `select_regions(..., nested=False)`
* - `rotifer.interval.utils.architecture`
  - `domain2architecture` for architecture tables, or `select_regions` followed by a groupby.
```
