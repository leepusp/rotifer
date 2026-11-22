# Search results

The `rotifer.io.search` module defines one standard table for sequence search
results (one row per matched region) and converts PSI-BLAST, PyHMMER and
hmmer2table layouts to and from it. It also reads BLAST tabular output
(one row per HSP, the HSP table) and merges HSPs into regions (the merged
table) the way blast2table does. A native layout is the table a tool or older
parser writes; after conversion, tables from different tools can be
concatenated. `rotifer.devel.beta.blast.psiblast` uses this table.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:table-alias` rotifer.io.search
:link: /api/generated/rotifer.io.search
:link-type: doc

The standard columns, validation, conversion from and to each tool's
layout, and HSP merging for BLAST tabular output.
:::
::::

## Quickstart

Every snippet below runs offline with only pandas and numpy, from the root of
the rotifer checkout.

::::{tab-set}

:::{tab-item} PSI-BLAST HSPs

Read a BLAST tabular report (format 7) into the standard table, one row per
HSP:

```python
import rotifer.io.search as st

df = st.read_blast_tabular("test/data/search/src3.fmt7.tsv", source="fixture")
print(len(df))          # 913
print(sorted(df.iteration.unique().tolist()))  # [1, 2]
```

:::

:::{tab-item} Merge HSPs

Collapse the HSP table to one row per region (the blast2table default mode),
then convert back to blast2table's column layout:

```python
import rotifer.io.search as st

df = st.read_blast_tabular("test/data/search/src3.fmt7.tsv", source="fixture")
merged = st.collapse_hsps(df)
print(len(merged))      # 465
blast_back = st.to_blast(merged)
print(list(blast_back.columns))
# ['hit', 'query', 'hitstart', 'hitend', 'evalue', 'querycoverage',
#  'querystart', 'queryend', 'iteration', 'bitscore', 'length']
```

:::

:::{tab-item} PyHMMER

Build a minimal in-memory frame with PyHMMER column names, convert to the
standard table, and validate it:

```python
import pandas as pd
import rotifer.io.search as st

pyhmmer_df = pd.DataFrame({
    "aln_target_name": ["p1", "p2"],
    "aln_hmm_name": ["PF1", "PF1"],
    "i_evalue": [1e-10, 3e-3],
    "score": [50.0, 10.0],
    "aln_target_from": [11, 2],
    "aln_target_to": [50, 40],
    "aln_target_length": [100, 40],
    "aln_hmm_length": [45, 45],
    "env_from": [10, 1],
    "env_to": [52, 40],
    "aln_hmm_from": [3, 5],
    "aln_hmm_to": [42, 43],
    "source": ["m.hmm", "m.hmm"],
})
std = st.from_pyhmmer(pyhmmer_df, "hmmscan")
validated = st.validate(std)
print(len(validated))   # 2
print(list(validated.columns[:len(st.COLUMNS)]))
# ['sequence', 'model', 'start', 'end', 'evalue', 'cov', 'qstart', 'qend',
#  'iteration', 'score', 'talilen', 'qlen', 'estart', 'eend', 'program',
#  'tlen', 'source']
```

:::
::::

## The standard table

One row per matched region, with the columns in this order. `validate`
returns a copy with these columns first and any extra columns after them. It
casts the numeric columns, keeps text columns as they come and raises
ValueError on bad data.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Column
  - Type
  - Meaning
* - sequence
  - text
  - Identifier of the matched sequence (BLAST hit, HMMER target)
* - model
  - text
  - Identifier of the query side (BLAST query, HMMER profile)
* - start
  - int64
  - First aligned residue in ``sequence`` (1-based, inclusive)
* - end
  - int64
  - Last aligned residue in ``sequence`` (1-based, inclusive)
* - evalue
  - float64
  - E-value of the region (HMMER: independent E-value)
* - cov
  - float64
  - Aligned fraction of the search query, 0 to 1
* - qstart
  - Int64 (nullable)
  - First aligned position in ``model`` (residue or HMM state)
* - qend
  - Int64 (nullable)
  - Last aligned position in ``model``
* - iteration
  - int64
  - Search round (PSI-BLAST, jackhmmer); 1 otherwise
* - score
  - float64
  - Bit score of the region
* - talilen
  - int64
  - Length of the region aligned in ``sequence``: ``end - start + 1``
* - qlen
  - Int64 (nullable)
  - Full length of the search query
* - estart
  - int64
  - First residue of the envelope in ``sequence`` (BLAST: same as start)
* - eend
  - int64
  - Last residue of the envelope in ``sequence`` (BLAST: same as end)
* - program
  - text
  - Program that produced the row (psiblast, hmmscan, hmmsearch, ...)
* - tlen
  - Int64 (nullable)
  - Full length of ``sequence``
* - source
  - text
  - Database or model file searched, when known
```

Coordinates are 1-based and inclusive: the region is ``seq[start - 1:end]``.
The search query is what the search was run with: the protein for PSI-BLAST
and hmmscan, the profile for hmmsearch. ``cov`` and ``qlen`` always refer to
it, so for hmmsearch ``cov`` is the fraction of the model that was aligned.

Differences from hmmer2table:
- ``talilen`` is ``end - start + 1``; hmmer2table writes ``end - start``.
- ``tlen`` and ``source`` are new columns.

## Conversions

Each ``from_*`` function returns the standard columns first, followed by any
input columns that have no standard equivalent. Those extra columns survive a
round trip back to the native layout. A blast2table-style table comes back
unchanged; the HMMER layouts may gain filled-in values, such as an iteration
of 1 where it was missing.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Native layout
  - To standard
  - Back
  - Notes
* - PSI-BLAST (blast2table columns, ``BLAST_COLUMNS``)
  - ``from_blast(df, program="psiblast", source=None)``
  - ``to_blast(df)``
  - ``qlen`` is missing (blast2table does not keep it); envelope equals aligned region
* - PyHMMER (``PYHMMER_COLUMNS``, renamed or raw)
  - ``from_pyhmmer(df, program, source=None)``
  - ``to_pyhmmer(df, rename=None)``
  - ``program`` is "hmmscan" or "hmmsearch"; ``qstart``/``qend`` need ``aln_hmm_from``/``aln_hmm_to``
* - hmmer2table (``HMMER2TABLE_COLUMNS``)
  - ``from_hmmer2table(df, source=None)``
  - ``to_hmmer2table(df)``
  - Corrects ``talilen`` to ``end - start + 1``; ``tlen`` is missing
* - BLAST tabular (format 7, ``BLAST_OUTFMT``)
  - ``read_blast_tabular(path_or_text, source=None)``
  - ``to_blast(df)``, after ``collapse_hsps`` if you want merged rows
  - One standard row per HSP; ``iteration`` from "# Iteration: N" comments
```

## Common patterns

Concatenate results from different tools after converting them:

```python
import pandas as pd
import rotifer.io.search as st

blast = pd.DataFrame({
    "hit": ["WP_1"], "query": ["Q"], "hitstart": [5], "hitend": [100],
    "evalue": [1e-20], "querycoverage": [0.96], "querystart": [1], "queryend": [96],
    "iteration": [1], "bitscore": [95.1], "length": [300],
})
hmmer2 = pd.DataFrame({
    "sequence": ["p1"], "model": ["PF1"], "start": [11], "end": [50],
    "evalue": [1e-10], "cov": [0.4], "qstart": [3], "qend": [42],
    "iteration": [1], "score": [50.0], "talilen": [39], "qlen": [100],
    "estart": [10], "eend": [52], "program": ["hmmscan"],
})
pyhmmer = pd.DataFrame({
    "aln_target_name": ["p2"], "aln_hmm_name": ["PF1"], "i_evalue": [3e-3],
    "score": [10.0], "aln_target_from": [2], "aln_target_to": [40],
    "aln_target_length": [40], "aln_hmm_length": [45], "env_from": [1],
    "env_to": [40], "aln_hmm_from": [5], "aln_hmm_to": [43],
})
combined = st.validate(pd.concat([
    st.from_blast(blast), st.from_hmmer2table(hmmer2),
    st.from_pyhmmer(pyhmmer, "hmmscan"),
], ignore_index=True))
print(len(combined))  # 3
print(sorted(combined.program.unique().tolist()))  # ['hmmscan', 'psiblast']
```

Round trips keep extra columns:

```python
import pandas as pd
import rotifer.io.search as st

native = pd.DataFrame({
    "hit": ["WP_1"], "query": ["Q"], "hitstart": [5], "hitend": [100],
    "evalue": [1e-20], "querycoverage": [0.96], "querystart": [1], "queryend": [96],
    "iteration": [1], "bitscore": [95.1], "length": [300], "note": ["a"],
})
back = st.to_blast(st.from_blast(native))
print(list(back.columns))
# ['hit', 'query', 'hitstart', 'hitend', 'evalue', 'querycoverage',
#  'querystart', 'queryend', 'iteration', 'bitscore', 'length', 'note']
```

`collapse_hsps` follows blast2table's choice of rows and coordinates, but
reports the real E-values, scores and iterations, which differ from the rounded
numbers blast2table prints. ``threshold`` is the inclusion E-value (default 0.01). ``overlap`` decides when
two HSPs are the same region: a float is a fraction of the longer HSP (default
0.4), an int is a number of shared residues.

```python
import rotifer.io.search as st

hsps = st.read_blast_tabular("test/data/search/src3.fmt7.tsv", source="fixture")
merged = st.collapse_hsps(hsps, threshold=0.001, overlap=0.5)
print(len(merged))  # 466
```