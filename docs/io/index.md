# Reading tool output

The `rotifer.io` package parses the reports of sequence analysis programs.
Each submodule reads one format; import the one you need, for example
`from rotifer.io import hmmer`.

The BLAST, HMMER and Phobius parsers share a shape: a stream function
(`iter_blast`, `parse_hmmer`, `parse_phobius`) yields rows of strings, which
the command line tools write, and `read_<tool>(source)` loads the same rows
into a pandas DataFrame. Their `read_*` functions take a file name, a list,
`"-"` for stdin, `.gz` or `.bz2` files, or open handles. The HH-suite reader
takes a report, a list or a directory, and the DALI readers have their own
inputs (see the table below).

The command line tools are described on the {doc}`CLI page </cli/index>`.
A standard table for search results, with conversions between tool
layouts, is on the {doc}`search results page </io/search>`.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:tools` rotifer.io
:link: /api/generated/rotifer.io
:link-type: doc

Shared helpers: `open_input`, `each_source`, `split_fields`,
`written_as_fraction` and `strip_suffix`.
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.blast
:link: /api/generated/rotifer.io.blast
:link-type: doc

BLAST+ pairwise text reports (blastp, psiblast, rpsblast). One row per
merged HSP; `swap=True` turns RPS-BLAST into a domain table.
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.hmmer
:link: /api/generated/rotifer.io.hmmer
:link-type: doc

HMMER 3 reports (hmmscan, hmmsearch, phmmer, jackhmmer). One row per
domain hit.
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.phobius
:link: /api/generated/rotifer.io.phobius
:link-type: doc

Phobius predictions (signal peptides and transmembrane segments).
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.hhsuite
:link: /api/generated/rotifer.io.hhsuite
:link-type: doc

HH-suite `.hhr` reports (hhblits, hhsearch, hhalign). One row per hit.
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.dali
:link: /api/generated/rotifer.io.dali
:link-type: doc

DALI server TSV export as a star alignment (needs
`rotifer.devel.beta.sequence`).
:::

:::{grid-item-card} {iconify}`tabler:file-text` rotifer.io.dalidat
:link: /api/generated/rotifer.io.dalidat
:link-type: doc

DaliLite `.dat` files (`dat2fa`).
:::

:::{grid-item-card} {iconify}`tabler:table-alias` rotifer.io.search
:link: /api/generated/rotifer.io.search
:link-type: doc

The standard columns, validation, conversion from and to each tool's
layout, and HSP merging for BLAST tabular output.
:::

:::{grid-item-card} {iconify}`tabler:book-2` Search results
:link: search
:link-type: doc

The curated page for the standard search table and its conversions.
:::
::::

## Quickstart

Every snippet below runs offline from the root of the rotifer checkout, on
reports bundled with the tests.

::::{tab-set}

:::{tab-item} HMMER

Read an hmmscan report into a domain table, cleaning model versions:

```python
from rotifer.io.hmmer import read_hmmer

df = read_hmmer("test/data/golden/inputs/hmmscan.out", clean={"model": "version"})
print(df.iloc[:2, :6])
#       sequence          model  start  end        evalue   cov
# 0  NP_414542.1     Leader_Thr      1   21  1.200000e-10  1.00
# 1  NP_414543.1  Homoserine_dh    296  353  7.300000e+03  0.07
```
:::

:::{tab-item} BLAST

Read a blastp report (one row per merged HSP):

```python
from rotifer.io.blast import read_blast

df = read_blast("test/data/golden/inputs/blastp.out")
print(df.head(3).to_string())
#      hit_name   query_name  hit_start  hit_end        evalue  query_coverage  query_start  query_end  iteration    bits  hit_length
# 0  NP_414542.1  NP_414542.1          1       21  1.000000e-09           1.000            1         21          1    38.1          21
# 1  NP_414543.1  NP_414543.1          1      820  0.000000e+00           1.000            1        820          1  1687.0         820
# 2  NP_418414.1  NP_414543.1        601      661  2.200000e+00           0.072          627        685          1    22.3        1342
```

Read an rpsblast report as a domain table (`swap=True` exchanges query
and hit; `clean=["cdd"]` extracts domain names from descriptions):

```python
from rotifer.io.blast import read_blast

domains = read_blast("test/data/golden/inputs/rpsblast.out", swap=True, clean=["cdd"])
print(domains.loc[:1, ["hit_name", "query_name", "hit_start", "hit_end", "evalue"]].to_string())
#    hit_name                      query_name  hit_start  hit_end  evalue
# 0  NP_414543.1                            thrA          1      820     0.0
# 1  NP_414543.1  Aspartokinase_Aspartate_kinase          3      461     0.0
```
:::

:::{tab-item} Phobius

Phobius gives no E-value, so `evalue` sets the value written in that column.
Pick one that passes the E-value filter you will use later:

```python
from rotifer.io.phobius import read_phobius

df = read_phobius("test/data/golden/inputs/phobius.out", evalue=0.0101)
print(df.head(2).to_string())
#       sequence domain  start  end  evalue
# 0  NP_414546.1     SP      1   23  0.0101
# 1  NP_414548.1     TM     12   30  0.0101
```
:::

:::{tab-item} HH-suite

Read an `.hhr` report:

```python
from rotifer.io.hhsuite import read_hhr

df = read_hhr("test/io/data/sample.hhr")
print(df[["hit", "qstart", "qend", "probability"]].to_string(index=False))
#       hit  qstart  qend  probability
# PF00001.1       3    47         99.9
# PF00002.2      10    29         50.0
```
:::
::::

## Supported formats

The table is sortable and searchable.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Format
  - Reader
  - Stream function
  - Command
  - Output
* - BLAST+ pairwise (blastp, psiblast, rpsblast)
  - `rotifer.io.blast.read_blast`
  - `rotifer.io.blast.iter_blast`
  - `blast2table`
  - One row per merged HSP (`COLUMNS`); `swap=True` gives a domain table
* - HMMER 3 (hmmscan, hmmsearch, phmmer, jackhmmer)
  - `rotifer.io.hmmer.read_hmmer`
  - `rotifer.io.hmmer.parse_hmmer`
  - `hmmer2table`
  - One row per domain hit (`COLUMNS`)
* - Phobius (default long format)
  - `rotifer.io.phobius.read_phobius`
  - `rotifer.io.phobius.parse_phobius`
  - `phobius2table`
  - One row per SP or TM segment (`COLUMNS`: sequence, domain, start, end, evalue)
* - HH-suite `.hhr` (hhblits, hhsearch, hhalign)
  - `rotifer.io.hhsuite.read_hhr` / `hhsuite2pandas`
  - `rotifer.io.hhsuite.parse_hhr` (one DataFrame per report)
  - `hhsuite2table`
  - One row per hit (`default_columns`)
* - DALI server TSV
  - `rotifer.io.dali.dali_tsv`
  - (none)
  - (none)
  - Star alignment (`rotifer.devel.beta.sequence` object)
* - DaliLite `.dat`
  - (none)
  - `rotifer.io.dalidat.parse_dat` (identifier, description, sequence tuples), `rotifer.io.dalidat.dat2fasta` (FASTA records)
  - `dat2fa`
  - FASTA, one record per chain
```

## Common patterns

The stream functions yield one row of strings at a time, the same rows a
`read_*` function collects. HMMER and Phobius rows come out as the report is
read. `iter_blast` merges HSPs by default, so it reads the whole report first;
with `nonredundant=False` it yields every HSP as it goes:

```python
from rotifer.io.blast import iter_blast

with open("test/data/golden/inputs/blastp.out") as handle:
    for row in iter_blast(handle):
        print(row)
        break
# ['NP_414542.1', 'NP_414542.1', '1', '21', '1.000e-09', '1.000', '1', '21', '1', '38.1', '21']
```

To build domain architectures from several tools at once, concatenate their
tables and pass the result to
`rotifer.interval.architecture.domain2architecture`. It finds the columns by
name (`sequence`, `domain` or `model`, `start`, `end`, `evalue`), so HMMER and
Phobius tables go in as they are. An RPS-BLAST table built with `swap=True`
has the protein in `hit_name` and the domain in `query_name`, so rename those:

```python
import pandas as pd
from rotifer.io.blast import read_blast
from rotifer.io.hmmer import read_hmmer
from rotifer.io.phobius import read_phobius
from rotifer.interval.architecture import domain2architecture

hmmer = read_hmmer("test/data/golden/inputs/hmmscan.out")
phobius = read_phobius("test/data/golden/inputs/phobius.out")
rpsblast = read_blast("test/data/golden/inputs/rpsblast.out", swap=True, clean=["cdd"])
rpsblast = rpsblast.rename(columns={"hit_name": "sequence", "query_name": "domain",
                                    "hit_start": "start", "hit_end": "end"})

regions = pd.concat([hmmer, phobius, rpsblast], ignore_index=True)
arch = domain2architecture(regions, evalue=0.01)
print(arch[arch.architecture.str.contains("TM|SP")][["ID", "architecture"]].to_string(index=False))
#          ID      architecture
# NP_414546.1        SP+DUF2502
# NP_414551.1 TM+TM+TM+TM+TM+TM
```

Tune HSP merges in `read_blast`. `threshold` (default 0.01) is the
inclusion E-value: a kept region below it is widened by overlapping HSPs
instead of replaced. `overlap` (default `"0.4"`) is the fraction of each
HSP's length (if written with a decimal point) or a residue count that two
HSPs must share on both sequences to be merged. `nonredundant=False` yields
every HSP without merging.

Clean identifiers with `clean=`. `read_hmmer` accepts
`clean={"model": "version", "sequence": "profiledb"}` (or functions) to
strip version suffixes. `read_blast` accepts
`clean=["cdd", "genbank2acc"]` (or functions) for query and hit names.
`read_hhr` does not clean identifiers; apply your own function to the
`query` and `hit` columns after reading.

## Removed

`rotifer.io.base` and `rotifer.io.fileinput` are gone. Nothing used them. The
{doc}`changelog </changelog>` lists every change to `rotifer.io`.

```{toctree}
:hidden:

search
```