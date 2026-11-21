# Genomes and neighbourhoods

The `rotifer.genome` package turns annotated genomes into tables. It
reads GenBank and GFF files, writes one row per feature into a
`NeighborhoodDF` (a pandas DataFrame with gene neighbourhood methods),
and works with blocks: the stretches of a genome around the proteins
you are studying. The NCBI cursors in [Data access](../db/index.md)
return the same tables, so everything on this page also applies to
neighbourhoods downloaded with `rotifer.db`.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:table` NeighborhoodDF
:link: /api/generated/rotifer.genome.data
:link-type: doc

The feature table and its methods: find neighbourhoods, trim blocks,
compact views, block reports, networks and block clustering.
:::

:::{grid-item-card} {iconify}`tabler:file-import` Genome files
:link: /api/generated/rotifer.genome.io
:link-type: doc

Read GenBank and GFF files into Biopython records, with assembly
identifiers set from the file name or a rule you give.
:::

:::{grid-item-card} {iconify}`tabler:table-import` Records to tables
:link: /api/generated/rotifer.genome.utils
:link-type: doc

Convert Biopython records into a NeighborhoodDF, one row per feature,
with custom column extractors.
:::

:::{grid-item-card} {iconify}`tabler:layout-list` Block presets
:link: /api/generated/rotifer.genome.blocks
:link-type: doc

The column choices behind `to_blocks()` and `to_yaml()` for annotated
neighbourhood tables.
:::

:::{grid-item-card} {iconify}`tabler:file-text` Block reports guide
:link: ../pandas/blocks
:link-type: doc

How block reports are laid out, read back and customised.
:::

:::{grid-item-card} {iconify}`tabler:topology-star-3` Networks
:link: /api/generated/rotifer.graph.network
:link-type: doc

The generic edge builders that `to_network()` and `jaccard()` call.
:::
::::

## Quickstart

The first tab downloads live data, so run it on a machine with
network access. The other tabs work on the table it returns, `df`.

::::{tab-set}

:::{tab-item} Fetch

Collect the genomic region around the genes that encode two proteins:

```python
from rotifer.db import ncbi

gnc = ncbi.GeneNeighborhoodCursor(before=5, after=5, progress=True)
df = gnc.fetchall(["WP_012291365.1", "WP_013208129.1"])
print(df[df["query"] == 1].filter(["pid", "block_id", "assembly"]))
```
:::

:::{tab-item} Local files

Read GenBank files from disk. The assembly accession comes from NCBI
style file names (`GCF_000005845.2_ASM584v2_genomic.gbff`):

```python
from rotifer.genome import io as rgio
from rotifer.genome.utils import seqrecords_to_dataframe

records = rgio.parse(["GCF_000005845.2_ASM584v2_genomic.gbff"], "genbank")
genome = seqrecords_to_dataframe(records)
genome["query"] = genome.pid.isin(["NP_414542.1"]).astype(int)
```
:::

:::{tab-item} Neighbours

Search whole replicons for the features around each target, counting
only features of the target's type and wrapping around circular
chromosomes:

```python
hood = genome.neighbors(["query == 1"], before=3, after=3)
same = genome.neighbors(["query == 1"], before=3, after=3, strand="same")
```

Shrink blocks you already have to the rows around a target:

```python
core = df.trim_blocks(["query == 1"], before=2, after=2)
```
:::

:::{tab-item} Compact

One line per block, with arrows for strands and `*` on the query:

```python
print(df.compact(annotation="locus").compact.head())
```
:::

:::{tab-item} Reports

Write blocks as text tables or as a YAML document. The presets expect
an annotated table (sequence clusters, Pfam domains, taxonomy); see
the [block reports guide](../pandas/blocks.md):

```python
annotated.to_blocks(buf="blocks.txt")
annotated.to_yaml(buf="blocks.yaml")
```
:::

:::{tab-item} Networks

Count which domains sit next to each other, and group blocks that
share protein clusters. Both need annotation columns (`pfam`, `c80e3`):

```python
edges = annotated.to_network("pfam")
groups = annotated.jaccard("c80e3", min_jaccard=0.3)
```
:::
::::

## The NeighborhoodDF columns

The table is sortable and searchable. `neighbors()` rewrites
`block_id`, `query`, `origin`, `is_fragment` and `internal_id`;
`trim_blocks()` rewrites only `query`. The numbering described for
`feature_order` and `internal_id` is the one
`seqrecords_to_dataframe()` writes.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Column
  - Group
  - Meaning
* - nucleotide
  - nucleotide
  - Sequence (contig, chromosome or plasmid) accession
* - start, end
  - location
  - 1-based, inclusive feature coordinates; a feature that crosses the origin of a circular replicon has one row per part
* - strand
  - location
  - 1 or -1
* - nlen
  - nucleotide
  - Sequence length in nucleotides
* - block_id
  - block
  - Neighbourhood identifier; `neighbors()` writes `nucleotide:start-end`
* - query
  - block
  - 1 for the rows that anchored the neighbourhood, 0 otherwise
* - pid
  - feature
  - Protein accession (CDS only)
* - type
  - feature
  - Feature type (CDS, PSE for pseudogenes, tRNA, gene, ...)
* - plen
  - feature
  - Protein length in residues for CDS, otherwise feature length
* - locus
  - feature
  - Locus tag
* - seq_type
  - nucleotide
  - chromosome or plasmid
* - assembly
  - nucleotide
  - Genome assembly accession
* - gene, product
  - feature
  - Gene name and product description
* - origin
  - feature
  - 1 when the feature or block runs through the origin of replication; `neighbors()` writes 2 when a block covers a whole circular replicon
* - topology
  - nucleotide
  - circular or linear
* - taxid, organism, lineage, classification
  - taxonomy
  - NCBI taxonomy identifier, organism name, preferred-taxa lineage and full classification
* - feature_order
  - feature
  - Position among features of the same type on the sequence, from 0
* - internal_id
  - feature
  - Position among all features on the sequence, from 0; after `neighbors()`, the row number in the result
* - is_fragment
  - nucleotide
  - True when a block does not cover its whole sequence
* - replaced
  - feature
  - Set by the gene neighbourhood cursors: the query protein each row's protein stands for
```

## Choosing a method

:::{dropdown} neighbors() or trim_blocks()?

`neighbors()` starts from a whole genome table. It counts features of
the same type as each target (so `before=3` means three CDS, however
many tRNA or gene rows lie between), wraps around circular replicons,
merges windows closer than `min_block_distance` and names the new
blocks `nucleotide:start-end`. The NCBI cursors call it.

`trim_blocks()` starts from blocks you already have. It counts every
row, never leaves a block, keeps `block_id` and only rewrites
`query`. Use it to look closer at part of a neighbourhood table.
:::

:::{dropdown} Assembly identifiers

GenBank records carry no assembly accession unless their DBLINK says
so. `parse()` takes it from the file name by default
(`GenbankAssemblyFromFilename`); pass a string to set one for every
record, any callable taking `(path, record)` to compute it, or `None`
to keep what the file says.
:::
