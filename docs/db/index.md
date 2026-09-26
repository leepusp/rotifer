# Data access

The `rotifer.db` package gives every biological data source the same
cursor interface: dictionary-style access for one entry
(`cursor[accession]`), lazy iteration for large queries
(`fetchone`), and bulk retrieval (`fetchall`). Cursors track which
identifiers could not be resolved, retry recoverable failures, and
can delegate to one another, so a query first checks fast local
resources and only then reaches remote services. This page is the
curated entry point; the generated API pages linked below hold the
full reference and are the template for documenting other ROTIFER
modules.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:dna-2` NCBI
:link: /api/generated/rotifer.db.ncbi
:link-type: doc

High level delegators for sequences, genomes, feature tables, gene
neighborhoods, IPG reports and taxonomy.
:::

:::{grid-item-card} {iconify}`tabler:cloud-download` NCBI Entrez
:link: /api/generated/rotifer.db.ncbi.entrez
:link-type: doc

E-utilities backends: EFetch and ELink access to every NCBI
database.
:::

:::{grid-item-card} {iconify}`tabler:server-2` NCBI FTP
:link: /api/generated/rotifer.db.ncbi.ftp
:link-type: doc

Checksum verified genome downloads from the NCBI FTP site.
:::

:::{grid-item-card} {iconify}`tabler:folders` Genome mirror
:link: /api/generated/rotifer.db.ncbi.mirror
:link-type: doc

The same genome cursors, served from a local mirror of the NCBI
genomes repository.
:::

:::{grid-item-card} {iconify}`tabler:search` Local FASTA
:link: /api/generated/rotifer.db.local.easel
:link-type: doc

Fast sequence retrieval from FASTA files indexed with Easel's
esl-sfetch.
:::

:::{grid-item-card} {iconify}`tabler:binary-tree-2` Taxonomy
:link: /api/generated/rotifer.db.local.ete3
:link-type: doc

NCBI Taxonomy lineages from the local ETE toolkit database.
:::

:::{grid-item-card} {iconify}`tabler:database` SQLite3
:link: /api/generated/rotifer.db.sql.sqlite3
:link-type: doc

Gene neighborhoods and IPG reports from local SQLite3 stores.
:::

:::{grid-item-card} {iconify}`tabler:arrows-exchange` UniProt
:link: /api/generated/rotifer.db.uniprot.webapi.idmapping
:link-type: doc

Identifier mapping through the UniProt REST API, with AlphaFold DB
and PDB structure links.
:::

:::{grid-item-card} {iconify}`tabler:stack-2` Cursor framework
:link: /api/generated/rotifer.db.core
:link-type: doc

The base classes behind every cursor: batching, retries, delegation
and missing entry bookkeeping.
:::
::::

## Quickstart

Every snippet below downloads live data, so run them on a machine
with network access. Cursors accept single accessions or any
iterable of accessions.

::::{tab-set}

:::{tab-item} Sequences

Fetch annotated protein records from NCBI:

```python
from rotifer.db import ncbi

sc = ncbi.SequenceCursor(database="protein")
records = sc.fetchall(["YP_009724395.1", "WP_063732599.1"])
for rec in records:
    print(rec.id, rec.description)
```
:::

:::{tab-item} Genomes

Download whole genomes, or just their feature tables, by assembly
accession:

```python
from rotifer.db import ncbi

gfc = ncbi.GenomeFeaturesCursor(progress=True)
df = gfc.fetchall(["GCA_018744545.1", "GCA_901308185.1"])
print(df.groupby("assembly").pid.nunique())
```
:::

:::{tab-item} Neighborhoods

Fetch the genomic region around the genes encoding target proteins:

```python
from rotifer.db import ncbi

gnc = ncbi.GeneNeighborhoodCursor(before=5, after=5, progress=True)
df = gnc.fetchall(["WP_012291365.1", "WP_013208129.1"])
print(df[df["query"] == 1].filter(["pid", "block_id", "assembly"]))
```
:::

:::{tab-item} Taxonomy

Resolve NCBI Taxonomy identifiers to named lineages, using the
local ETE toolkit database first:

```python
from rotifer.db import ncbi

tc = ncbi.TaxonomyCursor()
taxa = tc.fetchall([2599, 562])
print(taxa.filter(["taxid", "organism", "lineage"]))
```
:::

:::{tab-item} Local FASTA

Fetch sequences from an esl-sfetch indexed FASTA file, without
touching the network:

```python
from rotifer.db.local import easel

fc = easel.FastaCursor("/databases/fadb/nr/nr")
seqs = fc.fetchall(["WP_063732599.1"])
```
:::

:::{tab-item} UniProt

Map GenBank CDS accessions to UniProtKB and collect AlphaFold DB
structure links:

```python
from rotifer.db.uniprot.webapi import idmapping

links = idmapping.AF_link(["BAE76179.1", "AAC73502.1"])
print(links.filter(["From", "Entry", "urlAF"]))
```
:::
::::

## Supported databases

The table is sortable and searchable. "Network" marks resources
that contact a remote service at query time.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Resource
  - Module
  - Entry points
  - Identifier types
  - Network
  - Cached
* - NCBI Entrez (E-utilities)
  - {py:mod}`rotifer.db.ncbi.entrez`
  - SequenceCursor, FastaCursor, IPGCursor, TaxonomyCursor, NucleotideFeaturesCursor, GeneNeighborhoodCursor, elink
  - protein, nucleotide and taxonomy accessions
  - yes
  - no
* - NCBI FTP genomes
  - {py:mod}`rotifer.db.ncbi.ftp`
  - connection, GenomeCursor, GenomeFeaturesCursor, GeneNeighborhoodCursor
  - assembly accessions (GCA_, GCF_)
  - yes
  - downloads pass through the cache directory
* - NCBI genome mirror
  - {py:mod}`rotifer.db.ncbi.mirror`
  - GenomeCursor, GenomeFeaturesCursor, GeneNeighborhoodCursor
  - assembly accessions (GCA_, GCF_)
  - no
  - reads a local mirror tree
* - Indexed FASTA files (Easel)
  - {py:mod}`rotifer.db.local.easel`
  - FastaCursor
  - sequence accessions
  - no
  - builds .ssi indices on first use
* - NCBI Taxonomy (ETE toolkit)
  - {py:mod}`rotifer.db.local.ete3`
  - TaxonomyCursor
  - taxonomy identifiers
  - no, after the first ete3 setup
  - local SQLite copy of the taxonomy dump
* - SQLite3 stores
  - {py:mod}`rotifer.db.sql.sqlite3`
  - GeneNeighborhoodCursor, IPGCursor
  - protein accessions
  - no
  - the database file itself
* - UniProt REST API
  - {py:mod}`rotifer.db.uniprot.webapi.idmapping`
  - genbank_to_uniprot, AF_link, af_to_seq
  - EMBL/GenBank/DDBJ CDS and UniProtKB accessions
  - yes
  - no
* - NCBI delegators
  - {py:mod}`rotifer.db.ncbi`
  - SequenceCursor, FastaCursor, GenomeCursor, GenomeFeaturesCursor, GeneNeighborhoodCursor, IPGCursor, TaxonomyCursor, assemblies
  - mixed, backend dependent
  - only when local backends miss
  - through the local backends
```

## Tuning

:::{dropdown} Batching, threads and retries

Every remote cursor accepts `batch_size` (accessions per request),
`threads` (parallel workers) and `tries` (attempts per batch).
Entrez cursors cap `threads` at 3 without an NCBI API key and at 10
with one; set the `NCBI_API_KEY` environment variable to raise the
cap. When `batch_size` is unset, the input is divided evenly among
the workers.
:::

:::{dropdown} Giving up on unrecoverable errors

Cursors keep a `giveup` set of error message substrings. An error
matching any of them marks the affected accessions as unrecoverable
and stops further attempts. Gene neighborhood cursors give up on
HTTP 400 responses, on proteins without IPG reports and, unless
`eukaryotes=True`, on eukaryotic genomes.
:::

:::{dropdown} Local mirrors and SQLite acceleration

`ncbi.GeneNeighborhoodCursor` accepts `mirror` (one path or a list
of paths to local NCBI genome mirrors) and `save` (a SQLite3
database of previously processed neighborhoods). Both are queried
before any remote source. The `mirror` entry of the NCBI
configuration file makes the mirror the default for every genome
cursor.
:::

:::{dropdown} Inspecting failures

After any fetch, `cursor.missing` is a dataframe of unresolved
accessions with the last error message and whether a retry could
succeed. `missing_ids(retry=True)` lists the accessions worth
retrying; `remove_missing()` clears the registry.
:::

## Common patterns

Fetch protein sequences with a local first, remote second strategy:

```python
import rotifer.db as rdb

seqs = rdb.proteins(
    ["WP_063732599.1", "YP_009724395.1"],
    methods=["esl_sfetch", "entrez"],
)
```

Collect gene neighborhoods for many proteins, reusing one IPG
download:

```python
from rotifer.db import ncbi

query = ["WP_012291365.1", "WP_013208129.1", "WP_122330970.1"]
ic = ncbi.IPGCursor(progress=True)
ipgs = ic.fetchall(query)

gnc = ncbi.GeneNeighborhoodCursor(progress=True)
df = gnc.fetchall(query, ipgs=ipgs)
```

Stream large genome sets instead of holding them all in memory:

```python
from rotifer.db import ncbi

gfc = ncbi.GenomeFeaturesCursor(progress=True)
for table in gfc.fetchone(assembly_list):
    table[table.type == "CDS"].to_csv("cds.tsv", sep="\t", mode="a")
```

Annotate an assembly report with taxonomy and check what failed:

```python
from rotifer.db import ncbi

reports = ncbi.assemblies(targets=["refseq"])
tc = ncbi.TaxonomyCursor()
taxa = tc.fetchall(reports.taxid.unique().tolist())
print(tc.missing)
```

Attach AlphaFold models to a sequence object:

```python
from rotifer.db.uniprot.webapi import idmapping

with_structures = idmapping.af_to_seq(seqobj)
```
