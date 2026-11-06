# Command line tools

Every ROTIFER command line tool exists twice: as a **library
function** that returns a pandas DataFrame (or yields rows) and as a
**shell command** of the same name that writes a tab-separated table
to standard output. The command only parses options and calls the
function, so whatever a command does can also be done from Python,
without subprocesses or temporary files. Commands read the files named on the command line or,
when none is given, standard input, so they chain with pipes; plain,
gzip and bzip2 files are read the same way. Messages go to standard
error, never into the table. This page is the curated entry point;
the cards below link to the generated reference of each library
module.

## Modules

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {iconify}`tabler:search` HMMER
:link: /api/generated/rotifer.io.hmmer
:link-type: doc

`hmmer2table`: domain tables from hmmscan and hmmsearch reports.
:::

:::{grid-item-card} {iconify}`tabler:file-database` BLAST+
:link: /api/generated/rotifer.io.blast
:link-type: doc

`blast2table`: HSP tables from BLASTP, PSI-BLAST and RPS-BLAST text
reports, with redundant HSPs merged.
:::

:::{grid-item-card} {iconify}`tabler:dna-2` Phobius
:link: /api/generated/rotifer.io.phobius
:link-type: doc

`phobius2table`: signal peptides and transmembrane segments as
domain table rows.
:::

:::{grid-item-card} {iconify}`tabler:stack-2` Architectures
:link: /api/generated/rotifer.interval.architecture
:link-type: doc

`domain2architecture` and `architecture2table`: select non
overlapping regions and write them as domain architectures, and back.
:::

:::{grid-item-card} {iconify}`tabler:tools` Domain annotation
:link: /api/generated/rotifer.tools.domains
:link-type: doc

`acc2pfam`, `acc2profiledb` and the `dom_*` counters: annotate
proteins by accession and summarize their architectures.
:::

:::{grid-item-card} {iconify}`tabler:binary-tree-2` Alignments and tables
:link: /api/generated/rotifer.tools.alignments
:link-type: doc

`aln_order_by_tree` reorders an alignment by a tree; `tdesc`
({py:mod}`rotifer.tools.tables`) numbers the columns of a table.
:::

:::{grid-item-card} {iconify}`tabler:world-download` Database mirrors
:link: /api/generated/rotifer.tools.rrsw
:link-type: doc

`rrsw`: download and unpack mirrors of public databases over FTP or
rsync.
:::

:::{grid-item-card} {iconify}`tabler:folders` Structures and profiles
:link: /api/generated/rotifer.io.dalidat
:link-type: doc

`dat2fa` extracts sequences from DaliLite files; `hhsuite2table`
({py:mod}`rotifer.io.hhsuite`) parses HH-suite results.
:::

:::{grid-item-card} {iconify}`tabler:terminal-2` Command layer
:link: /api/generated/rotifer.cli
:link-type: doc

How commands are built on the library, and the helpers they share.
:::
::::

## Quickstart

Load the module (`module load rotifer`) or add the repository's
`bin/` directory to `PATH`. Every command prints its options with
`-h`.

::::{tab-set}

:::{tab-item} Domain tables

Turn a hmmscan report into a domain table, keeping hits with
E-value up to 0.01:

```bash
hmmscan Pfam-A.hmm proteins.faa > proteins.hmmscan
hmmer2table -y -e 0.01 proteins.hmmscan > proteins.domains.tsv
```

RPS-BLAST reports need `-s`, so that the profile is reported as the
domain and the protein as the sequence:

```bash
blast2table -s -c cdd -y proteins.rpsblast > proteins.cdd.tsv
```
:::

:::{tab-item} Architectures

Choose the best non overlapping domains of each protein and write
them as an architecture string:

```bash
hmmer2table proteins.hmmscan | domain2architecture > proteins.arch
```

```text
ID           architecture                              regions      evalue
NP_414543.1  AA_kinase+ACT_9+ACT_9+NAD_binding_3+...   2..284&...   1.2e-45&...
```

`architecture2table` turns architectures back into one row per
region.
:::

:::{tab-item} Annotate proteins

Annotate proteins with Pfam, signal peptides and transmembrane
segments. Sequences come from the local sequence database first,
then from NCBI; `--fasta` uses a file instead:

```bash
acc2pfam WP_063732599.1 YP_009724395.1
acc2profiledb --fasta proteins.faa
```

Summarize an architecture table:

```bash
dom_count proteins.arch               # regions per domain and protein
dom_per_prot proteins.arch            # regions per protein
distinct_dom_per_prot proteins.arch   # distinct domains per protein
```
:::

:::{tab-item} Mirrors

List the configured mirrors, show what an update would download,
then run it:

```bash
rrsw -l
rrsw -n uniprot
rrsw uniprot
```
:::

:::{tab-item} From Python

The same work without subprocesses. The parsers return DataFrames
that chain directly:

```python
from rotifer.io.hmmer import read_hmmer
from rotifer.interval.architecture import domain2architecture, architecture2table

domains = read_hmmer("proteins.hmmscan", evalue=0.01)
architectures = domain2architecture(domains)
regions = architecture2table(architectures)
```
:::
::::

## Commands

The table is sortable and searchable. "Reads" is the input expected
on the command line or standard input; every command writes a
tab-separated table to standard output unless noted.

```{list-table}
:class: sphinx-datatable
:header-rows: 1

* - Command
  - Reads
  - Writes
  - Library function
* - hmmer2table
  - HMMER 3 reports (hmmscan, hmmsearch), or directories of them
  - one row per domain hit
  - {py:func}`rotifer.io.hmmer.read_hmmer`
* - blast2table
  - BLAST+ pairwise text reports
  - one row per HSP, redundant HSPs merged
  - {py:func}`rotifer.io.blast.read_blast`
* - phobius2table
  - Phobius reports (default format)
  - one row per signal peptide or transmembrane segment
  - {py:func}`rotifer.io.phobius.read_phobius`
* - dat2fa
  - DaliLite .dat files
  - FASTA
  - {py:func}`rotifer.io.dalidat.dat2fasta`
* - hhsuite2table
  - HH-suite .hhr results
  - one row per hit
  - {py:func}`rotifer.io.hhsuite.hhsuite2pandas`
* - domain2architecture
  - domain tables (the output of the commands above)
  - one architecture per sequence
  - {py:func}`rotifer.interval.architecture.domain2architecture`
* - architecture2table
  - architecture tables
  - one row per region
  - {py:func}`rotifer.interval.architecture.architecture2table`
* - acc2pfam
  - protein accessions or FASTA
  - Pfam and Phobius regions, padded for reading
  - {py:func}`rotifer.tools.domains.acc2pfam`
* - acc2profiledb
  - protein accessions or FASTA
  - profile database and Phobius regions, padded for reading
  - {py:func}`rotifer.tools.domains.acc2profiledb`
* - dom_count
  - an architecture table
  - regions of each domain in each protein
  - {py:func}`rotifer.tools.domains.dom_count`
* - dom_per_prot
  - an architecture table
  - regions in each protein
  - {py:func}`rotifer.tools.domains.dom_per_prot`
* - distinct_dom_per_prot
  - an architecture table
  - distinct domains in each protein
  - {py:func}`rotifer.tools.domains.distinct_dom_per_prot`
* - tdesc
  - any delimited table
  - column numbers and names
  - {py:func}`rotifer.tools.tables.tdesc`
* - aln_order_by_tree
  - an alignment and a tree
  - the alignment, in the order of the tree's leaves
  - {py:func}`rotifer.tools.alignments.aln_order_by_tree`
* - rrsw
  - names of targets in the mirror table
  - files under the data directory; progress on standard error
  - {py:func}`rotifer.tools.rrsw.rrsw`
```

## How a command runs

:::{dropdown} From bin/ to the library

The executables in `bin/` are symbolic links to one dispatcher,
`bin/_rotifer_cli`. It looks at the name it was called by and runs
`rotifer.cli.<name>`: `bin/hmmer2table` runs
`rotifer.cli.hmmer2table.main()`. The dispatcher uses the library of
the checkout it lives in, so several checkouts can coexist. The same
command can be run without `bin/` as
`python -m rotifer.cli.hmmer2table`.
:::

:::{dropdown} Inside a command module

Each module in {py:mod}`rotifer.cli` defines `add_arguments(parser)`,
`run(args)` and `main(argv=None)`, using only the standard library's
`argparse`. Command modules import the library lazily and never
import pandas at load time, so a command starts in about a tenth of
a second: `hmmer2table` and `blast2table` run faster than the Perl
tools they replace. To add a command, write the library function
first, then a module in `lib/rotifer/cli/` and a link in `bin/` to
`_rotifer_cli`.
:::

:::{dropdown} Removed commands

The Perl tools without a Python replacement were removed. Their
names in `bin/` now print what to use instead and exit with an
error, so scripts that still call them fail loudly:

```text
$ tjoin a.tsv b.tsv
tjoin was removed from ROTIFER together with its Perl code: use pandas.merge, or join(1).
```

These placeholders will be deleted in a later release; see the
{doc}`changelog </changelog>` for the full list.
:::

:::{dropdown} Output stability

The column order and formatting of each command match the Perl tool
it replaced, byte for byte, on the recorded examples in
`test/data/golden`. Those examples run as tests, so a change that
alters the output of a command fails the test suite.
:::

## Common patterns

Annotate a proteome with Pfam and count domains per protein:

```bash
hmmscan --cpu 8 Pfam-A.hmm proteome.faa > proteome.hmmscan
hmmer2table proteome.hmmscan | domain2architecture > proteome.arch
dom_count proteome.arch > proteome.dom_count.tsv
```

Combine Pfam, CDD and Phobius regions in one architecture:

```bash
{ hmmer2table proteins.hmmscan
  blast2table -s -c cdd proteins.rpsblast
  phobius2table proteins.phobius
} | domain2architecture -e 0.01 > proteins.arch
```

In Python, keeping everything in memory. Tables are matched by column
name (sequence, domain or model, start, end, evalue), so the HMMER and
Phobius tables can be concatenated as they are:

```python
import pandas as pd
from rotifer.io.hmmer import read_hmmer
from rotifer.io.phobius import read_phobius
from rotifer.interval.architecture import domain2architecture

regions = pd.concat([
    read_hmmer("proteins.hmmscan"),
    read_phobius("proteins.phobius"),
])
architectures = domain2architecture(regions, evalue=0.01)
```

Find the columns of an unfamiliar table before filtering it:

```bash
tdesc proteins.domains.tsv
```
