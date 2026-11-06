# Changelog

All notable changes to ROTIFER are documented on this page. The
format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/):
the newest release comes first, and the changes of each release are
grouped as **Added**, **Changed**, **Deprecated**, **Removed** and
**Fixed**. Changes that can break existing code or scripts are listed
again, with what to do, under **Upgrade notes**. Click a
version to fold or unfold its notes; the newest is unfolded.

:::{dropdown} v.x.x {bdg-warning}`RELEASE TEMPLATE TO BE DEFINED`
:open:
:name: v-x-x

```{rubric} Highlights
```

- **ROTIFER no longer contains Perl.** The Perl tools still in use
  were rewritten in Python, under the same names and options and with
  the same output; the rest were removed. Every rewritten tool is
  both a library function and a shell command, described in
  {doc}`/cli/index` ({pr}`138`).
- **ROTIFER works outside a git checkout.** Installed copies (pip,
  conda, containers, shared servers) now find their configuration and
  data files ({pr}`130`).
- **Log messages go to standard error, once.** Warnings no longer end
  up in the tables that command line tools write to standard output
  ({pr}`134`).

```{rubric} Upgrade notes
```

- The Perl library (`perl/`, the `Rotifer::` modules) and the
  `PERL5LIB` setting of the Lmod module are gone. Perl scripts outside
  this repository that use `Rotifer::` modules stop working ({pr}`138`).
- 51 command line tools were removed (listed under **Removed**).
  Calling one prints what to use instead and exits with status 1
  ({pr}`138`).
- `gi2context`, `runprofiles` and `rbashpipe` were removed without a
  replacement. They expected GI numbers and an SGE cluster; use
  `bin/scanseqs`, or HMMER and BLAST directly ({pr}`138`).
- `acc2pfam`, `acc2profiledb`, `dom_count`, `dom_per_prot`,
  `distinct_dom_per_prot`, `tdesc` and `aln_order_by_tree` are now
  commands instead of shell functions defined by the Lmod module.
  Shell startup files that sourced the copies in `etc/profile.d` must
  stop doing so ({pr}`138`).
- `dom_count`, `dom_per_prot` and `distinct_dom_per_prot` print a
  real header line. Scripts that skipped the old, malformed first line
  should check it ({pr}`138`).
- `rrsw` messages are timestamped log lines on standard error.
  Scripts that parse the output of `rrsw` (such as `update_uniprot`
  logs) should be checked ({pr}`138`).
- `rotifer.core.loadpath` was removed: use `shutil.which` to find
  executables, or `rotifer.config['base']` for the installation
  directory ({pr}`131`).
- ROTIFER's log handler writes to standard error instead of standard
  output. Code that captured warnings from standard output must read
  standard error ({pr}`134`).

```{rubric} Added
```

- Python versions of `hmmer2table`, `blast2table`, `phobius2table`,
  `dat2fa`, `domain2architecture`, `architecture2table` and `rrsw`,
  each as a library function and a command of the same name:
  {py:mod}`rotifer.io.hmmer`, {py:mod}`rotifer.io.blast`,
  {py:mod}`rotifer.io.phobius`, {py:mod}`rotifer.io.dalidat`,
  {py:mod}`rotifer.interval.architecture` and
  {py:mod}`rotifer.tools.rrsw` ({pr}`138`).
- {py:mod}`rotifer.tools.domains`, {py:mod}`rotifer.tools.tables` and
  {py:mod}`rotifer.tools.alignments`: the former Lmod shell functions
  as library functions, with commands of the same names in `bin/`
  ({pr}`138`).
- {py:mod}`rotifer.cli`, the command layer: every command in `bin/`
  is a link to `bin/_rotifer_cli`, which runs
  `rotifer.cli.<command>`. Commands also run as
  `python -m rotifer.cli.<command>` ({pr}`138`).
- `domain2architecture()` accepts DataFrames, so
  `read_hmmer()`, `read_phobius()` and `read_blast()` chain into it
  without files or subprocesses;
  {py:func}`rotifer.interval.architecture.select_regions` returns the
  selected regions with their coordinates ({pr}`138`).
- {py:func}`rotifer.io.open_input` and
  {py:func}`rotifer.io.each_source`: read files, standard input, gzip
  and bzip2 the same way in every parser ({pr}`138`).
- Recorded output of the Perl tools on fixed inputs
  (`test/data/golden`), run as tests against the Python commands
  ({pr}`138`).
- `ROTIFER_BASE`, `ROTIFER_USER` and `ROTIFER_CACHE` environment
  variables to choose the installation, user and cache directories
  ({pr}`130`).
- `rotifer.io` is a package, and the HH-suite parser
  (`hhsuite2pandas`, `parse_result`, `read_result`,
  `HHRFormatError`) lives in {py:mod}`rotifer.io.hhsuite`
  ({pr}`131`).
- This changelog and the {doc}`/cli/index` page ({pr}`138`).

```{rubric} Changed
```

- `hmmer2table` and `blast2table` are faster: about 0.1 s per call,
  against 0.2 s and 0.5 s for the Perl versions; a 2 MB RPS-BLAST
  report takes 0.17 s instead of 2.9 s ({pr}`138`).
- `acc2pfam` and `acc2profiledb` fetch sequences with
  `rotifer.db.proteins` (the local sequence database first, then
  NCBI) instead of EDirect's `efetch`; `--fasta` skips fetching
  ({pr}`138`).
- `hmmer2table` reads standard input only when no file is given (the
  Perl version also read it whenever it was not a terminal) ({pr}`138`).
- The Lmod module (`rotifer/0.3`) no longer defines shell functions or
  sets `PERL5LIB`, and prints a short notice when loaded ({pr}`138`).
- The `rrsw` configuration patterns are Python regular expressions,
  and `eggNOG.yml` no longer uses a Perl-only YAML tag ({pr}`138`).
- Functions in `rotifer.devel` that ran the Perl tools in a
  subprocess call the library instead (`gian_func.split_by_model`,
  `gian_func.hmm_models_vs_sequence_db`, `gian_func.psiblast`,
  `malu.extract_by_hmm`, and `psiblast` in `beta.blast` and
  `alpha.rodolfo`) ({pr}`138`).
- Configuration and data files ship as package data under
  `lib/rotifer/_data`; `etc/rotifer` and `share/rotifer` are links to
  it ({pr}`130`).
- `pyyaml` is a declared dependency ({pr}`130`).
- `bin/hhsuite2table` imports its parser from
  {py:mod}`rotifer.io.hhsuite` ({pr}`131`).
- ROTIFER's log handler is attached once, to the `rotifer` logger, and
  writes to standard error ({pr}`134`).

```{rubric} Deprecated
```

- The names of the removed tools in `bin/` (links to
  `bin/_rotifer_removed`) will be deleted in a later release
  ({pr}`138`).

```{rubric} Removed
```

- The Perl library (`perl/`, linked as `lib/rotifer/site_perl`) and
  the Perl dependencies in `README.md` and `requirements.txt`
  ({pr}`138`).
- Tools without a remaining use, now links to `bin/_rotifer_removed`:
  `accession2gi`, `aln2seqrows`, `alntool`, `architecture2svg`,
  `blastclust2table`, `blocks2table`, `boundtable2table`,
  `cleanblastclust2table`, `context2network`, `dali`,
  `dali2seqrows`, `fasta2description`, `fasta2id`, `ft2seq`,
  `genbank2table`, `getipgs`, `gi2pp`, `histogram`, `mblastclust`,
  `neighborhood2svg`, `neighborhood_explorer`, `nt2aa`,
  `ontology2hierarchy`, `operon2arch`, `operon2table`,
  `operon_uniquize`, `padtable`, `pdb2seq`, `protein2dna`,
  `rcomplement`, `rconvertseq`, `rotifer_perl`, `runionFind`,
  `runprofiles`, `seq2taxonomy`, `seqsplit`, `seqstats`,
  `table2blocks`, `taxid2taxonomy`, `taxutil`, `tfilter`, `tgroup`,
  `tjoin`, `tpivot`, `tree2distance`, `treeutil`, `ttranspose`,
  `tungroup`, `tupdate`, `txt2tsv` and `unpadtable` ({pr}`138`).
- The `gi2context` pipeline, `bin/rbashpipe` and `runprofiles`, with
  their configuration files ({pr}`138`).
- Scripts that only wrapped Perl tools: `bin/rsearch`,
  `bin/rsearch2`, `bin/rsearch_old`, `bin/phmmscan`, `bin/rgcontext`
  (superseded by `bin/scanseqs`), `rotifer.tools.search`, and an old
  copy of `scanseqs` in `share/scripts` ({pr}`138`).
- The shell function copies in `etc/profile.d` ({pr}`138`).
- `rrsw`: the `kegg` target (KEGG closed its FTP site) and the
  `unpack_kegg`, `unpack_blast`, `dump_blast` and `index_fasta`
  processors, which no target used; callbacks written in Perl
  ({pr}`138`).
- Options that nothing used: Perl code in `hmmer2table -c` and in
  `blast2table -c` and `-r`; `blast2table --mode m8` and
  `--mode orthomcl`, input formats other than BLAST+ text, `-i`,
  `--type` and `-sq` ({pr}`138`).
- `rotifer.core.loadpath` ({pr}`131`).

```{rubric} Fixed
```

- `process_blast_dbs.sh` replaces dashes in `db/FASTA` sequences by
  `X`: its Perl one-liner never ran, because it lost its quotes when
  evaluated ({pr}`138`).
- `gian_func.psiblast` ran `blast2table` from a hard-coded path that
  does not exist on the cluster ({pr}`138`).
- Configuration and data lookups came back empty in installed copies,
  containers and shared servers ({issue}`99`, {issue}`101`,
  {issue}`119`, {issue}`120`, {pr}`130`).
- `~/.rotifer/lib/python` is added to the module search path only if
  it exists, and ROTIFER falls back to a private temporary directory
  when the user directory is not writable ({pr}`130`).
- `rotifer.io.hhsuite` found its parser by cutting its own path at
  the first "rotifer", which failed in checkouts under another name
  and in installed copies, and loaded another checkout's script in
  worktrees ({pr}`131`).
- Log messages were printed twice, and to standard output
  ({pr}`134`).
- `loadConfig` raised `TypeError` or `NameError` on unparsable YAML
  files; it now logs the error and the file ({pr}`134`).
- `rotifer.genome.data` passed an invalid `file=` argument to the
  logger ({pr}`134`).
:::
