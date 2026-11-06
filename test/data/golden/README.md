# Golden data: behaviour of the former Perl tools

This directory freezes the output of ROTIFER's Perl command line tools
(`bin/*`, `perl/lib`) on fixed inputs, recorded **before** they were
replaced by Python. The Python replacements must reproduce these files
byte for byte; see `test/golden.py` for the test helper.

## Layout

| Path | Content |
|---|---|
| `inputs/` | raw reports from external programs, plus intermediate tables produced by the Perl tools |
| `<command>/<case>.cmd` | command line, run from this directory |
| `<command>/<case>.out` | the Perl tool's standard output for that command line |
| `lmod/<function>.cmd` | the pipeline behind each shell function of `etc/lmod/modules/rotifer/0.3.lua`, with `efetch`, `hmmscan`, `rpsblast` and `phobius` replaced by the recorded reports in `inputs/` |

## Inputs

* `set.faa`: 20 *E. coli* K-12 proteins (NC_000913.3), chosen to include
  multi-domain, membrane and secreted proteins; `multi.faa`: the first 6.
* `hmmscan.out`: `hmmscan Pfam-A.hmm set.faa` (HMMER 3.4, Pfam 38.2).
* `hmmsearch.out`: Response_reg, HisKA and HATPase_c models against `set.faa`.
* `rpsblast.out`: `rpsblast -db cdd -evalue 0.01 -max_target_seqs 6` (BLAST+ 2.16.0, CDD 2025-05).
* `psiblast.out`: 3 PSI-BLAST rounds of NP_414555.1 against a 47-protein database.
* `blastp.out`: `multi.faa` against the same database (several HSPs per hit).
* `phobius.out`: Phobius 1.01, default (long) output.
* `dali/*.dat`: DaliLite DAT files from `dalidb/DAT`.
* `aln.fa`, `tree.nwk`: FAMSA 2.2.3 alignment of `multi.faa` and a hand-written tree.

## How it was recorded

Perl 5.32.1 with BioPerl 1.7.8 and the CPAN modules required by `perl/lib`
(conda-forge/bioconda). Each case was run as `bash -c "$(cat <case>.cmd)" </dev/null`
with `bin/` (the Perl tools) first in `PATH`.

## Known quirks preserved on purpose

* `hmmer2table` prints no header unless `-y` is given.
* `blast2table` reports every PSI-BLAST round as iteration 1 and merges HSPs
  repeated across rounds (BioPerl behaviour).
* `domain2architecture -i all` may print empty brackets, e.g. `PRK01212[]`,
  when a region contains others that are attributed to a later container.

The `lmod/dom_*` outputs include a spurious first row (`tgroup` aggregated the
header line as data); the Python commands print a real header instead, so
tests compare only the data rows of those three files.
