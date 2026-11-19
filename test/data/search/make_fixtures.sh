#!/usr/bin/env bash
# How the files in this directory were made (issue #141). Not run by the tests:
# they only read the committed files.
#
#   fixture.fa.gz       Swiss-Prot subset: every sequence that psiblast (3 rounds,
#                       E <= 0.01) finds for SRC_HUMAN (P12931) and ArcA (P0A9Q1)
#                       in the full database, plus 400 random sequences as filler
#   src.fa, arca.fa     the queries
#   <case>.asn          not kept: psiblast -outfmt 11 archive of each case
#   <case>.fmt7.tsv    blast_formatter -outfmt 7 with BLAST_OUTFMT's fields
#                       (rotifer.io.search.BLAST_OUTFMT)
#   <case>.txt.gz       blast_formatter -outfmt 0, the text report blast2table reads
#   <case>.b2t.tsv      bin/blast2table <case>.txt, default options
#
# Cases: src1 (1 round), src3 (up to 3 rounds; SRC converges after 2), arca3 (3 rounds).
# -evalue 10 -inclusion_ethresh 0.01 puts weak hits (E > 0.01) in the report so that
# the inclusion threshold of blast2table (0.01) is crossed between rounds.
#
# Usage: BLASTBIN=/path/to/blast+ SPROT=/path/to/sprot (BLAST db of Swiss-Prot,
#        made with -parse_seqids) SPROT_FASTA=... PERL5=... ./make_fixtures.sh OUTDIR
#        BioPerl must be loadable by the perl on PATH (blast2table needs Bio::SearchIO).
set -euo pipefail

OUT=${1:?usage: make_fixtures.sh OUTDIR}
BLASTBIN=${BLASTBIN:-$HOME/.local/envs/blast/bin}
SPROT=${SPROT:?BLAST database of Swiss-Prot}
SPROT_FASTA=${SPROT_FASTA:?Swiss-Prot FASTA}
ROTIFER=${ROTIFER:-$(cd "$(dirname "$0")/../../.." && pwd)}
FIELDS='qaccver saccver sstart send evalue qstart qend bitscore slen qlen'
export PATH=$BLASTBIN:$PATH
mkdir -p "$OUT"; cd "$OUT"

# Queries
blastdbcmd -db "$SPROT" -entry P12931 -outfmt %f | sed '1s/.*/>P12931/' > src.fa
blastdbcmd -db "$SPROT" -entry P0A9Q1 -outfmt %f | sed '1s/.*/>P0A9Q1/' > arca.fa

# Fixture database: homologs of both queries in the full database + random filler
for q in src arca; do
    psiblast -query $q.fa -db "$SPROT" -num_iterations 3 -evalue 0.01 -num_alignments 300 \
             -num_threads 8 -outfmt 7 -out $q.full.tsv
done
cat src.full.tsv arca.full.tsv | grep -v '^#' | cut -f2 | sort -u > hits.txt
python3 - "$SPROT_FASTA" <<'PY'
import random, sys
ids = [l[1:].split()[0].split('|')[1] for l in open(sys.argv[1]) if l[0] == '>']
random.seed(141)
have = set(open('hits.txt').read().split())
sample = [i for i in random.sample(ids, 400) if i not in have]
open('filler.txt', 'w').write('\n'.join(sample) + '\n')
PY
blastdbcmd -db "$SPROT" -entry_batch hits.txt   -outfmt %f  > fixture.fa
blastdbcmd -db "$SPROT" -entry_batch filler.txt -outfmt %f >> fixture.fa
mkdir -p fdb; makeblastdb -in fixture.fa -dbtype prot -parse_seqids -out fdb/fx > /dev/null

# Searches
run() {   # name query rounds
    psiblast -query $2.fa -db fdb/fx -num_iterations $3 -evalue 10 -inclusion_ethresh 0.01 \
             -num_alignments 1000 -num_threads 4 -outfmt 11 -out $1.asn
    blast_formatter -archive $1.asn -outfmt "7 $FIELDS" -num_alignments 1000 -out $1.fmt7.tsv
    blast_formatter -archive $1.asn -outfmt 0 -num_alignments 1000 -out $1.txt
    # < /dev/null: blast2table otherwise waits for stdin after reading the report
    perl "$ROTIFER/bin/blast2table" $1.txt < /dev/null > $1.b2t.tsv
}
run src1 src 1
run src3 src 3
run arca3 arca 3

gzip -9 -f fixture.fa src1.txt src3.txt arca3.txt
rm -rf fdb *.asn *.full.tsv hits.txt filler.txt
