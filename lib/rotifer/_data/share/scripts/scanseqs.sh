#!/bin/bash
# File: seqscan

# Basename
target=$1
input=$2

# Cluster sequences using hashclust and 100% coverage + 100% identity
# Table mapping all sequences to non-redundant representatives is mergednrdb.tsv
# Non-redundant FASTA file is mergednrdb_rep.fasta
mmseqs createdb ${input} ${target}
mmseqs clusthash ${target} ${target}nr --threads 20 --min-seq-id 1
mmseqs clust ${target} ${target}nr ${target}.c100i100 --threads 20 --cluster-mode 1 --max-iterations 1000000
mmseqs createtsv ${target} ${target} ${target}.c100i100 ${target}.c100i100.tsv
mmseqs createsubdb ${target}.c100i100 ${target} ${target}.c100i100.seqs
mmseqs convert2fasta ${target}.c100i100.seqs ${target}.c100i100.fa
mmseqs easy-cluster ${target}.c100i100.fa ${target}.c80i70 tmp --min-seq-id 0.7 -c 0.8 --threads 36 > ${target}.c80i70.log 2>&1
ln -s ${target}.c80i70_cluster.tsv ${target}.c80i70.tsv
ln -s ${target}.c80i70_rep_seq.fasta ${target}.c80i70.fa
mmseqs easy-cluster ${target}.c80i70.fa ${target}.c80e3 tmp --threads 36 > ${target}.c80e3.log 2>&1
ln -s ${target}.c80e3_cluster.tsv ${target}.c80e3.tsv
rm -fr tmp

# Detecting sequence signatures in the new proteins
input=${target}.c100i100.fa

# Phobius
cut -f1 -d " " ${input} | parallel -N1 -j10 --pipe --recstart ">" phobius > ${target}.phobius.out 2> ${target}.phobius.err
phobius2table -e 0.0101 ${target}.phobius.out > ${target}.phobius.tsv

# Pfam: hmmscan
cut -f1 -d " " ${input} \
| parallel -N1 -j10 --pipe --recstart ">" hmmscan /databases/pfam/Pfam-A.hmm - \
> ${target}.pfam.hmmscan.out \
2> ${target}.pfam.hmmscan.err
hmmer2table ${target}.pfam.hmmscan.out > ${target}.pfam.hmmscan.tsv
domain2architecture -e 0.0101 ${target}.phobius.tsv ${target}.pfam.hmmscan.tsv > ${target}.pfam.hmmscan.arch
architecture2table ${target}.pfam.hmmscan.arch > ${target}.pfam.hmmscan.arch.tsv
ln -s ${target}.pfam.hmmscan.arch ${target}.pfam.scan.arch
ln -s ${target}.pfam.hmmscan.arch.tsv ${target}.pfam.scan.arch.tsv

# Aravind: hmmscan
cut -f1 -d " " ${input} \
| parallel --pipe -N1 -j36 --recstart '>' hmmscan --cpu 1 /databases/profiledb/hmmdb/aravindDB.hmm - \
> ${target}.aravind.hmmscan.out \
2> ${target}.aravind.hmmscan.err
hmmer2table -c model=version ${target}.aravind.hmmscan.out > ${target}.aravind.hmmscan.tsv

# Aravind: rpsblast
cut -f1 -d " " ${input} \
| parallel -N1 -j36 --pipe --recstart '>' rpsblast -db /databases/profiledb/rpsdb/allprofiles \
> ${target}.aravind.rpsblast.out \
2> ${target}.aravind.rpsblast.err
cat ${target}.aravind.rpsblast.out \
| parallel -N1 -j36 --pipe --recstart RPSBLAST blast2table -s -c profiledb \
> ${target}.aravind.rpsblast.tsv

# Aravind: domain2architecture
(cat ${target}.phobius.tsv ${target}.aravind.hmmscan.tsv; awk '$5<=0.1{print}' ${target}.aravind.rpsblast.tsv) \
| domain2architecture -e 0.0101 > ${target}.aravind.scan.arch
architecture2table ${target}.aravind.scan.arch > ${target}.aravind.scan.arch.tsv

# Main output files
#ls *.tsv *arch
