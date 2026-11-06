tmp=$(mktemp); treeutil -of leaves inputs/tree.nwk > $tmp; aln2seqrows -r "\t" inputs/aln.fa | tjoin -r "\n" -i1 1 -i2 1 -1 1 -f1 '$F[1] ne "ID"' -f '$F[0]=">$F[0]";1' -p $tmp -; rm -f $tmp
