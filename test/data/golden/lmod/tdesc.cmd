head -n1 inputs/table.tsv | ttranspose -s "\t" | tfilter -f 'unshift(@F,$H{c}++);1'
