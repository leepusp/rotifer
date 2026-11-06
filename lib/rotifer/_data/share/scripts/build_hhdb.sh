#!/bin/bash
# Code by G.G.Nicastro

# Input
MSA_DIR=$1
BASE=$2

# REF: https://github.com/soedinglab/hh-suite/wiki#building-customized-databases
# You need to have a directory with all alignments you want 1 dir up to local.

# indexing alinments
cd ${MSA_DIR}/
ffindex_build -s ../${BASE}_msa.ff{data,index} .
cd ..

# create db
ffindex_apply ${BASE}_msa.ff{data,index} \
  -i ${BASE}_a3m.ffindex -d ${BASE}_a3m.ffdata \
    -- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0
rm ${BASE}_msa.ff{data,index}

# computing hhm
ffindex_apply ${BASE}_a3m.ff{data,index} -i ${BASE}_hhm.ffindex -d ${BASE}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0

# compute contexts to prefiltering
cstranslate -f -x 0.3 -c 4 -I a3m -i ${BASE}_a3m -o ${BASE}_cs219 

# Putting everything together
sort -k3 -n -r ${BASE}_cs219.ffindex | cut -f1 > sorting.dat

ffindex_order sorting.dat ${BASE}_hhm.ff{data,index} ${BASE}_hhm_ordered.ff{data,index}
mv ${BASE}_hhm_ordered.ffindex ${BASE}_hhm.ffindex
mv ${BASE}_hhm_ordered.ffdata ${BASE}_hhm.ffdata

ffindex_order sorting.dat ${BASE}_a3m.ff{data,index} ${BASE}_a3m_ordered.ff{data,index}
mv ${BASE}_a3m_ordered.ffindex ${BASE}_a3m.ffindex
mv ${BASE}_a3m_ordered.ffdata ${BASE}_a3m.ffdata
