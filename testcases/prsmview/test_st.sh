#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
../../bin/test_prsm_view tmp/st.fasta tmp/st_no_digestion.msalign -c C57 -d -t FDR| tee tmp/log 
