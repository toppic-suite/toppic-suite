#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
#../../bin/test_one_ptm_filter tmp/st.fasta tmp/st_no_digestion.msalign -c C57 | tee tmp/log 

