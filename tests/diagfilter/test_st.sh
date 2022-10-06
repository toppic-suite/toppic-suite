#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
#../../bin/test_diag_filter tmp/st.fasta tmp/st_no_digestion.msalign -c C57 | tee tmp/log 

cp ../data/st_no_digestion_java_0.7.table ./tmp
./find_ptm_ids.py -i tmp/st_no_digestion_java_0.7.table
./compare.py -c ./tmp/st_no_digestion.FILTER_TABLE -j tmp/st_no_digestion_java_0.7.table.PTM -o tmp/filter_diff
