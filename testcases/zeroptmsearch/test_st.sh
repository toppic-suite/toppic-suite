#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
#../../bin/test_zero_ptm tmp/st.fasta tmp/st_no_digestion.msalign -c C57 | tee tmp/log 

cp ../data/st_no_digestion_java_0.7_result ./tmp
./find_semi_type_ids.py -i  tmp/st_no_digestion_java_0.7_result -n 0
