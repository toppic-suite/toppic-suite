#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
# ../../bin/test_zero_ptm tmp/st.fasta tmp/st_no_digestion.msalign -c C57 | tee tmp/log 

cp ../data/st_no_digestion_java_0.7.table ./tmp
./find_semi_type_ids.py -i  tmp/st_no_digestion_java_0.7.table -n 0 |tee tmp/find_type.log
./compare.py -c ./tmp/st_no_digestion.ZERO_COMPLETE_TABLE -j tmp/st_no_digestion_java_0.7.table.COMPLETE -o tmp/complete_diff
./compare.py -c ./tmp/st_no_digestion.ZERO_PREFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.PREFIX -o tmp/prefix_diff
./compare.py -c ./tmp/st_no_digestion.ZERO_SUFFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.SUFFIX -o tmp/suffix_diff
./compare.py -c ./tmp/st_no_digestion.ZERO_INTERNAL_TABLE -j tmp/st_no_digestion_java_0.7.table.INTERNAL -o tmp/internal_diff
