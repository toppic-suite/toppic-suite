#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
#../../bin/test_ptm_search tmp/st.fasta tmp/st_no_digestion.msalign -c C57 | tee tmp/log 

cp ../data/st_no_digestion_java_0.7.table ./tmp
./find_semi_type_ids.py -i  tmp/st_no_digestion_java_0.7.table -n 1 |tee tmp/find_type.log
./find_semi_type_ids.py -i  tmp/st_no_digestion_java_0.7.table -n 2 |tee tmp/find_type.log
./compare.py -c ./tmp/st_no_digestion.PTM_1_COMPLETE_TABLE -j tmp/st_no_digestion_java_0.7.table.1_COMPLETE -o tmp/complete_1_diff
./compare.py -c ./tmp/st_no_digestion.PTM_1_PREFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.1_PREFIX -o tmp/prefix_1_diff
./compare.py -c ./tmp/st_no_digestion.PTM_1_SUFFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.1_SUFFIX -o tmp/suffix_1_diff
./compare.py -c ./tmp/st_no_digestion.PTM_1_INTERNAL_TABLE -j tmp/st_no_digestion_java_0.7.table.1_INTERNAL -o tmp/internal_1_diff
./compare.py -c ./tmp/st_no_digestion.PTM_2_COMPLETE_TABLE -j tmp/st_no_digestion_java_0.7.table.2_COMPLETE -o tmp/complete_2_diff
./compare.py -c ./tmp/st_no_digestion.PTM_2_PREFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.2_PREFIX -o tmp/prefix_2_diff
./compare.py -c ./tmp/st_no_digestion.PTM_2_SUFFIX_TABLE -j tmp/st_no_digestion_java_0.7.table.2_SUFFIX -o tmp/suffix_2_diff
./compare.py -c ./tmp/st_no_digestion.PTM_2_INTERNAL_TABLE -j tmp/st_no_digestion_java_0.7.table.2_INTERNAL -o tmp/internal_2_diff
