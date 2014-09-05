#!/bin/sh

mkdir tmp
cd tmp
ln -s  ../../data/st.fasta ./ 
ln -s  ../../data/st_no_digestion.msalign ./
cd ..
../../bin/test_evalue tmp/st.fasta tmp/st_no_digestion.msalign -c C57 -d -t FDR| tee tmp/log 
cp ../data/st_no_digestion_java_0.7.table ./tmp
awk -F '\t' '{if ($20 < 0.01) print $4}' tmp/st_no_digestion_java_0.7.table  > tmp/java_fdr_0.01
sort -n tmp/java_fdr_0.01 > tmp/java_fdr_0.01_sorted
awk -F '\t' '{print $5}' tmp/st_no_digestion.OUTPUT_TABLE  > tmp/c_fdr_0.01
diff -y tmp/java_fdr_0.01_sorted tmp/c_fdr_0.01
