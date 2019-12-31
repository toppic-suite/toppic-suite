#!/bin/sh
set -x
if [ "$#" -ne 1 ];then
  echo "Usage: $0 version_number" >&2
  exit 1
fi

rm -rf toppic-linux-${1}
mkdir toppic-linux-${1}
mkdir toppic-linux-${1}/example_files
cp ../bin/toppic toppic-linux-${1}
cp ../bin/topfd toppic-linux-${1} 
cp ../bin/topmg toppic-linux-${1} 
cp ../bin/topmerge toppic-linux-${1} 
cp ../LICENSE toppic-linux-${1}

cp ./example_files/common_mods.txt toppic-linux-${1}/example_files
cp ./example_files/fixed_mod.txt toppic-linux-${1}/example_files
cp ../testcases/data/mzxml_test.fasta toppic-linux-${1}/example_files
cp ../testcases/data/mzxml_test.mzXML toppic-linux-${1}/example_files

cp -r ../toppic_resources toppic-linux-${1}

mkdir toppic-linux-${1}/toppic_resources/lib

cp /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_thread.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_serialization.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_iostreams.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_chrono.so.1.65.1 toppic-linux-${1}/toppic_resources/lib/

cp /usr/lib/x86_64-linux-gnu/libxerces-c-3.2.so toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libicuuc.so.60 toppic-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libicudata.so.60 toppic-linux-${1}/toppic_resources/lib/
cp /lib/x86_64-linux-gnu/libbz2.so.1.0 toppic-linux-${1}/toppic_resources/lib/

zip -r toppic-linux-${1}.zip toppic-linux-${1}
