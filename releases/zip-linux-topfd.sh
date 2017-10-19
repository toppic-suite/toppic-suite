#!/bin/sh
if [ "$#" -ne 1 ];then
  echo "Usage: $0 version_number" >&2
  exit 1
fi

rm -rf topfd-linux-${1}
mkdir topfd-linux-${1}
mkdir topfd-linux-${1}/example_files
cp ../bin/topfd topfd-linux-${1} 
cp ../LICENSE topfd-linux-${1}
cp ./example_files/common_mods.txt topfd-linux-${1}/example_files
cp ./example_files/fixed_mod.txt topfd-linux-${1}/example_files
cp -r ../toppic_resources topfd-linux-${1}

mkdir topfd-linux-${1}/toppic_resources/lib

cp /usr/lib/x86_64-linux-gnu/libxerces-c-3.1.so topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/libpwiz.so.3 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_system.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_thread.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_serialization.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_iostreams.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libboost_chrono.so.1.58.0 topfd-linux-${1}/toppic_resources/lib/
cp /lib/x86_64-linux-gnu/libbz2.so.1.0 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libicudata.so.55 topfd-linux-${1}/toppic_resources/lib/
cp /usr/lib/x86_64-linux-gnu/libicuuc.so.55 topfd-linux-${1}/toppic_resources/lib/

zip -r topfd-linux-${1}.zip topfd-linux-${1}

