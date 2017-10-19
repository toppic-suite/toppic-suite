#!/bin/sh
if [ "$#" -ne 1 ];then
  echo "Usage: $0 version_number" >&2
  exit 1
fi

rm -rf toppic-linux-${1}
mkdir toppic-linux-${1}
mkdir toppic-linux-${1}/example_files
cp ../bin/toppic toppic-linux-${1} 
cp ../LICENSE toppic-linux-${1}
cp ./example_files/common_mods.txt toppic-linux-${1}/example_files
cp ./example_files/fixed_mod.txt toppic-linux-${1}/example_files
cp -r ../toppic_resources toppic-linux-${1}

zip -r toppic-linux-${1}.zip toppic-linux-${1}

