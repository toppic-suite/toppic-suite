#!/bin/sh
if [ "$#" -ne 1 ];then
  echo "Usage: $0 version_number" >&2
  exit 1
fi

mkdir toppic-linux-${1}
cp ../bin/toppic toppic-linux-${1} 
cp ../LICENSE toppic-linux-${1}
cp ../var_mods.txt toppic-linux-${1}
cp ../fix_mods.txt toppic-linux-${1}
cp -r ../toppic_resources toppic-linux-${1}

zip -r toppic-linux-${1}.zip toppic-linux-${1}

