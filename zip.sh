#!/bin/sh
if [ "$#" -ne 1 ];then
  echo "Usage: $0 version_number" >&2
  exit 1
fi

cp ./bin/toppic .
zip -r toppic-$1-linux.zip ./toppic_resources/* toppic LICENSE

