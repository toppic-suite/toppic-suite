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
cp ../bin/topindex toppic-linux-${1} 
cp ../bin/topfd toppic-linux-${1} 
cp ../bin/topmg toppic-linux-${1} 
cp ../bin/topdiff toppic-linux-${1} 
cp ../bin/topconvert toppic-linux-${1} 
cp ../LICENSE toppic-linux-${1}

cp ./example_files/common_mods.txt toppic-linux-${1}/example_files
cp ./example_files/fixed_mod.txt toppic-linux-${1}/example_files
cp ../testcases/data/mzxml_test.fasta toppic-linux-${1}/example_files
cp ../testcases/data/mzxml_test.mzXML toppic-linux-${1}/example_files

cp -r ../toppic_resources toppic-linux-${1}
rm ./toppic-linux-${1}/toppic_resources/envcnn_models/envcnn_5_block_model.json

mkdir toppic-linux-${1}/toppic_resources/lib

cp /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_system.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_thread.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_serialization.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_iostreams.so.1.71.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libboost_chrono.so.1.71.0 toppic-linux-${1}

cp /usr/lib/x86_64-linux-gnu/libxerces-c-3.2.so toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libicuuc.so.66 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libicudata.so.66 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libcurl-gnutls.so.4 toppic-linux-${1}

cp /usr/lib/x86_64-linux-gnu/libnghttp2.so.14 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/librtmp.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libssh.so.4 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libpsl.so.5 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libgssapi_krb5.so.2 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libldap_r-2.4.so.2 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libbrotlidec.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/liblber-2.4.so.2 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libcrypto.so.1.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libkrb5.so.3 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libk5crypto.so.3 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libkrb5support.so.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libsasl2.so.2 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libgssapi.so.3 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libbrotlicommon.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libkeyutils.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libheimntlm.so.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libkrb5.so.26 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libasn1.so.8 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libhcrypto.so.4 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libroken.so.18 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libwind.so.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libheimbase.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libhx509.so.5 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libsqlite3.so.0 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libheimbase.so.1 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libnettle.so.7 toppic-linux-${1}
cp /usr/lib/x86_64-linux-gnu/libhogweed.so.5 toppic-linux-${1}

zip -r toppic-linux-${1}.zip toppic-linux-${1}
