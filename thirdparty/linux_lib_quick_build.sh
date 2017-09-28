#!/bin/bash

BOOST_VERSION=1_65_0

XERCES_C_VERSION=3.2.0

XALAN_C_VERSION=1.11

PWIZ_VERSION=3_0_11392

echo "boost version: " ${BOOST_VERSION}
echo "xerces-c version: " ${XERCES_C_VERSION}
echo "xalan-c version: " ${XALAN_C_VERSION}
echo "Proteowizard version:" ${PWIZ_VERSION}

INSTALL_ROOT=$(pwd)

echo "Install into " ${INSTALL_ROOT}

if [ ! -e "pwiz-src-without-v-${PWIZ_VERSION}.tar.bz2" ]; then
  echo "pwiz-src-without-v-${PWIZ_VERSION}.tar.bz2 does not exist"
  exit
fi 

echo "Building boost"

wget https://dl.bintray.com/boostorg/release/1.65.0/source/boost_${BOOST_VERSION}.tar.gz

tar zxvf boost_${BOOST_VERSION}.tar.gz

cd boost_${BOOST_VERSION}

./bootstrap.sh

./b2 --with-filesystem --with-system --with-program_options --with-thread --with-regex --with-chrono --with-iostreams --prefix=${INSTALL_ROOT} toolset=gcc variant=release link=shared threading=multi runtime-link=shared install

echo "Building xerces-c"

cd ${INSTALL_ROOT}

wget http://ftp.wayne.edu/apache//xerces/c/3/sources/xerces-c-${XERCES_C_VERSION}.tar.gz

tar zxvf xerces-c-${XERCES_C_VERSION}.tar.gz

cd xerces-c-${XERCES_C_VERSION}

./configure CFLAGS=-O3 CXXFLAGS=-O3 --prefix=${INSTALL_ROOT} --enable-transcoder-iconv --disable-network

cd src

make install

echo "Buildling xalan-c"

cd ${INSTALL_ROOT}

wget http://apache.cs.utah.edu/xalan/xalan-c/sources/xalan_c-${XALAN_C_VERSION}-src.tar.gz

tar zxvf xalan_c-${XALAN_C_VERSION}-src.tar.gz

cd xalan-c-${XALAN_C_VERSION}/c

export XERCESCROOT=${INSTALL_ROOT}

export XALANCROOT=${INSTALL_ROOT}/xalan-c-${XALAN_C_VERSION}/c/

./runConfigure -p linux -c gcc -x g++ -b 64 -P ${INSTALL_ROOT}

make install

echo "Building Proteowizard"

mkdir -p pwiz

tar xvjf pwiz-src-without-v-${PWIZ_VERSION}.tar.bz2 -C pwiz

cd pwiz

bash quickbuild.sh

echo "All libraries are in "${INSTALL_ROOT}"/lib"
