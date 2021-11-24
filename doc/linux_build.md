For most Linux distribution, you should be able to install the dependency through the package manager. This document describes the configuration we used to build all the required libraries from source.

# Boost

We don't have enough resource to test the minimum version for this project, a recent version of boost should work well. `boost_system`, `boost_filesystem`, `boost_program_options`, `boost_thread`, `boost_chrono`, `boost_iostreams` are required.

```sh
wget https://dl.bintray.com/boostorg/release/1.65.0/source/boost_1_65_0.tar.gz

tar zxvf boost_1_65_0.tar.gz

cd boost_1_65_0

./bootstrap.sh

./b2 --with-filesystem --with-system --with-program_options --with-thread --with-chrono --with-iostreams --prefix=${INSTALL_DIR} toolset=gcc variant=release link=shared threading=multi runtime-link=shared install
```

# Xerces-C++

[Xerces-C++](https://xerces.apache.org/xerces-c/) is used to process XML files.

```sh
wget hhttp://ftp.wayne.edu/apache//xerces/c/3/sources/xerces-c-3.2.0.tar.gz

tar zxvf xerces-c-3.2.0.tar.gz

cd xerces-c-3.2.0

./configure CFLAGS=-O3 CXXFLAGS=-O3 --prefix=${INSTALL_DIR} --enable-transcoder-iconv --disable-network

cd src

make install
```

# Xalan-C++

[Xalan-C++](https://xml.apache.org/xalan-c/) is used for transforming XML documents into HTML.

```sh
wget http://apache.cs.utah.edu/xalan/xalan-c/sources/xalan_c-1.11-src.tar.gz

tar zxvf xalan_c-1.11-src.tar.gz

cd xalan-c-1.11/c

export XERCESCROOT=${INSTALL_DIR}

export XALANCROOT=<path>/xalan-c-1.11/c/

./runConfigure -p linux -c gcc -x g++ -b 64 -P ${INSTALL_DIR}

make install
```

# ProteoWizard

[ProteoWizard](http://proteowizard.sourceforge.net) is used to handle different file formats in proteomics data analysis.

```sh
mkdir -p pwiz

tar xvjf pwiz-src-without-v-3_0_11392.tar.bz2 -C pwiz

cd pwiz

bash quickbuild.sh
```
