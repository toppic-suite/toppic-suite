This document describes how to build the command line and GUI tools in TopPIC suite.

# MSYS2

[MSYS2](http://www.msys2.org) provides an excellent building platform for Windows. Please follow the instructions on MSYS2 website for installation and update the package database and core system packages. Then we can install other required packages using the commands below in MSYS2 shell:

```sh
pacman -S git

pacman -S mingw-w64-x86_64-gcc

pacman -S mingw-w64-x86_64-make

pacman -S mingw-w64-x86_64-cmake

pacman -S mingw-w64-x86_64-boost

pacman -S mingw-w64-x86_64-qt5
```

Thanks to MSYS2, we can install them very easily. Now we have `C:\msys64\mingw64\include` and `C:\msys64\mingw64\lib` for include and linking.

After installing, please add `C:\msys64\mingw64\bin` into your PATH environmental variable. In the following instructions, we will use both the MSYS2 shell and Windows CMD.

# Xerces-C++

[Xerces-C++](https://xerces.apache.org/xerces-c/) is used to process XML files.

**We use Windows CMD to build the Xerces-c.**

```sh
cd xerces-c-3.2.0

mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=C:/msys64/mingw64 -G "MinGW Makefiles" ..

mingw32-make install
```

# Xalan-C++

[Xalan-C++](https://xml.apache.org/xalan-c/) is used for transforming XML documents into HTML.

**We use MSYS2 shell to build the Xerces-c.**

```sh
cd xalan-c-1.11/c

export XALANCROOT=<path>/xalan-c-1.11/c

./runConfigure -p mingw-msys -c gcc -x g++ -r none -b 64 

mingw32-make
```

# ProteoWizard

[ProteoWizard](http://proteowizard.sourceforge.net) is used to handle different file formats in proteomics data analysis.

**We use Windows CMD to build the Xerces-c.**

Since we only use a small part of the ProteoWizard library, we use our own Makefile to build a patched version. Please see the [Makefile.pwiz](Makefile.pwiz) and [pwiz.patch](pwiz.patch) in this folder.

# TopPIC suite

**We use Windows CMD to build the Xerces-c.**

```sh
git clone https://github.com/toppic-suite/toppic-suite.git
cd toppic-suite
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```
