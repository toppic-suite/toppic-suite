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

pacman -S mingw-w64-x86_64-xerces-c

pacman -S mingw-w64-x86_64-xalan-c
```

Thanks to MSYS2, we can install them very easily. Now we have `C:\msys64\mingw64\include` and `C:\msys64\mingw64\lib` for include and linking.

After installing, please add `C:\msys64\mingw64\bin` into your PATH environmental variable. In the following instructions, we will use both the MSYS2 shell and Windows CMD.

**Please make sure you can use `g++` in both MSYS2 shell and Windows CMD.**

# TopPIC suite

**We use Windows CMD to build TopPIC suite.**

```sh
git clone https://github.com/toppic-suite/toppic-suite.git
cd toppic-suite
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```
