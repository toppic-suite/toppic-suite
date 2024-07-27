This document describes how to build the command line and GUI tools in TopPIC suite.

# MSYS2

[MSYS2](http://www.msys2.org) provides an excellent building platform for Windows. Follow the instructions on the MSYS2 website for installing MSYS2 and updating the package database and core system packages. 

# Install required packages and download the source code using an MSYS2 shell

Open an MSYS2 shell and use the following commands for installing packages and
downloading the source code.

```sh
pacman -S git

pacman -S mingw-w64-x86_64-gcc

pacman -S mingw-w64-x86_64-make

pacman -S mingw-w64-x86_64-cmake

pacman -S mingw-w64-x86_64-boost

pacman -S mingw-w64-x86_64-xerces-c

pacman -S mingw-w64-x86_64-qt5-base

git clone https://github.com/toppic-suite/toppic-suite.git
```

# Use Windows Terminal to build TopPIC suite.

Add `C:\msys64\mingw64\bin` into your PATH environmental variable. 

**Make sure you can use `g++` in Windows Terminal.**

Open a Windows Terminal

```sh
cd C:\msys64\home\$your_user_name
cd toppic-suite
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

**Move  the folder toppic_resources to the folder bin**
```sh
cd ..\bin
move ..\resources .
```
