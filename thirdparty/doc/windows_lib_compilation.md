##Windows:
* download mingw64 with x86_64-w64-mingw32 and msys64 
* set Environment Variablse add <mingw64 path>/bin to PATH

##Boost: 
* download boost 1.56.0
* cd boost_1.56.0/tools/build/v2
* bootstrap.sh mingw
* cd boost_1.56.0
* tools/build/v2/b2.exe --prefix=<destination path> --toolset=gcc
target-os=windows link=static runtime-link=static architecture=x86
address-model=64 threading=single (threading=single maybe removed)

##Xcerces
* download xerces-c-3.1.1
* unzip xerces to msys64 and get into the path xerces-c-3.1.1
```sh
./configure host=x86_64-w64-mingw64 build=x86_64-w64-mingw64 CFLAGS=-O3 
```
CXXFLAGS=-O3 --prefix=<the path you want to save the include and lib of xerces>
--disable-sse2 (--disable-sse2 maybe removed)
* enter src and run make clean make

##Xalan
* export XERCESCROOT=<path>/proteomics_cpp/thirdparty/xerces-c-3.1.1/
* export XALANCROOT=<path>/proteomics_cpp/thirdparty/xalan-c-1.11/c/
* ./runConfigure -p mingw-msys -c gcc -x g++ -r none -z -march=athlon64 -z -m64
* edit src/xalanc/Makefile 
* edit src/xalanc/Utils/XalanMsgLig/Makefile
* edit Makefile.incl and remove "/" in prefix
   see the Makefile in the directory
* make make install

##Download links:


| Software | Version   |
|----------|-----------|
|[mingw64](http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-posix/seh/x86_64-4.8.2-release-posix-seh-rt_v3-rev4.7z/download) | x86_64-4.8.2 |
|[msys64](http://sourceforge.net/projects/msys2/files/Base/x86_64/msys2-x86_64-20140910.exe/download)| 20140910   |
|[boost](http://sourceforge.net/projects/boost/files/boost/1.56.0/)| 1.56.0|
|[xerces-c](http://xerces.apache.org/xerces-c/download.cgi)|3.1.1|
|[xalan](http://xalan.apache.org/xalan-c/download.html)|1.11|
