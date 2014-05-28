Windows:
1?	download mingw64 with x86_64-w64-mingw32 and msys64 
2?	set Environment Variablse add <mingw64’s path>/bin to PATH

Xcerces
3?	download xerces-c-3.1.1
4?	unzip xerces to msys64 and get into the path xerces-c-3.1.1
5?	run ./configure --disable-sse2 –host=x86_64-w64-mingw64 –build=x86_64-w64-mingw64 CFLAGS=-O3 CXXFLAGS=-O3 --prefix=<the path you want to save the include and lib of xerces>
6?	enter src and run “make clean” “make”
7?	copy the result to mingw64 

Linux
./configure --disable-sse2

Boost: 
1.  download boost 1.55.0
2.  cd boost_1.55.0/tools/build/v2
3.  bootstrap.sh mingw
4.  cd boost_1.55.0
5.  tools/build/v2/b2.exe --prefix=<destination path> --toolset=gcc target-os=windows link=static runtime-link=static architecture=x86 address-model=64 threading=single

Xalan
  export XERCESCROOT=~/git_cpp/proteomics_cpp/thirdparty/xerces-c-3.1.1/
  export XALANCROOT=~/git_cpp/proteomics_cpp/thirdparty/xalan-c-1.11/c/
  ./runConfigure -p mingw-msys -c gcc -x g++ -r none -z -march=athlon64 -z -m64
  edit src/xalanc/Makefile: all extra line for TestXPath 
  edit Makefile.incl and remove "/" in prefix
  see the Makefile in the directory
  make 
  make install

8?	copy src of proteomics_cpp and Release into msys and enter Release and run “make clean” “make” and it will build proteomics_cpp.exe and copy it to /MsAlignPlus/WebContent/tool
9?	download tomcat and create MsAlignPlus folder in \webapps and copy the files in /MsAlignPlus/WebContent/ to \webapps\ MsAlignPlus
10?	run tomcat and access it with address http://localhost:8080/MsAlignPlus/

Download address:
mingw64: http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-posix/seh/
msys64: http://sourceforge.net/projects/msys2/files/Base/x86_64/
xerces-c-3.1.1 : http://xerces.apache.org/xerces-c/download.cgi
tomcat: http://tomcat.apache.org/download-70.cgi
