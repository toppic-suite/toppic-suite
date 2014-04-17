Windows:
1?	download mingw64 with x86_64-w64-mingw32 and msys64 
2?	set Environment Variablse add <mingw64’s path>/bin to PATH
3?	download xerces-c-3.1.1
4?	unzip xerces to msys64 and get into the path xerces-c-3.1.1
5?	run ./configure –host=x86_64-w64-mingw32 –build=x86_64-w64-mingw32 --prefix=<the path you want to save the include and lib of xerces>
6?	enter src and run “make clean” “make”
7?	copy the result to mingw64 
8?	copy src of proteomics_cpp and Release into msys and enter Release and run “make clean” “make” and it will build proteomics_cpp.exe and copy it to /MsAlignPlus/WebContent/tool
9?	download tomcat and create MsAlignPlus folder in \webapps and copy the files in /MsAlignPlus/WebContent/ to \webapps\ MsAlignPlus
10?	run tomcat and access it with address http://localhost:8080/MsAlignPlus/

Download address:
mingw64: http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-posix/seh/
msys64: http://sourceforge.net/projects/msys2/files/Base/x86_64/
xerces-c-3.1.1 : http://xerces.apache.org/xerces-c/download.cgi
tomcat: http://tomcat.apache.org/download-70.cgi
