rmdir /S topfd-win-cmd-%1
mkdir topfd-win-cmd-%1
copy ..\bin\topfd.exe topfd-win-cmd-%1 
copy ..\LICENSE topfd-win-cmd-%1
xcopy /S ..\toppic_resources topfd-win-cmd-%1\toppic_resources

copy C:\msys64\mingw64\bin\libwinpthread-1.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libstdc++-6.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\zlib1.dll topfd-win-cmd-%1

copy C:\msys64\mingw64\bin\libboost_chrono-mt.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_iostreams-mt.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_program_options-mt.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll topfd-win-cmd-%1

copy C:\msys64\mingw64\bin\libxerces-c.dll topfd-win-cmd-%1

copy C:\msys64\mingw64\bin\libpwiz.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libicuuc58.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libicudt58.dll topfd-win-cmd-%1
copy C:\msys64\mingw64\bin\libbz2-1.dll topfd-win-cmd-%1

mkdir topfd-win-cmd-%1\example_files
copy ..\testcases\data\mzxml_test.mzXML topfd-win-cmd-%1\example_files
