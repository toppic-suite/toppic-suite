rmdir /S toppic-win-cmd-%1
mkdir toppic-win-cmd-%1
copy ..\bin\toppic.exe toppic-win-cmd-%1 
copy ..\LICENSE toppic-win-cmd-%1
xcopy /S ..\toppic_resources toppic-win-cmd-%1\toppic_resources

copy C:\msys64\mingw64\bin\libwinpthread-1.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libstdc++-6.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\zlib1.dll toppic-win-cmd-%1

copy C:\msys64\mingw64\bin\libboost_program_options-mt.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll toppic-win-cmd-%1

copy C:\msys64\mingw64\bin\libxerces-c.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libxalan-c1_11_0.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libicuuc58.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libxalanMsg1_11_0.dll toppic-win-cmd-%1
copy C:\msys64\mingw64\bin\libicudt58.dll toppic-win-cmd-%1

mkdir toppic-win-cmd-%1\example_files
copy ..\testcases\data\mzxml_test.msalign toppic-win-cmd-%1\example_files
copy ..\testcases\data\mzxml_test.fasta toppic-win-cmd-%1\example_files
