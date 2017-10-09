rmdir /S topmg-win-cmd-%1
mkdir topmg-win-cmd-%1
copy ..\bin\topmg.exe topmg-win-cmd-%1 
copy LICENSE topmg-win-cmd-%1
xcopy /S ..\toppic_resources topmg-win-cmd-%1\toppic_resources

copy C:\msys64\mingw64\bin\libwinpthread-1.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll topmg-win-cmd-%1
copy "C:\msys64\mingw64\bin\libstdc++-6.dll" topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\zlib1.dll topmg-win-cmd-%1

copy C:\msys64\mingw64\bin\libboost_program_options-mt.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll topmg-win-cmd-%1

copy C:\msys64\mingw64\bin\libxerces-c.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libxalan-c1_11_0.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libicuuc58.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libxalanMsg1_11_0.dll topmg-win-cmd-%1
copy C:\msys64\mingw64\bin\libicudt58.dll topmg-win-cmd-%1

mkdir topmg-win-cmd-%1\example_files
copy .\example_files\mass_graph_mods.txt topmg-win-cmd-%1\example_files
