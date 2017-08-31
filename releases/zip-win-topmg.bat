rmdir /S topmg-win-%1
mkdir topmg-win-%1
copy ..\bin\topmg.exe topmg-win-%1 
copy ..\LICENSE topmg-win-%1
xcopy /S ..\toppic_resources topmg-win-%1\toppic_resources

copy C:\mingw64\bin\libgcc_s_seh-1.dll topmg-win-%1
copy "C:\mingw64\bin\libstdc++-6.dll" topmg-win-%1
copy C:\mingw64\bin\libwinpthread-1.dll topmg-win-%1

mkdir topmg-win-%1\example_files
copy .\example_files\mass_graph_mods.txt topmg-win-%1\example_files
