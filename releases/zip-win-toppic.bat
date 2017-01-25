rmdir /S toppic-win-%1
mkdir toppic-win-%1
copy ..\bin\toppic.exe toppic-win-%1 
copy ..\LICENSE toppic-win-%1
xcopy /S ..\toppic_resources toppic-win-%1\toppic_resources

copy C:\mingw64\bin\libgcc_s_seh-1.dll toppic-win-%1
copy "C:\mingw64\bin\libstdc++-6.dll" toppic-win-%1
copy C:\mingw64\bin\libwinpthread-1.dll toppic-win-%1

mkdir toppic-win-%1\example_files
copy ..\testcases\data\mzxml_test.msalign toppic-win-%1\example_files
copy ..\testcases\data\mzxml_test.fasta toppic-win-%1\example_files
copy ..\localization_mods.txt toppic-win-%1\example_files
