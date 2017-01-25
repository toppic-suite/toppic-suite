rmdir /S topfd-win-%1
mkdir topfd-win-%1
copy ..\bin\topfd.exe topfd-win-%1 
copy ..\LICENSE topfd-win-%1
xcopy /S ..\toppic_resources topfd-win-%1\toppic_resources
copy C:\mingw64\bin\libgcc_s_seh-1.dll topfd-win-%1
copy "C:\mingw64\bin\libstdc++-6.dll" topfd-win-%1
copy C:\mingw64\bin\libwinpthread-1.dll topfd-win-%1

mkdir topfd-win-%1\example_files
copy ..\testcases\data\mzxml_test.mzXML topfd-win-%1\example_files
