mkdir toppic-win-%1
copy ..\bin\toppic.exe toppic-win-%1 
copy ..\LICENSE toppic-win-%1
xcopy /S ..\toppic_resources toppic-win-%1\toppic_resources
copy .\inno_setup_script\toppic.ico toppic-win-%1
copy C:\mingw64\bin\libgcc_s_seh-1.dll toppic-win-%1
copy "C:\mingw64\bin\libstdc++-6.dll" toppic-win-%1
copy C:\mingw64\bin\libwinpthread-1.dll toppic-win-%1
