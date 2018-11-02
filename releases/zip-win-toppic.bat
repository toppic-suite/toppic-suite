rmdir /S toppic-windows-%1
mkdir toppic-windows-%1

copy ..\bin\toppic_gui.exe toppic-windows-%1
copy ..\bin\toppic.exe toppic-windows-%1

copy ..\bin\topfd_gui.exe toppic-windows-%1
copy ..\bin\topfd.exe toppic-windows-%1

copy ..\bin\topmg_gui.exe toppic-windows-%1
copy ..\bin\topmg.exe toppic-windows-%1

copy ..\LICENSE toppic-windows-%1
xcopy /S ..\toppic_resources toppic-windows-%1\toppic_resources

copy C:\msys64\mingw64\bin\libwinpthread-1.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll toppic-windows-%1
copy "C:\msys64\mingw64\bin\libstdc++-6.dll" toppic-windows-%1
copy C:\msys64\mingw64\bin\zlib1.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\libboost_program_options-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_chrono-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_iostreams-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\libxerces-c.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libxalan-c.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libicuuc58.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libxalanMsg.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libicudt58.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\Qt5Core.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\Qt5Gui.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\Qt5Widgets.dll toppic-windows-%1
mkdir toppic-windows-%1\platforms
copy C:\msys64\mingw64\share\qt5\plugins\platforms\qwindows.dll  toppic-windows-%1\platforms

copy C:\msys64\mingw64\bin\libicuin58.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libpng16-16.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libharfbuzz-0.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libgraphite2.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libfreetype-6.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libbz2-1.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libglib-2.0-0.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libintl-8.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libiconv-2.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libpcre-1.dll toppic-windows-%1

mkdir toppic-windows-%1\example_files
copy ..\testcases\data\mzxml_test.msalign toppic-windows-%1\example_files
copy ..\testcases\data\mzxml_test.fasta toppic-windows-%1\example_files
copy ..\testcases\data\mzxml_test.mzXML toppic-windows-%1\example_files
