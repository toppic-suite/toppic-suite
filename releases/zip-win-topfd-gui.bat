rmdir /S topfd-win-gui-%1
mkdir topfd-win-gui-%1
copy ..\bin\topfd_gui.exe topfd-win-gui-%1\topfd.exe 
copy ..\LICENSE topfd-win-gui-%1
xcopy /S ..\toppic_resources topfd-win-gui-%1\toppic_resources

copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libwinpthread-1.dll topfd-win-gui-%1
copy "C:\msys64\mingw64\bin\libstdc++-6.dll" topfd-win-gui-%1
copy C:\msys64\mingw64\bin\zlib1.dll topfd-win-gui-%1

copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libboost_chrono-mt.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libboost_iostreams-mt.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll topfd-win-gui-%1

copy C:\msys64\mingw64\bin\libxerces-c.dll topfd-win-gui-%1

copy C:\msys64\mingw64\bin\libpwiz.dll topfd-win-gui-%1

copy C:\msys64\mingw64\bin\Qt5Core.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\Qt5Gui.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\Qt5Widgets.dll topfd-win-gui-%1
mkdir topfd-win-gui-%1\platforms
copy C:\msys64\mingw64\share\qt5\plugins\platforms\qwindows.dll  topfd-win-gui-%1\platforms

copy C:\msys64\mingw64\bin\libicuin58.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libicuuc58.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libharfbuzz-0.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libpng16-16.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libicudt58.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libfreetype-6.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libglib-2.0-0.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libgraphite2.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libbz2-1.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libintl-8.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libpcre-1.dll topfd-win-gui-%1
copy C:\msys64\mingw64\bin\libiconv-2.dll topfd-win-gui-%1

mkdir topfd-win-gui-%1\example_files
copy ..\testcases\data\mzxml_test.mzXML topfd-win-gui-%1\example_files
