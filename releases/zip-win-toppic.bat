rmdir /S toppic-windows-%1
mkdir toppic-windows-%1

copy ..\bin\toppic_gui.exe toppic-windows-%1
copy ..\bin\toppic.exe toppic-windows-%1

copy ..\bin\topfd_gui.exe toppic-windows-%1
copy ..\bin\topfd.exe toppic-windows-%1

copy ..\bin\topmg_gui.exe toppic-windows-%1
copy ..\bin\topmg.exe toppic-windows-%1

copy ..\bin\topdiff_gui.exe toppic-windows-%1
copy ..\bin\topdiff.exe toppic-windows-%1

copy ..\bin\topindex_gui.exe toppic-windows-%1
copy ..\bin\topindex.exe toppic-windows-%1

copy ..\bin\topconvert.exe toppic-windows-%1

copy ..\LICENSE toppic-windows-%1
xcopy /S ..\resources toppic-windows-%1\resources

copy C:\msys64\mingw64\bin\libwinpthread-1.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libgcc_s_seh-1.dll toppic-windows-%1
copy "C:\msys64\mingw64\bin\libstdc++-6.dll" toppic-windows-%1
copy C:\msys64\mingw64\bin\zlib1.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\libboost_program_options-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_system-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_filesystem-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_chrono-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_iostreams-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_serialization-mt.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libboost_thread-mt.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\libxerces-c-3-2.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libicuuc69.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libicudt69.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\libzstd.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libcurl-4.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\liblzma-5.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libdouble-conversion.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libbrotlidec.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libbrotlicommon.dll toppic-windows-%1

copy C:\msys64\mingw64\bin\Qt5Core.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\Qt5Gui.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\Qt5Widgets.dll toppic-windows-%1
mkdir toppic-windows-%1\platforms
copy C:\msys64\mingw64\share\qt5\plugins\platforms\qwindows.dll  toppic-windows-%1\platforms

copy C:\msys64\mingw64\bin\libicuin69.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libpng16-16.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libharfbuzz-0.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libgraphite2.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libfreetype-6.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libbz2-1.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libglib-2.0-0.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libintl-8.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libiconv-2.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libpcre2-16-0.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libpcre-1.dll toppic-windows-%1
copy C:\msys64\mingw64\bin\libmd4c.dll toppic-windows-%1

