@echo off
setlocal EnableExtensions

rem scripts\win_release.bat - package the Windows build of the TopPIC Suite.
rem
rem Collects the twelve executables from bin\, every MSYS2 UCRT64 DLL they
rem import (walked transitively with objdump), the Qt plugins the GUI tools
rem need, the runtime resources in res\ and the LICENSE into a staging folder
rem toppic-win-<version>\ and zips it as toppic-win-<version>.zip, both in
rem the repository root (the parent of this scripts folder). The staging
rem folder is removed afterwards.
rem
rem Usage:   scripts\win_release.bat [version]   (default: the VERSION below)
rem Needs:   a finished build (docs\build\windows_build.md) and MSYS2 with
rem          the UCRT64 packages from that document, installed at
rem          %MSYS2_ROOT% (default C:\msys64).

set "VERSION=1.9.0.0"
if not "%~1"=="" set "VERSION=%~1"

if "%MSYS2_ROOT%"=="" set "MSYS2_ROOT=C:\msys64"
set "UCRT_BIN=%MSYS2_ROOT%\ucrt64\bin"
set "QT_PLUGINS=%MSYS2_ROOT%\ucrt64\share\qt6\plugins"

rem Repository root = the parent of the directory holding this script.
for %%I in ("%~dp0..") do set "ROOT=%%~fI"
set "BIN=%ROOT%\bin"
set "NAME=toppic-win-%VERSION%"
set "STAGE=%ROOT%\%NAME%"
set "ZIP=%ROOT%\%NAME%.zip"
set "TOOLS=topfd topdia topindex toppic topmg topdiff"

rem objdump (binutils, installed with the UCRT64 gcc) walks the DLL imports.
set "PATH=%UCRT_BIN%;%PATH%"

rem ---- checks -------------------------------------------------------------
if not exist "%UCRT_BIN%\objdump.exe" (
  echo ERROR: %UCRT_BIN%\objdump.exe not found.
  echo        Install MSYS2 UCRT64 with the packages listed in docs\build\windows_build.md,
  echo        or set MSYS2_ROOT to the MSYS2 installation directory.
  exit /b 1
)
if not exist "%QT_PLUGINS%\platforms\qwindows.dll" (
  echo ERROR: %QT_PLUGINS%\platforms\qwindows.dll not found ^(mingw-w64-ucrt-x86_64-qt6-base^).
  exit /b 1
)
for %%T in (%TOOLS%) do (
  if not exist "%BIN%\%%T.exe" (
    echo ERROR: %BIN%\%%T.exe not found. Build the suite first.
    exit /b 1
  )
  if not exist "%BIN%\%%T_gui.exe" (
    echo ERROR: %BIN%\%%T_gui.exe not found. Build the suite first.
    exit /b 1
  )
)
if not exist "%ROOT%\res\base_data" (
  echo ERROR: %ROOT%\res\base_data not found.
  exit /b 1
)

echo Packaging TopPIC Suite %VERSION% into %ZIP%

rem ---- staging folder -----------------------------------------------------
if exist "%STAGE%" rmdir /s /q "%STAGE%"
if exist "%ZIP%" del /q "%ZIP%"
mkdir "%STAGE%" || exit /b 1

echo Copying executables and resources ...
for %%T in (%TOOLS%) do (
  copy /y "%BIN%\%%T.exe" "%STAGE%\" >nul || exit /b 1
  copy /y "%BIN%\%%T_gui.exe" "%STAGE%\" >nul || exit /b 1
)
xcopy /e /i /q /y "%ROOT%\res" "%STAGE%\res" >nul || exit /b 1
copy /y "%ROOT%\LICENSE" "%STAGE%\" >nul

rem ---- Qt plugins (loaded at run time by the *_gui tools, not imported) -----
rem Qt looks for these subdirectories next to the executable.
echo Copying Qt plugins ...
mkdir "%STAGE%\platforms"
copy /y "%QT_PLUGINS%\platforms\qwindows.dll" "%STAGE%\platforms\" >nul || exit /b 1
if exist "%QT_PLUGINS%\styles\*.dll" (
  xcopy /i /q /y "%QT_PLUGINS%\styles\*.dll" "%STAGE%\styles\" >nul
)
if exist "%QT_PLUGINS%\imageformats\*.dll" (
  xcopy /i /q /y "%QT_PLUGINS%\imageformats\*.dll" "%STAGE%\imageformats\" >nul
)

rem ---- DLLs ---------------------------------------------------------------
rem Walk the import tables of everything in the staging folder and copy each
rem imported DLL that exists in the UCRT64 bin directory; repeat until a pass
rem adds nothing, so the dependencies of the copied DLLs are picked up too.
rem Windows system DLLs are not in UCRT64\bin and are skipped.
rem Note: onnxruntime.dll must be bundled even though Windows ships one in
rem System32 (an older build); the loader prefers the copy next to the exe.
echo Collecting DLLs from %UCRT_BIN% ...
:walk
set "ADDED="
for /r "%STAGE%" %%F in (*.exe *.dll) do (
  for /f "tokens=3" %%D in ('objdump -p "%%F" ^| findstr /c:"DLL Name:"') do (
    if not exist "%STAGE%\%%D" if exist "%UCRT_BIN%\%%D" (
      copy /y "%UCRT_BIN%\%%D" "%STAGE%\" >nul
      set "ADDED=1"
    )
  )
)
if defined ADDED goto walk

rem onnxruntime.dll loads this one with LoadLibrary, so objdump cannot see it.
if exist "%UCRT_BIN%\onnxruntime_providers_shared.dll" (
  copy /y "%UCRT_BIN%\onnxruntime_providers_shared.dll" "%STAGE%\" >nul
)

rem Warn about imports that are neither bundled nor Windows system DLLs.
for /r "%STAGE%" %%F in (*.exe *.dll) do (
  for /f "tokens=3" %%D in ('objdump -p "%%F" ^| findstr /c:"DLL Name:"') do (
    call :check_import "%%D"
  )
)

set "DLL_COUNT=0"
for %%F in ("%STAGE%\*.dll") do set /a DLL_COUNT+=1
echo Bundled %DLL_COUNT% DLLs.

rem ---- zip ----------------------------------------------------------------
echo Creating %NAME%.zip ...
pushd "%ROOT%"
if exist "%SystemRoot%\System32\tar.exe" (
  "%SystemRoot%\System32\tar.exe" -a -c -f "%NAME%.zip" "%NAME%"
) else (
  powershell -NoProfile -Command "Compress-Archive -Path '%NAME%' -DestinationPath '%NAME%.zip' -Force"
)
set "RC=%ERRORLEVEL%"
popd
if not "%RC%"=="0" (
  echo ERROR: creating the zip file failed.
  exit /b 1
)
rmdir /s /q "%STAGE%"

for %%Z in ("%ZIP%") do echo Done: %%~fZ ^(%%~zZ bytes^)
exit /b 0

rem ---- subroutines --------------------------------------------------------
:check_import
set "IMP=%~1"
if /i "%IMP:~0,11%"=="api-ms-win-" exit /b 0
if exist "%STAGE%\%IMP%" exit /b 0
if exist "%SystemRoot%\System32\%IMP%" exit /b 0
echo WARNING: %IMP% is imported but neither bundled nor a Windows system DLL.
exit /b 0
