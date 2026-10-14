# Building TopPIC Suite on Windows

This document explains how to build the command-line and GUI tools of the
TopPIC Suite on Windows with [MSYS2](https://www.msys2.org).

## System requirements

* The MSYS2 UCRT64 environment
* GCC (the MSYS2 `mingw-w64-ucrt-x86_64-gcc` package), not clang: the
  MSYS2 libraries are built with GCC, and mixing compilers fails at link time
* CMake >= 3.16 and Ninja
* Boost >= 1.74 (filesystem, iostreams, thread, chrono, serialization,
  program_options, and nowide, which the ProteoWizard code needs on Windows)
* pugixml, SQLite3, zlib and zstd development libraries
* ONNX Runtime from MSYS2 (the vendored library under `ext/onnx/` is for
  Linux only)
* Qt6 (Core, Gui, Widgets) for the GUI tools
* Git LFS, for the model files under `res/`

The packages in step 2 provide all of these.

## 1. Install MSYS2

Follow the instructions on the MSYS2 website to install MSYS2 and to update
the package database and the core system packages (`pacman -Syu`, restart
the shell, `pacman -Su`).

All of the following steps are done in the **MSYS2 UCRT64** shell (the
"MSYS2 UCRT64" entry in the Start menu, or `C:\msys64\ucrt64.exe`). The UCRT64
environment is used because MSYS2 packages ONNX Runtime, which TopFD and
TopDIA need, only for its UCRT64 and CLANG64 environments, not for MINGW64.

## 2. Install the required packages

```sh
# version control (Git LFS is needed for the model files under res/; the
# vendored Linux ONNX Runtime library under ext/onnx/ is LFS-tracked too but
# not used on Windows)
pacman -S git mingw-w64-ucrt-x86_64-git-lfs

# compiler and build tools
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja

# libraries
pacman -S mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-pugixml \
          mingw-w64-ucrt-x86_64-sqlite3 mingw-w64-ucrt-x86_64-zlib \
          mingw-w64-ucrt-x86_64-zstd \
          mingw-w64-ucrt-x86_64-onnxruntime

# Qt6 for the GUI tools
pacman -S mingw-w64-ucrt-x86_64-qt6-base
```

The other third-party code (htslib, ProteoWizard) is vendored under `ext/`
and built automatically.

## 3. Clone the source code

```sh
git lfs install
git clone https://github.com/toppic-suite/toppic-suite.git
cd toppic-suite
```

If the clone was made before installing Git LFS, run `git lfs pull` inside the
repository to replace the LFS pointer files with the real model files.

## 4. Configure and build

```sh
mkdir -p build
cd build
cmake -G Ninja ..
ninja
```

The Ninja generator is used because CMake's "MinGW Makefiles" generator
refuses to run inside an MSYS2 shell (it rejects `sh.exe` on the `PATH`). On
Windows the build uses the environment's default compiler, which is GCC in
the UCRT64 shell (on Linux and macOS it defaults to clang). Do not build with
clang in the UCRT64 shell: the UCRT64 packages, including the static Boost
libraries, are built with GCC, and linking clang-compiled objects against them
fails with "duplicate section ... has different size" warnings followed by
"multiple definition" errors. To use clang, use the MSYS2 CLANG64 shell with
the matching `mingw-w64-clang-x86_64-*` packages instead. The executables
are placed in the repository's `bin/` directory: the ten command-line tools
(`topfd.exe`, `topdia.exe`, `topindex.exe`, `toppic.exe`, `topmg.exe`,
`topdiff.exe`, `toplib_build.exe`, `toplib_search.exe`, `topcalib.exe` and
`topzstd.exe`) and the six GUI tools (`topfd_gui.exe`, `topdia_gui.exe`,
`topindex_gui.exe`, `toppic_gui.exe`, `topmg_gui.exe` and `topdiff_gui.exe`).

## 5. Run the tools

The tools can be run directly from `bin/`: they find the runtime resources in
the repository's `res/` directory (they look in `res` next to the executable
first, then in `..\res`).

The executables depend on the DLLs of the UCRT64 environment (Boost, Qt6,
ONNX Runtime, ...), which live in `C:\msys64\ucrt64\bin`. Inside the UCRT64
shell that directory is already on the `PATH`; to run the tools from a
Windows Terminal or PowerShell instead, add `C:\msys64\ucrt64\bin` to the
`PATH` environment variable.

## 6. (Optional) Install and uninstall

`cmake --install .` (or `ninja install`) run in the build directory installs
the executables and the `res` directory next to them under the install
prefix, which defaults to the UCRT64 prefix (`C:\msys64\ucrt64`). Pass
`-DCMAKE_INSTALL_PREFIX=<dir>` to `cmake` at configure time to install
somewhere else. `ninja uninstall` from the same build directory removes the
installed files again.

## 7. (Optional) Package a release zip

`scripts\win_release.bat` bundles a finished build into a self-contained
`toppic-windows-<version>.zip` in the repository root (default version
`2.1.0`; pass another as the first argument). Run it from a Windows
command prompt or PowerShell in the repository root after building:

```bat
scripts\win_release.bat
scripts\win_release.bat 2.1.1
```

It copies the twelve executables of the six main tools and their GUIs
from `bin\` (not `toplib_build.exe`, `toplib_search.exe`, `topcalib.exe` or
`topzstd.exe`), walks their import tables
with `objdump` and copies every UCRT64 DLL they need (Boost is linked
statically; the runtime DLLs are the GCC runtime, pugixml, SQLite, zlib, zstd,
ONNX Runtime and, for the GUIs, Qt6 and its dependencies), the Qt
`platforms`, `styles` and `imageformats` plugins, the `res` directory and
the `LICENSE` into `toppic-windows-<version>\`, zips that folder with the
built-in `tar.exe` (or `Compress-Archive`) and deletes the staging folder.
The extracted folder runs without MSYS2 on the `PATH`. If MSYS2 is not
installed at `C:\msys64`, set the `MSYS2_ROOT` environment variable first.
