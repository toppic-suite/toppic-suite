# Building TopPIC Suite on Windows

This document explains how to build the command-line and GUI tools of the
TopPIC Suite on Windows with [MSYS2](https://www.msys2.org).

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
# version control (Git LFS is needed for the model files under res/)
pacman -S git mingw-w64-ucrt-x86_64-git-lfs

# compiler and build tools
pacman -S mingw-w64-ucrt-x86_64-clang mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja

# libraries
pacman -S mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-pugixml \
          mingw-w64-ucrt-x86_64-sqlite3 mingw-w64-ucrt-x86_64-zlib \
          mingw-w64-ucrt-x86_64-onnxruntime

# Qt5 for the GUI tools
pacman -S mingw-w64-ucrt-x86_64-qt5-base
```

The other third-party code (htslib, ProteoWizard) is vendored under `ext/`
and built automatically.

## 3. Clone the source code

```sh
git lfs install
git clone https://github.com/toppic_suite/toppic_suite.git
cd toppic_suite
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
refuses to run inside an MSYS2 shell (it rejects `sh.exe` on the `PATH`). The
build defaults to clang when it is found, otherwise to gcc. The executables
(`topfd.exe`, `topdia.exe`, `topindex.exe`, `toppic.exe`, `topmg.exe`,
`topdiff.exe` and their `*_gui.exe` counterparts) are placed in the
repository's `bin/` directory.

## 5. Run the tools

The tools look for the runtime resources in a `res` directory next to the
executable. To run them from `bin/`, copy the repository's `res/` there once:

```sh
cd ..
cp -r res bin/
```

The executables depend on the DLLs of the UCRT64 environment (Boost, Qt5,
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
