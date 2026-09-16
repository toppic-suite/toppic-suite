---
name: compile
description: Build, run, install and uninstall the TopPIC Suite (topfd, topdia, topindex, toppic, topmg, topdiff and their Qt GUIs) from source on Ubuntu, Red Hat, Windows (MSYS2) or macOS. Use when asked to build or rebuild the suite or one tool, to fix a configure/build failure, to install or uninstall it, or to set up a machine for building.
---

# Building the TopPIC Suite

The build is CMake-based (C++17, `CMakeLists.txt` at the repository root).
All library code compiles into one OBJECT library, `toppic_common`, which is
linked into each executable; no shared library is produced. The executables
land in `<repo>/bin/` and run from there without installing, because they
find the runtime resources in `<repo>/res/` (`res` next to the executable,
then `../res`, then the installed share directory).

The authoritative, user-facing instructions are `README.md` (Ubuntu) and
`docs/build/{redhat,windows,macos}_build.md`. Keep them in step with
`CMakeLists.txt` when the build changes.

## Quick path (an existing checkout on Linux or macOS)

```sh
cd <repo>
mkdir -p build && cd build
cmake ..                 # Release by default; clang preferred when found
make -j$(nproc)          # macOS: make -j$(sysctl -n hw.ncpu)
```

- Build one tool: `make -j$(nproc) topfd` (targets: `topfd`, `topdia`,
  `topindex`, `toppic`, `topmg`, `topdiff`, and `topfd_gui`, `topdia_gui`,
  `topindex_gui`, `toppic_gui`, `topmg_gui`, `topdiff_gui`).
- Generator-independent form: `cmake --build build --target toppic -j8`.
- Verify: `bin/toppic -h` (or the tool you built) prints the option list.
- A `build/` directory usually already exists in this checkout; reuse it
  instead of configuring a second one.

## Requirements (all platforms)

- C++17 compiler: Clang >= 7 (default) or GCC >= 8. Force one with
  `-DCMAKE_CXX_COMPILER=...`.
- CMake >= 3.16.
- Boost >= 1.74: filesystem, iostreams, thread, chrono, system,
  serialization, program_options.
- pugixml, SQLite3, zlib development packages.
- Qt6 (Core, Gui, Widgets) for the `*_gui` targets only.
- Git LFS: the EnvCNN/ECScore model files and the isotope table under
  `res/` are LFS objects. Run `git lfs install` once per machine before
  cloning; if the clone predates it, run `git lfs pull` in the repository.
  Symptom of missing LFS: small text pointer files under `res/`, and the
  tools failing to load models at run time.
- Vendored, nothing to install: htslib (`ext/htslib`), ProteoWizard
  (`ext/pwiz`), and on Linux x86-64 the ONNX Runtime shared library
  (`ext/onnx/libonnxruntime.so.1.14.1`).

## Platform notes

### Ubuntu

```sh
sudo apt-get update
sudo apt-get install build-essential cmake clang git git-lfs \
    zlib1g-dev libsqlite3-dev libpugixml-dev \
    libboost-filesystem-dev libboost-iostreams-dev libboost-thread-dev \
    libboost-chrono-dev libboost-system-dev libboost-serialization-dev \
    libboost-program-options-dev \
    qt6-base-dev
```

`libboost-all-dev` works in place of the individual Boost packages. Then
follow the quick path above.

### Red Hat Enterprise Linux 10 (Rocky, AlmaLinux)

`pugixml-devel` comes from EPEL, which needs CodeReady Builder:

```sh
sudo subscription-manager repos --enable codeready-builder-for-rhel-10-$(arch)-rpms
# Rocky/AlmaLinux instead: sudo dnf config-manager --set-enabled crb
sudo dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-10.noarch.rpm
sudo dnf install cmake gcc-c++ make clang git git-lfs
sudo dnf install zlib-ng-compat-devel sqlite-devel pugixml-devel boost-devel
sudo dnf install qt6-qtbase-devel
```

`zlib-ng-compat-devel` is RHEL 10's zlib; `boost-devel` is Boost 1.83.
Then follow the quick path.

### Windows (MSYS2 UCRT64)

Use the **UCRT64** shell (ONNX Runtime is packaged only for UCRT64 and
CLANG64, not MINGW64).

```sh
pacman -S git mingw-w64-ucrt-x86_64-git-lfs
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja
pacman -S mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-pugixml \
          mingw-w64-ucrt-x86_64-sqlite3 mingw-w64-ucrt-x86_64-zlib \
          mingw-w64-ucrt-x86_64-onnxruntime
pacman -S mingw-w64-ucrt-x86_64-qt6-base
mkdir -p build && cd build
cmake -G Ninja ..
ninja
```

- Ninja is required: the "MinGW Makefiles" generator refuses to run inside
  an MSYS2 shell.
- The compiler must be **GCC** in UCRT64. `CMakeLists.txt` deliberately
  keeps CMake's default compiler on Windows hosts; do not pass clang. The
  UCRT64 static Boost is GCC-built, and mixing clang objects with it fails
  at link time with "duplicate section ... has different size" then
  "multiple definition" errors. For clang use the CLANG64 shell with the
  `mingw-w64-clang-x86_64-*` packages.
- The system ONNX Runtime (`libonnxruntime.dll.a`) is located with
  `find_library`; override with `-DONNXRUNTIME_LIBRARY=...`.
- Running outside the UCRT64 shell needs `C:\msys64\ucrt64\bin` on `PATH`
  for the Boost/Qt6/ONNX DLLs.
- Install: `cmake --install .` or `ninja install`; `res/` is installed next
  to the executables. `ninja uninstall` reverses it.

### macOS (Xcode CLT + Homebrew)

```sh
xcode-select --install
brew install cmake git git-lfs boost pugixml sqlite zlib onnxruntime
brew install qt            # Qt6
mkdir -p build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

- Nothing ONNX-related is vendored for macOS; Homebrew `onnxruntime` is
  required and found in the Homebrew prefix (`-DONNXRUNTIME_LIBRARY=` to
  override). Apple Silicon Homebrew lives in `/opt/homebrew`; put its `bin`
  on `PATH` if `brew` is not found.
- If Qt6 is not found: `-DCMAKE_PREFIX_PATH="$(brew --prefix qt)"`.
- Install puts resources in `<prefix>/bin/res`, not `share/toppic`, and
  does not install an ONNX library.

## Configure options

| Option | Effect |
|---|---|
| `-DCMAKE_BUILD_TYPE=Debug` | Debug build (default is `Release`). |
| `-DCMAKE_CXX_COMPILER=g++` | Choose the compiler (default clang when found, except on Windows hosts). |
| `-DCMAKE_INSTALL_PREFIX=<dir>` | Install prefix (default `/usr/local`; MSYS2's prefix on Windows). |
| `-DONNXRUNTIME_LIBRARY=<file>` | ONNX Runtime library on macOS/Windows when `find_library` misses it. |
| `-DTOPPIC_ENABLE_IWYU=ON` | Run include-what-you-use during compilation (slow; needs the tool installed). `iwyu_tool -p build` is the standalone alternative. |

## Install and uninstall (Linux)

```sh
cd build
sudo make install
sudo make uninstall      # same build directory; reads install_manifest.txt
```

`make install` puts the twelve executables in `<prefix>/bin`, the vendored
ONNX Runtime library in `<prefix>/lib/toppic` (the executables carry an
`$ORIGIN/../lib/toppic` rpath, so the tree is relocatable), and `res/` in
`<prefix>/share/toppic`. `make uninstall` deletes every file in the
manifest and then any directory left empty; without the manifest (build
directory deleted) it errors instead of guessing, and the files must be
removed by hand.

## Troubleshooting

- **"ONNX Runtime library not found"** (macOS/Windows): install the system
  package (`brew install onnxruntime`, `pacman -S
  mingw-w64-ucrt-x86_64-onnxruntime`) or pass `-DONNXRUNTIME_LIBRARY`.
- **Link errors about duplicate sections / multiple definitions on
  Windows**: clang objects against GCC Boost. Rebuild in a clean `build/`
  with the UCRT64 GCC.
- **Qt6 not found**: install `qt6-base-dev` / `qt6-qtbase-devel` /
  `mingw-w64-ucrt-x86_64-qt6-base` / `brew install qt`, or point
  `CMAKE_PREFIX_PATH` at it. Only the `*_gui` targets need it.
- **Boost too old**: >= 1.74 is required; check `cmake` output for the
  version it found.
- **Model files are tiny text files / tools cannot load models**: Git LFS
  objects were not fetched; run `git lfs install && git lfs pull`.
- **Stale configuration** after changing compiler or packages: delete
  `build/CMakeCache.txt` (or the whole `build/`) and configure again.
- **Adding an executable target**: append it to the `TOPPIC_EXECUTABLES`
  list at the end of `CMakeLists.txt` so it is installed and gets the
  rpath, then update `README.md` and the `docs/build/` documents that list
  the tools.
