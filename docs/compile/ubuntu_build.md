# Building TopPIC Suite on Ubuntu Linux

This document explains how to build the command-line and GUI tools of the
TopPIC Suite on Ubuntu Linux (22.04 or later). The same steps work on other
Debian-based distributions with the matching package names.

## System requirements

* A C++17 compiler: Clang >= 7 (the default when installed) or GCC >= 8
* CMake >= 3.16
* Boost >= 1.74 (filesystem, iostreams, thread, chrono, serialization,
  program_options; Boost.System is header-only and not needed as a library)
* pugixml, SQLite3 and zlib development libraries
* Qt6 (Core, Gui, Widgets) for the GUI tools
* Git LFS, for the model files under `res/` and the vendored ONNX Runtime
  library under `ext/onnx/`; the vendored library is a Linux x86-64 binary,
  so the build needs an x86-64 machine

The packages in step 1 provide all of these on Ubuntu 22.04 and later.

## 1. Install the required packages

```sh
sudo apt-get update
# build tools and libraries (Git LFS is needed for the model files under res/
# and the ONNX Runtime library under ext/onnx/)
sudo apt-get install build-essential cmake clang git git-lfs \
    zlib1g-dev libsqlite3-dev libpugixml-dev \
    libboost-filesystem-dev libboost-iostreams-dev libboost-thread-dev \
    libboost-chrono-dev libboost-serialization-dev \
    libboost-program-options-dev \
    qt6-base-dev
```

The same packages can be installed with the script
[`scripts/install_ubuntu_packages.sh`](../../scripts/install_ubuntu_packages.sh)
(run `scripts/install_ubuntu_packages.sh` from the repository root).

Notes:

- The build defaults to **clang/clang++** when they are found; if clang is not
  installed, CMake falls back to the system default compiler (g++ works too).
  To force a compiler, pass `-DCMAKE_CXX_COMPILER=...` at configure time.
- **Boost ≥ 1.74** is required. The Ubuntu packages above are sufficient
  (`libboost-all-dev` also works if you prefer one package).
- **Qt6** (`qt6-base-dev`) is needed for the GUI tools (`topfd_gui`, etc.).
- Other third-party code (htslib, ProteoWizard, ONNX Runtime) is vendored
  under `ext/` and built/linked automatically — no packages needed.

## 2. Clone the source code

```sh
# register the Git LFS hooks/filters for your user (one time per machine)
git lfs install
git clone https://github.com/toppic-suite/toppic-suite.git
cd toppic-suite
```

If the clone was made before installing Git LFS, run `git lfs pull` inside the
repository to replace the LFS pointer files with the real model files and the
ONNX Runtime library (`ext/onnx/libonnxruntime.so.1.23.2`); without it the
build fails to link.

## 3. Configure and build

```sh
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

The build type defaults to `Release`. The executables (`topfd`, `topdia`,
`topindex`, `toppic`, `topmg`, `topdiff` and their `*_gui` counterparts) are
placed in the repository's `bin/` directory. To build just one tool, e.g.
TopFD:

```sh
make -j$(nproc) topfd
```

## 4. Run the tools

The tools can be run directly from `bin/` without installing: they find the
runtime resources in the repository's `res/` directory (they look in `res`
next to the executable first, then in `../res`, then in the installed shared
directory).

## 5. (Optional) Install

```sh
sudo make install
```

This installs the binaries to `/usr/local/bin`, the ONNX Runtime shared
library they use to `/usr/local/lib/toppic`, and the runtime resources (model
files, isotope tables, ...) to `/usr/local/share/toppic`. Use
`cmake -DCMAKE_INSTALL_PREFIX=<dir> ..` at configure time for a different
prefix.

## 6. (Optional) Uninstall

To remove everything that `make install` placed under the prefix (the
binaries, the ONNX Runtime library and the runtime resources), run the
following from the **same build directory** that was used for `make install`:

```sh
sudo make uninstall
```

The command reads the `install_manifest.txt` written into the build directory
by `make install`, deletes every file listed there, and removes any directories
(e.g. `/usr/local/share/toppic`) that are left empty. If the build directory
has been deleted, there is no manifest and `make uninstall` will report an
error; in that case remove the installed files by hand.
