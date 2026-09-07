# Building TopPIC Suite on macOS

This document explains how to build the command-line and GUI tools of the
TopPIC Suite on macOS. Xcode Command Line Tools provide the compiler (clang),
and [Homebrew](https://brew.sh/) provides everything else. Both Apple Silicon
(Homebrew in `/opt/homebrew`) and Intel Macs (Homebrew in `/usr/local`) are
supported.

## 1. Install the compiler and Homebrew

```sh
# Xcode Command Line Tools (clang, make, git)
xcode-select --install
```

Then install Homebrew by following the instructions on
[brew.sh](https://brew.sh/) if it is not installed yet.

## 2. Install the required packages

```sh
# build tools and libraries (Git LFS is needed for the model files under res/)
brew install cmake git git-lfs boost pugixml sqlite zlib onnxruntime

# Qt6 for the GUI tools
brew install qt
```

Notes:

- No ONNX Runtime library is vendored for macOS (the copy under `ext/onnx` is
  a Linux x86-64 binary), so the Homebrew `onnxruntime` package is required.
  CMake locates it with `find_library` in the Homebrew prefix; pass
  `-DONNXRUNTIME_LIBRARY=/path/to/libonnxruntime.dylib` if it is installed
  somewhere else.
- Homebrew's `qt` formula is Qt6 and is linked into the Homebrew prefix, so
  CMake finds it without further hints. If CMake still cannot find Qt6, pass
  `-DCMAKE_PREFIX_PATH="$(brew --prefix qt)"` to `cmake`.
- The other third-party code (htslib, ProteoWizard) is vendored under `ext/`
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
cmake ..
make -j$(sysctl -n hw.ncpu)
```

The build type defaults to `Release` and the compiler to clang. The
executables (`topfd`, `topdia`, `topindex`, `toppic`, `topmg`, `topdiff` and
their `*_gui` counterparts) are placed in the repository's `bin/` directory.
To build just one tool, e.g. TopFD:

```sh
make -j$(sysctl -n hw.ncpu) topfd
```

## 5. Run the tools

The tools can be run directly from `bin/` without installing: they find the
runtime resources in the repository's `res/` directory (they look in `res`
next to the executable first, then in `../res`).

## 6. (Optional) Install and uninstall

```sh
sudo make install
```

This installs the binaries to `/usr/local/bin`. On macOS the runtime resources
are installed next to the executables, in `/usr/local/bin/res`, rather than in
`share/toppic`, because that is where the tools look for them on this
platform. The ONNX Runtime library is not installed; the executables use the
Homebrew one. Pass `-DCMAKE_INSTALL_PREFIX=<dir>` to `cmake` at configure
time for a different prefix.

To remove everything that `make install` placed under the prefix, run the
following from the **same build directory** that was used for `make install`:

```sh
sudo make uninstall
```

The command reads the `install_manifest.txt` written into the build directory
by `make install`, deletes every file listed there, and removes any directories
that are left empty. If the build directory has been deleted, there is no
manifest and `make uninstall` reports an error; in that case remove the
installed files by hand.
