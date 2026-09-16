# Building TopPIC Suite on Red Hat Enterprise Linux

This document explains how to build the command-line and GUI tools of the
TopPIC Suite on Red Hat Enterprise Linux 10 (RHEL 10). The steps also apply
to RHEL-compatible distributions such as Rocky Linux and AlmaLinux 10, except
that the CodeReady Builder repository is enabled with `dnf config-manager`
instead of `subscription-manager` (see step 1).

## 1. Enable the extra repositories

`pugixml-devel` is not in RHEL 10 itself and comes from EPEL (Extra Packages
for Enterprise Linux), which in turn needs the CodeReady Builder repository,
so enable both first.

```sh
# CodeReady Builder (RHEL)
sudo subscription-manager repos --enable codeready-builder-for-rhel-10-$(arch)-rpms
# CodeReady Builder (Rocky Linux / AlmaLinux) -- use this line instead
# sudo dnf config-manager --set-enabled crb

# EPEL
sudo dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-10.noarch.rpm
```

## 2. Install the required packages

```sh
# compilers and build tools (Git LFS is needed for the model files under res/)
sudo dnf install cmake gcc-c++ make clang git git-lfs

# libraries (zlib-ng-compat-devel is RHEL 10's replacement for zlib-devel)
sudo dnf install zlib-ng-compat-devel sqlite-devel pugixml-devel boost-devel

# Qt6 for the GUI tools
sudo dnf install qt6-qtbase-devel
```

Notes:

- The build defaults to **clang/clang++** when they are found; without clang
  it falls back to the system compiler (`g++` works too). To force a compiler,
  pass `-DCMAKE_CXX_COMPILER=...` at configure time.
- `boost-devel` on RHEL 10 provides Boost 1.83, which satisfies the >= 1.74
  requirement.
- The other third-party code (htslib, ProteoWizard, ONNX Runtime) is vendored
  under `ext/` and built or linked automatically; no packages are needed.

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
make -j$(nproc)
```

The build type defaults to `Release`. The executables (`topfd`, `topdia`,
`topindex`, `toppic`, `topmg`, `topdiff` and their `*_gui` counterparts) are
placed in the repository's `bin/` directory. To build just one tool, e.g.
TopFD:

```sh
make -j$(nproc) topfd
```

## 5. Run the tools

The tools can be run directly from `bin/` without installing: they find the
runtime resources in the repository's `res/` directory (they look in `res`
next to the executable first, then in `../res`, then in the installed shared
directory).

## 6. (Optional) Install and uninstall

```sh
sudo make install
```

This installs the binaries to `/usr/local/bin`, the ONNX Runtime shared
library they use to `/usr/local/lib/toppic`, and the runtime resources (model
files, isotope tables, ...) to `/usr/local/share/toppic`. Pass
`-DCMAKE_INSTALL_PREFIX=<dir>` to `cmake` at configure time for a different
prefix.

To remove everything that `make install` placed under the prefix, run the
following from the **same build directory** that was used for `make install`:

```sh
sudo make uninstall
```

The command reads the `install_manifest.txt` written into the build directory
by `make install`, deletes every file listed there, and removes any directories
(e.g. `/usr/local/share/toppic`) that are left empty. If the build directory
has been deleted, there is no manifest and `make uninstall` reports an error;
in that case remove the installed files by hand.
