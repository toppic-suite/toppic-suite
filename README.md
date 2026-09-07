# toppic_claude

## Git LFS is required

Some of the runtime resources under `res/` are large binary/data blobs
and are stored with [Git LFS](https://git-lfs.com/) rather than in the normal
Git history:

- `res/envcnn_models/*.onnx`, `res/ecscore_models/*.onnx` — the
  EnvCNN / ECScore neural-network models.
- `res/base_data/theo_patt.txt` — the theoretical isotope-pattern table.

You **must** have Git LFS installed to get the real files. Without it, a plain
`git clone` leaves small text *pointer* files in their place, and the build /
the tools will fail to load the models and tables.

### Install Git LFS (one time per machine)

```sh
# Debian/Ubuntu
sudo apt-get install git-lfs
# macOS (Homebrew)
brew install git-lfs
# Then register the Git hooks/filters for your user
git lfs install
```

### Clone

With Git LFS installed, a normal clone fetches the LFS files automatically:

```sh
git clone https://github.com/liuxiaowen/toppic_claude.git
```

### Already cloned without Git LFS?

If you cloned before installing Git LFS (so the `.onnx` / `theo_patt.txt`
files are pointer text), install it as above and then pull the real blobs:

```sh
git lfs install
git lfs pull
```

## Building on Ubuntu Linux

### 1. Install the build dependencies

```sh
sudo apt-get update
sudo apt-get install build-essential cmake clang git git-lfs \
    zlib1g-dev libsqlite3-dev libpugixml-dev \
    libboost-filesystem-dev libboost-iostreams-dev libboost-thread-dev \
    libboost-chrono-dev libboost-system-dev libboost-serialization-dev \
    libboost-program-options-dev \
    qtbase5-dev
```

Notes:

- The build defaults to **clang/clang++** when they are found; if clang is not
  installed, CMake falls back to the system default compiler (g++ works too).
  To force a compiler, pass `-DCMAKE_CXX_COMPILER=...` at configure time.
- **Boost ≥ 1.74** is required. The Ubuntu packages above are sufficient
  (`libboost-all-dev` also works if you prefer one package).
- **Qt5** (`qtbase5-dev`) is needed for the GUI tools (`topfd_gui`, etc.).
- Other third-party code (htslib, ProteoWizard, ONNX Runtime) is vendored
  under `ext/` and built/linked automatically — no packages needed.

### 2. Clone (with Git LFS — see above)

```sh
git clone https://github.com/liuxiaowen/toppic_claude.git
cd toppic_claude
```

### 3. Configure and build

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

### 4. (Optional) Install

```sh
sudo make install
```

This installs the binaries to `/usr/local/bin`, the shared library directory
to `/usr/local/lib/toppic`, and the runtime resources (model files, isotope
tables, ...) to `/usr/local/share/toppic`. Use
`cmake -DCMAKE_INSTALL_PREFIX=<dir> ..` at configure time for a different
prefix.

To run the tools from `bin/` **without** installing, they need to find the
runtime resources in a `res` directory next to the executable; create a
symlink to the repository's `res/` once:

```sh
ln -s ../res bin/res
```
