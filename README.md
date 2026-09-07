## TopPIC Suite

TopPIC Suite consists of six software tools for the analysis of top-down mass spectrometry-based proteomics data. 

* **TopFD** (Top-down mass spectral Feature Detection) is a software tool for top-down mass spectral deconvolution. It groups top-down spectral peaks into isotopic envelopes and converts isotopic envelopes to monoisotopic neutral masses. In addition, it extracts proteoform features from LC-MS or CE-MS data.

* **TopDIA** is a tool for demultiplexing top-down data independent acquisition mass spectrometry (TD-DIA-MS) data. It processes TD-DIA-MS data to generate demultiplexed pseudo-MS/MS spectra, which are subsequently searched against a protein database to identify proteoforms.

* **TopIndex** (Top-down protein sequence database Indexing) generates index files for protein sequence databases. The index files are used in TopPIC and TopMG to speed up proteoform identification by database search.

* **TopPIC** (Top-down mass spectrometry-based Proteoform Identification and Characterization) identifies and characterizes proteoforms at the proteome level by searching top-down tandem mass spectra against a protein sequence database. It efficiently identifies proteoforms with post-translational modificatons (PTMs) and unexpected alterations, such as mutations, accurately estimates the statistical significance of identifications, and characterizes reported proteoforms with unknown mass shifts. It uses several techniques, such as indexes, spectral alignment, generating function methods, and the modification identification score (MIScore), to increase the speed, sensitivity, and accuracy.

* **TopMG** (Top-down mass spectrometry-based proteoform identification using Mass Graphs) is a software tool for identifying highly modified proteoforms by searching top-down tandem mass spectra against a protein sequence database. It is capable of identifying proteoforms with multiple variable PTMs and unexpected alterations, such as histone proteoforms and phosphorylated ones. It uses mass graphs, which efficiently represent candidate proteoforms with multiple variable PTMs, to increase the speed and sensitivity in proteoform identification. In addition, approximate spectrum-based filtering methods are employed for protein sequence filtering, and a Markov chain Monte Carlo method (TopMCMC) is used for estimating the statistical significance of identifications.

* **TopDiff** (Top-down mass spectrometry-based identification of Differentially expressed proteoforms) compares the abundances of proteoforms and finds differentially expressed proteoforms by using identifications of top-down mass spectrometry data of several protein samples.

* **TopDIA** is a software tool for top-down data-independent-acquistion mass spectrometry (TD-DIA-MS) data analysis. It generates demultiplexed pseudo MS/MS spectra from TD-DIA-MS data, which are then searched against a protein sequence database using TopPIC or TopMG for proteoform identification.

**For manuals, tutorials, and publications, please visit https://www.toppic.org/software/toppic/.** 

## System requirements

* A C++17 compiler: Clang >= 7 (default) or GCC >= 8
* CMake version >= 3.16
* Boost version >= 1.74 (filesystem, iostreams, thread, chrono, system, serialization, program_options)
* pugixml, SQLite3 and zlib development libraries
* Qt5 (Core, Gui, Widgets) for the GUI tools
* Git LFS (see below)

### Git LFS is required

Some of the runtime resources under `res/` are large binary/data blobs
and are stored with [Git LFS](https://git-lfs.com/) rather than in the normal
Git history:

- `res/envcnn_models/*.onnx`, `res/ecscore_models/*.onnx` — the
  EnvCNN / ECScore neural-network models.
- `res/base_data/theo_patt.txt` — the theoretical isotope-pattern table.

You **must** have Git LFS installed to get the real files. Without it, a plain
`git clone` leaves small text *pointer* files in their place, and the build /
the tools will fail to load the models and tables.

#### Install Git LFS (one time per machine)

```sh
# Debian/Ubuntu
sudo apt-get install git-lfs
# macOS (Homebrew)
brew install git-lfs
# Then register the Git hooks/filters for your user
git lfs install
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
git clone https://github.com/toppic_suite/toppic_suite.git
cd toppic_suite
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

### 4. Install

```sh
sudo make install
```

This installs the binaries to `/usr/local/bin`, the ONNX Runtime shared
library they use to `/usr/local/lib/toppic`, and the runtime resources (model
files, isotope tables, ...) to `/usr/local/share/toppic`. Use
`cmake -DCMAKE_INSTALL_PREFIX=<dir> ..` at configure time for a different
prefix.

The tools can also be run directly from `bin/` **without** installing: they
find the runtime resources in the repository's `res/` directory (the tools
look in `res` next to the executable first, then in `../res`, then in the
installed shared directory).

### 5. Uninstall

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

## Building on Red Hat Enterprise Linux 9

The steps are the same as for Ubuntu above; only the package installation
differs. `pugixml-devel` comes from EPEL (Extra Packages for Enterprise Linux),
which in turn needs the CodeReady Builder repository, so enable both first.

```sh
# enable CodeReady Builder and EPEL
sudo subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
sudo dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm

# compilers and build tools
sudo dnf install cmake gcc-c++ make clang git git-lfs

# libraries
sudo dnf install zlib-devel sqlite-devel pugixml-devel boost-devel

# Qt5 for the GUI tools
sudo dnf install qt5-qtbase-devel
```

`boost-devel` on RHEL 9 provides Boost 1.75, which satisfies the ≥ 1.74
requirement. Then continue with steps 2–5 of the Ubuntu instructions (clone
with Git LFS, configure and build, install, uninstall); they are identical.
On Red Hat the same commands apply, e.g.:

```sh
git lfs install
git clone https://github.com/toppic_suite/toppic_suite.git
cd toppic_suite
mkdir -p build
cd build
cmake ..
make -j$(nproc)
sudo make install
```

## Building on Windows

[MSYS2](http://www.msys2.org/) is used for building TopPIC Suite on Windows systems. Please follow the instructions from [here](doc/windows_build.md).


## Building on macOS

The steps are the same as for Ubuntu above; only the package installation and
two configure details differ. Xcode Command Line Tools provide clang, and
[Homebrew](https://brew.sh/) provides everything else. No ONNX Runtime library
is vendored for macOS, so the Homebrew package is used instead.

```sh
# compiler
xcode-select --install

# build tools and libraries
brew install cmake git git-lfs boost pugixml sqlite zlib onnxruntime

# Qt5 for the GUI tools
brew install qt@5
```

Then clone as in step 2 of the Ubuntu instructions and build:

```sh
mkdir -p build
cd build
cmake -DCMAKE_PREFIX_PATH="$(brew --prefix qt@5)" ..
make -j$(sysctl -n hw.ncpu)
sudo make install
```

Notes:

- `qt@5` is a keg-only Homebrew formula (it is not linked into the Homebrew
  prefix), so CMake cannot find it unless its location is passed via
  `CMAKE_PREFIX_PATH` as above.
- On macOS the runtime resources are installed next to the executables, in
  `<prefix>/bin/res`, rather than in `share/toppic`, because that is where the
  tools look for them on this platform. `make uninstall` works as on Linux.
