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

### System requirements

* Clang version >= 16.0.0 for C++17 support
* Boost version >= 1.74.0
* CMake version >= 3.5.0

### Building on Linux (Ubuntu 24.04)

```sh
# install compiling tools
sudo apt install build-essential cmake clang

# install dependencies
sudo apt install libboost-chrono-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev libboost-thread-dev libxerces-c-dev libsqlite3-dev zlib1g-dev 

# install Qt5 for GUI
sudo apt install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

### Building on Linux (Redhat 9)

```sh
# install Extra Packages for Enterprise Linux (EPEL)
sudo subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
sudo dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm

# install compiling tools
sudo dnf install cmake gcc-c++ make clang

# install dependencies
sudo dnf install boost-devel 
sudo dnf install xerces-c-devel
sudo dnf install sqlite3-devel 
sudo dnf install zlib-devel

# install Qt5 for GUI
sudo dnf install qt5-qtbase-devel

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

#### Language setting

On some Linux distributions, you might have the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
```

### Building on Windows

[MSYS2](http://www.msys2.org/) is used for building TopPIC Suite on Windows systems. Please follow the instructions from [here](doc/windows_build.md).

