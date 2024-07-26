# TopPIC: TOP-Down Mass Spectrometry Based Proteoform Identification and Characterization

For manual and reference, please visit https://www.toppic.org/software/toppic/

## System requirements

* GCC-C++ version higher than 11.4.0 for C++17 support
* CMake (>= 3.5)

### Linux (Ubuntu 24.04):

```sh
# install compiling tools
sudo apt install build-essential cmake

# install dependencies
sudo apt install libboost-chrono-dev 
sudo apt install libboost-filesystem-dev 
sudo apt install libboost-iostreams-dev 
sudo apt install libboost-program-options-dev 
sudo apt install libboost-thread-dev 
sudo apt install libxerces-c-dev  
sudo apt install zlib1g-dev 

# install Qt5 for GUI
sudo apt install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

### Linux (Redhat 9):

```sh
# install Extra Packages for Enterprise Linux (EPEL)
sudo subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
sudo dnf install https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm

# install compiling tools
sudo dnf install cmake gcc-c++ make

# install dependencies
sudo dnf install boost-devel 
sudo dnf install xerces-c-devel
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

### Language setting

On some Linux distributions, you might have the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
```

### Windows:

[MSYS2](http://www.msys2.org/) is used for Windows building. Please follow the instructions from [here](doc/windows_build.md).
