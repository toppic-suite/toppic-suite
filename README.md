# TopPIC: TOP-Down Mass Spectrometry Based Proteoform Identification and Characterization

For manual and reference, please visit https://www.toppic.org/software/toppic/

## System requirements

* GCC version higher than 5.5.0 for C++14 support
* CMake (>= 3.1)

### Linux (Ubuntu 22.04):

```sh
# install compiling tools
sudo apt install build-essential cmake

# install other dependencies
sudo apt install zlib1g-dev 
sudo apt install libboost-filesystem-dev 
sudo apt install libboost-program-options-dev 
sudo apt install libboost-system-dev 
sudo apt install libboost-thread-dev 
sudo apt install libboost-iostreams-dev 
sudo apt install libboost-chrono-dev 

sudo apt install libxerces-c-dev  
sudo apt install libeigen3-dev 
sudo apt install nlohmann-json3-dev


# install the catch unit test framework (https://github.com/philsquared/Catch)
sudo apt install catch

# Qt5 for GUI
sudo apt install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

### Linux (CentOS Stream 8):

```sh
# install Extra Packages for Enterprise Linux (EPEL)
sudo dnf install 'dnf-command(config-manager)'
sudo dnf config-manager --set-enabled powertools
sudo dnf install epel-release 

# install compiling tools
sudo dnf install gcc gcc-c++ make cmake

# install dependencies
sudo dnf install zlib-devel
sudo dnf install boost-devel 
sudo dnf install xerces-c-devel
sudo dnf install eigen3-devel
sudo dnf install json-devel

# Qt5 for GUI
sudo dnf install qt5-qtbase-devel

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

### Language setting

On some Linux distributions, you might meet the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
```

### Windows:

[MSYS2](http://www.msys2.org/) is used for Windows building. Please follow the instructions from [here](doc/windows_build.md).
