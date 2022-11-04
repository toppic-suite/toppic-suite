# TopPIC: TOP-Down Mass Spectrometry Based Proteoform Identification and Characterization

For manual and reference, please visit https://www.toppic.org/software/toppic/

## System requirements

* GCC version higher than 5.5.0 for C++14 support
* CMake (>= 3.1)

### Linux (Ubuntu):

```sh
# install compiling tools
sudo apt-get install build-essential cmake

# install other dependencies
sudo apt-get install zlib1g-dev 
sudo apt-get install libxerces-c-dev  
sudo apt-get install libboost-filesystem-dev 
sudo apt-get install libboost-program-options-dev 
sudo apt-get install libboost-system-dev 
sudo apt-get install libboost-thread-dev 
sudo apt-get install libboost-iostreams-dev 
sudo apt-get install libboost-chrono-dev 
sudo apt-get install libeigen3-dev 
sudo apt-get install nlohmann-json3-dev


# install the catch unit test framework (https://github.com/philsquared/Catch)
sudo apt-get install catch

# Qt5 for GUI
sudo apt-get install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)

cd ../bin
ln -s ../toppic_resources .
```

On some Linux distributions, you might meet the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
```

### Windows:

[MSYS2](http://www.msys2.org/) is used for Windows building. Please follow the instructions from [here](doc/windows_build.md).
