# TopPIC: TOP-Down Mass Spectrometry Based Proteoform Identification and Characterization

For manual and reference, please visit http://proteomics.informatics.iupui.edu/software/toppic/

## System requirements

* GCC version higher than 4.8.2 for C++11 support
* CMake (>= 2.8)

### Linux (Ubuntu):

```sh
# install compiling tools
sudo apt-get install build-essential cmake

# install the pwiz lib to read proteomics data
sudo apt-get install libpwiz-dev

# install other dependencies
sudo apt-get install zlib1g-dev libboost-filesystem-dev \
                     libboost-program-options-dev \
                     libboost-system-dev \
                     libboost-thread-dev \
                     libxalan-c-dev

mkdir build
cd build
cmake ..
make -j$(nproc)

cd bin
ln -s ../toppic_resources .
```

### Windows:

* Download [CMake](https://cmake.org/) and add it to `PATH` environment;

* Download [mingw64](http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-posix/seh/x86_64-4.8.2-release-posix-seh-rt_v3-rev4.7z/download) and add it to `PATH` environment.

Then

```sh
cd third_party
unzip windows_include.zip

cd ..
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

On some Linux distributions, you might meet the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
 ```

