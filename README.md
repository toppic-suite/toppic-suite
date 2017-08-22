# TopPIC: TOP-Down Mass Spectrometry Based Proteoform Identification and Characterization

For manual and refenrence, please visit http://proteomics.informatics.iupui.edu/software/toppic/

## System requirements
* GCC version higher than 4.8.2 for C++11 support
* CMake (>= 2.8)
* For TopPIC-server, Java, ant and Apache Tomcat server

### Linux (Ubuntu):

```sh
sudo apt-get install build-essential cmake unzip zlib1g-dev

cd third_party
unzip linux_include.zip
cd ..

mkdir build
cd build
cmake ..
make -j$(nproc)

cd bin
link -s ../conf ./
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

All the static libraries included were built using GCC 4.8. You might meet the incompatibility during the linking using GCC 5 or other versions of GCC. You can set the version of compiler in CMakeLists.txt as below:

```
SET(CMAKE_C_COMPILER gcc-4.8)
SET(CMAKE_CXX_COMPILER g++-4.8)
```
