## System requirements
* GCC version higher than 4.8.2 for C++11 support
* cmake
* For TopPIC-server, Java, ant and Apache Tomcat server

###Linux:

```sh
cd third_party
unzip linux_include.zip

mkdir build
cd build
cmake ..
make -j4 

cd bin
link -s ../conf ./
```

###Windows:

Download [mingw64](http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-posix/seh/x86_64-4.8.2-release-posix-seh-rt_v3-rev4.7z/download) and set it to `PATH` environment.

Then

```sh
cd third_party
unzip windows_include.zip

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
