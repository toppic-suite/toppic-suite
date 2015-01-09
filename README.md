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
make 

cd bin
link -s ../conf ./
```

###Windows:

```sh
cd third_party
unzip windows_include.zip

mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```
