Linux:

```sh
cd third_party
unzip linux_include.zip

cd build
cmake ..
make 

cd bin
link -s ../conf ./
```

Windows:

```sh
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```
