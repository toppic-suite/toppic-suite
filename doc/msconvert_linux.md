This document provides the instruction to run msconvert to convert RAW format on Linux.

`msconvert` from ProteoWizard is designed to run on Windows, however, it is possible to run it using `wine`. The following instructions were test on Ubuntu 14.04.

## install wine

```sh
sudo add-apt-repository ppa:ubuntu-wine/ppa
sudo apt-get -y update
sudo apt-get -y install wine
sudo apt-get -y install winetricks
```

## set environment variable

```sh
export WINEARCH=win32      ## force to run in 32-bit mode
export WINEPREFIX=~/.wine  ## change as you want
```

## install essential dlls

```sh
## VC runtime enviroment
winetricks vcrun2005 vcrun2008 vcrun2010   
## dotnet40
winetricks dotnet40
## for GUI program
winetricks gdiplus          
```

## download and unzip pwiz

```sh
## load the dll for RAW format
wine regsvr32 MSFileReader.XRawfile2.dll
## run msconvert.exe
wine msconvert.exe
```
