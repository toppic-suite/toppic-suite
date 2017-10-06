.PHONY : clean

CXX=g++

CPPFLAGS = -std=c++11 -Wall -O3 -DWIN32 -DWITHOUT_MZ5 -fPIC -I../libraries/boost_aux -I./ -I..

LDFLAGS= -shared -lboost_system-mt -lboost_filesystem-mt -lboost_thread-mt -lboost_chrono-mt -lboost_iostreams-mt -lz -lWs2_32

SOURCES = \
./utility/minimxml/SAXParser.cpp \
./utility/minimxml/XMLWriter.cpp \
./utility/misc/Filesystem.cpp \
./utility/misc/random_access_compressed_ifstream.cpp \
./utility/misc/SHA1Calculator.cpp \
./utility/misc/SHA1.cpp \
./utility/misc/Base64.cpp \
./utility/misc/IterationListener.cpp \
./utility/misc/IntegerSet.cpp \
./data/common/ParamTypes.cpp \
./data/common/cv.cpp \
./data/common/diff_std.cpp \
./data/common/CVTranslator.cpp \
./data/msdata/SpectrumInfo.cpp \
./data/msdata/Reader.cpp \
./data/msdata/MSData.cpp \
./data/msdata/MSDataFile.cpp \
./data/msdata/BinaryDataEncoder.cpp \
./data/msdata/Serializer_MGF.cpp \
./data/msdata/Serializer_MSn.cpp \
./data/msdata/Serializer_mzML.cpp \
./data/msdata/Serializer_mzXML.cpp \
./data/msdata/SpectrumList_MGF.cpp \
./data/msdata/SpectrumList_MSn.cpp \
./data/msdata/SpectrumList_mzML.cpp \
./data/msdata/SpectrumList_mzXML.cpp \
./data/msdata/DefaultReaderList.cpp \
./data/msdata/MSNumpress.cpp \
./data/msdata/SpectrumWorkerThreads.cpp \
./data/msdata/IO.cpp \
./data/msdata/Index_mzML.cpp \
./data/msdata/ChromatogramList_mzML.cpp \
./data/msdata/LegacyAdapter.cpp \
./data/msdata/References.cpp \
./data/msdata/Version.cpp \
./data/msdata/SpectrumList_BTDX.cpp

OBJECTS=$(SOURCES:.cpp=.o)

TARGET=libpwiz.dll

all:    $(TARGET)

$(TARGET):	$(OBJECTS)
	$(CXX) $(CFLAGS) $(OBJECTS) -o libpwiz.dll $(LDFLAGS)

clean:
	del $(OBJECTS) $(TARGET)
