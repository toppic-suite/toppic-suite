include ../Makefile.in

CPPFLAGS = -Wall -O3 -DPROT_LINUX -DWITHOUT_MZ5


TEST_OBJS =  $(addprefix $(OBJ_DIR)/, test_pwiz.o)


PWIZ_OBJ_ARCH = $(OBJ_DIR)

#
# define the executable file 
MAIN = $(BIN_DIR)/test_pwiz

PWIZ_OBJS = \
            $(PWIZ_OBJ_ARCH)/utility/minimxml/SAXParser.o \
            $(PWIZ_OBJ_ARCH)/utility/minimxml/XMLWriter.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/Filesystem.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/random_access_compressed_ifstream.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/SHA1Calculator.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/SHA1.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/Base64.o \
            $(PWIZ_OBJ_ARCH)/utility/misc/IterationListener.o \
            $(PWIZ_OBJ_ARCH)/data/common/ParamTypes.o \
      			$(PWIZ_OBJ_ARCH)/data/common/cv.o \
			      $(PWIZ_OBJ_ARCH)/data/common/diff_std.o \
            $(PWIZ_OBJ_ARCH)/data/common/CVTranslator.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumInfo.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Reader.o \
			      $(PWIZ_OBJ_ARCH)/data/msdata/MSData.o \
			      $(PWIZ_OBJ_ARCH)/data/msdata/MSDataFile.o \
					  $(PWIZ_OBJ_ARCH)/data/msdata/BinaryDataEncoder.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Serializer_MGF.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Serializer_MSn.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Serializer_mzML.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Serializer_mzXML.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumList_MGF.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumList_MSn.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumList_mzML.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumList_mzXML.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/DefaultReaderList.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/MSNumpress.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/IO.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Index_mzML.o \
						$(PWIZ_OBJ_ARCH)/data/msdata/ChromatogramList_mzML.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/LegacyAdapter.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/References.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/Version.o \
            $(PWIZ_OBJ_ARCH)/data/msdata/SpectrumList_BTDX.o \


#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
all:    $(MAIN)

$(MAIN): $(PWIZ_OBJS) $(TEST_OBJS)
	ar rcs libpwiz.a $(PWIZ_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(PWIZ_OBJS) $(TEST_OBJS) $(LFLAGS) $(LIBS)

#$(OBJ_DIR)/%.o: %.cpp
#	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

$(PWIZ_OBJ_ARCH)/%.o: %.cpp
	@mkdir -v -p $(dir $@)
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(TEST_OBJS) $(PWIZ_OBJS) $(MAIN)

depend: $(SRCS)
