# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../base
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c -llog4cxx
#
# define the C++ source files
SRCS = msalign_reader.cpp deconv_ms.cpp ms.cpp peak_list.cpp ms_header.cpp deconv_peak.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj
OBJS = $(addprefix $(OBJ_DIR)/, ion_type.o xml_dom.o xml_dom_document.o activation.o msalign_reader.o deconv_ms.o ms.o peak_list.o ms_header.o deconv_peak.o)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/spec
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

$(MAIN): $(OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(OBJS)

depend: $(SRCS)
