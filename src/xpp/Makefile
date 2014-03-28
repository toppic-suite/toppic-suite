# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../ -I../../lib/xerces-c/include -I../../lib/xalan-c/c/src -I../../lib/xalan-c/c/nls/include
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib -L/opt/local/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c 
#
# define the C++ source files
SRCS = xml_generator.cpp transformer.cpp 
#
# define the C object files 
#

OBJ_DIR = ../../obj

XPP_OBJS = $(addprefix $(OBJ_DIR)/, xml_generator.o transformer.o) 

OBJS = $(XPP_OBJS)

BIN_DIR = ../../bin


#
# define the executable file 
#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
all:    $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(XPP_OBJS)

depend: $(SRCS)
