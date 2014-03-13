# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I..
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib -L/opt/local/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c 
#
# define the C++ source files
SRCS = output_selector.cpp table_writer.cpp 
#
# define the C object files 
#

OBJ_DIR = ../../obj

OBJS = $(addprefix $(OBJ_DIR)/, output_selector.o table_writer.o)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/ptm_search
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
#	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(OBJS)

depend: $(SRCS)
