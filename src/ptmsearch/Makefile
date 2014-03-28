# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../ -I../../lib/xerces-c/include
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib -L/opt/local/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c 
#
# define the C++ source files
SRCS = main.cpp basic_diag_pair.cpp diagonal.cpp dp_pair.cpp ptm_main.cpp ptm_processor.cpp ptm_slow_filter.cpp comp_shift_low_mem.cpp diagonal_header.cpp ps_align.cpp ptm_mng.cpp ptm_searcher.cpp ptm_slow_match.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj
#OBJS = $(addprefix $(OBJ_DIR)/, basic_diag_pair.o diagonal.o dp_pair.o ptm_main.o ptm_processor.o ptm_slow_filter.o comp_shift_low_mem.o diagonal_header.o ps_align.o ptm_mng.o ptm_searcher.o ptm_slow_match.o)

OBJS = $(addprefix $(OBJ_DIR)/, ptm_processor.o ptm_slow_filter.o ptm_slow_match.o ps_align.o basic_diag_pair.o dp_pair.o diagonal.o diagonal_header.o ptm_mng.o comp_shift_low_mem.o )

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
