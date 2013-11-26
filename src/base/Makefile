# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../../include/base
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c -llog4cxx
#
# define the C++ source files
SRCS = activation.cpp neutral_loss.cpp fasta_reader.cpp residue_seq.cpp prot_mod.cpp acid.cpp component.cpp trunc.cpp ion_type.cpp ptm.cpp residue.cpp xml_dom.cpp xml_dom_document.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj
OBJS = $(addprefix $(OBJ_DIR)/, activation.o neutral_loss.o fasta_reader.o residue_seq.o prot_mod.o acid.o component.o trunc.o ion_type.o ptm.o residue.o xml_dom.o xml_dom_document.o)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/prot
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
