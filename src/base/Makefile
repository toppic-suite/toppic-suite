# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib -L/opt/local/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c 
#
# define the C++ source files
SRCS = main.cpp semi_align_type.cpp db_residue_seq.cpp extreme_value.cpp algorithm.cpp ion.cpp base_data.cpp fasta_reader.cpp proteoform.cpp bp_spec.cpp break_point.cpp residue_seq.cpp change.cpp acid.cpp ptm.cpp residue_freq.cpp residue.cpp xml_dom.cpp xml_dom_document.cpp ion_type.cpp neutral_loss.cpp trunc.cpp prot_mod.cpp activation.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj

OBJS = $(addprefix $(OBJ_DIR)/, main.o semi_align_type.o db_residue_seq.o extreme_value.o algorithm.o ion.o base_data.o fasta_reader.o proteoform.o bp_spec.o break_point.o residue_seq.o change.o acid.o ptm.o residue_freq.o residue.o xml_dom.o xml_dom_document.o ion_type.o neutral_loss.o trunc.o prot_mod.o activation.o)

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
