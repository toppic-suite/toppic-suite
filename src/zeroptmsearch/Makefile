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
LFLAGS = -L/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c -llog4cxx
#
# define the C++ source files
SRCS = zero_ptm_mng.cpp zero_ptm_console.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj
OBJS = $(addprefix $(OBJ_DIR)/, deconv_ms.o ms.o ms_header.o deconv_peak.o  msalign_reader.o residue_seq.o fasta_reader.o change.o break_point.o bp_spec.o proteoform.o base_data.o activation.o ion_type.o neutral_loss.o trunc.o prot_mod.o acid.o ptm.o residue.o zero_ptm_console.o zero_ptm_mng.o extend_sp_para.o peak_tolerance.o xml_dom.o xml_dom_document.o)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/zero_ptm
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
