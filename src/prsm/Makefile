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
SRCS = prsm_writer.cpp peak_ion_pair.cpp prsm.cpp simple_prsm_writer.cpp simple_prsm.cpp main.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj

BASE_OBJS = $(addprefix $(OBJ_DIR)/, neutral_loss.o semi_align_type.o db_residue_seq.o prot_mod.o trunc.o activation.o extreme_value.o algorithm.o break_point.o bp_spec.o ion.o ion_type.o change.o ptm.o acid.o residue.o residue_seq.o proteoform.o xml_dom_document.o xml_dom.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, theo_peak.o deconv_peak.o extend_sp_para.o  sp_para.o deconv_ms.o peak_tolerance.o extend_peak.o ms_header.o)

PRSM_OBJS = $(addprefix $(OBJ_DIR)/, main.o prsm_writer.o simple_prsm_writer.o peak_ion_pair.o prsm.o simple_prsm.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/prsm
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
