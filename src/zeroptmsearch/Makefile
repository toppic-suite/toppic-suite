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
SRCS = zero_ptm_slow_match.cpp zero_ptm_search.cpp zero_ptm_console.cpp zero_ptm_fast_match.cpp 
#
# define the C object files 
#

OBJ_DIR = ../../obj
BASE_OBJS = $(addprefix $(OBJ_DIR)/, db_residue_seq.o extreme_value.o algorithm.o ion.o residue_seq.o fasta_reader.o change.o break_point.o bp_spec.o proteoform.o base_data.o activation.o ion_type.o neutral_loss.o trunc.o prot_mod.o acid.o ptm.o residue.o xml_dom.o xml_dom_document.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, theo_peak.o support_peak_type.o support_peak.o extend_peak.o prm_peak.o spectrum_set.o sp_para.o deconv_ms.o ms_header.o deconv_peak.o msalign_reader.o extend_sp_para.o peak_tolerance.o)

PRSM_OBJS = $(addprefix $(OBJ_DIR)/, peak_ion_pair.o prsm.o prsm_writer.o)

ZERO_PTM_OBJS = $(addprefix $(OBJ_DIR)/, zero_ptm_slow_match.o zero_ptm_fast_match.o zero_ptm_search.o zero_ptm_console.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_OBJS)

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
	$(RM) *~ $(MAIN) $(ZERO_PTM_OBJS)

depend: $(SRCS)
