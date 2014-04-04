# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -g -std=c++11
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../ -I../../lib/xalan-c/c/src -I../../lib/xalan-c/c/nls/include
#
# define library paths in addition to /usr/lib
LFLAGS = -L/lib -L/opt/local/lib -L../../lib/xalan-c/c/lib
#
# define any libraries to link into executable:
LIBS = -lxerces-c-3.1 -lxalan-c -lxalanMsg 
#
# define the C++ source files
SRCS = xml_generator.cpp transformer.cpp 
#
# define the C object files 
#

OBJ_DIR = ../../obj
BASE_OBJS = $(addprefix $(OBJ_DIR)/, anno_residue.o xml_writer.o species.o ion.o algorithm.o residue_freq.o support_peak_type.o semi_align_type.o db_residue_seq.o extreme_value.o residue_seq.o fasta_reader.o change.o break_point.o bp_spec.o proteoform.o base_data.o activation.o ion_type.o neutral_loss.o trunc.o prot_mod.o acid.o ptm.o residue.o xml_dom.o xml_dom_document.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, theo_peak.o spectrum_set.o msalign_reader.o deconv_ms.o support_peak.o deconv_peak.o extend_peak.o prm_peak.o ms_header.o sp_para.o extend_sp_para.o peak_tolerance.o) 

PRSM_OBJS = $(addprefix $(OBJ_DIR)/, output_selector.o table_writer.o prsm_selector.o prsm_combine.o simple_prsm_writer.o simple_prsm.o cleavage.o peak_ion_pair.o pair.o prsm.o prsm_writer.o)

XPP_OBJS = $(addprefix $(OBJ_DIR)/, anno_view.o xml_generator.o transformer.o) 

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(XPP_OBJS)

BIN_DIR = ../../bin

MAIN = $(BIN_DIR)/xpp

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
#
#
all: $(MAIN)

$(MAIN): $(OBJS) 
#	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(XPP_OBJS)

depend: $(SRCS)
