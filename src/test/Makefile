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
LIBS = -lxerces-c-3.1 
#-lxalan-c -lxalanMsg 
#
# define the C++ source files
SRCS = test_console.cpp 
#
# define the C object files 
#

OBJ_DIR = ../../obj
BASE_OBJS = $(addprefix $(OBJ_DIR)/, anno_residue.o xml_writer.o species.o ion.o algorithm.o residue_freq.o support_peak_type.o semi_align_type.o db_residue_seq.o extreme_value.o residue_seq.o fasta_reader.o change.o break_point.o bp_spec.o proteoform.o base_data.o activation.o ion_type.o neutral_loss.o trunc.o prot_mod.o acid.o ptm.o residue.o xml_dom.o xml_dom_document.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, theo_peak.o spectrum_set.o msalign_reader.o deconv_ms.o support_peak.o deconv_peak.o extend_peak.o prm_peak.o ms_header.o sp_para.o extend_sp_para.o peak_tolerance.o) 

PRSM_OBJS = $(addprefix $(OBJ_DIR)/, prsm_fdr.o prsm_species.o output_selector.o table_writer.o prsm_selector.o prsm_combine.o simple_prsm_writer.o simple_prsm.o cleavage.o peak_ion_pair.o pair.o prsm.o prsm_writer.o)

ZERO_PTM_OBJS = $(addprefix $(OBJ_DIR)/, zero_ptm_slow_match.o zero_ptm_fast_match.o zero_ptm_search.o)

FILTER_OBJS = $(addprefix $(OBJ_DIR)/, comp_shift_hi_mem.o ptm_fast_filter_block.o ptm_fast_filter_hi_mem.o ptm_fast_filter_processor.o) 

PTM_OBJS = $(addprefix $(OBJ_DIR)/, ptm_processor.o ptm_slow_filter.o ptm_slow_match.o ps_align.o basic_diag_pair.o dp_pair.o diagonal.o diagonal_header.o comp_shift_low_mem.o )

TDGF_OBJS = $(addprefix $(OBJ_DIR)/, evalue_processor.o comp_prob_value.o count_test_num.o comp_pvalue_array.o) 

#XPP_OBJS = $(addprefix $(OBJ_DIR)/, anno_view.o xml_generator.o transformer.o) 

TEST_OBJS = $(addprefix $(OBJ_DIR)/, test_console.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_OBJS) $(FILTER_OBJS) $(PTM_OBJS) $(TDGF_OBJS) $(XPP_OBJS) $(TEST_OBJS)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/test_console
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
	$(RM) *~ $(MAIN) $(TEST_OBJS)

depend: $(SRCS)
