# define the g++ compiler to use
CC = g++
#
# define any compile-time flags
CPPFLAGS = -Wall -std=c++11 -O3 -pg
FILECPPFLAGS = -Wall -O3 -pg
#
# define any directories containing header files other than /usr/include
#
INCLUDES = -I../../thirdparty/linux_include  -I../ 
#
# define library paths in addition to /usr/lib
LFLAGS = -L../../thirdparty/linux_lib -L/lib -L/opt/local/lib
#
# define any libraries to link into executable:
LIBS = -static -lxalan-c -lxalanMsg -lxerces-c -lboost_program_options -lboost_filesystem -lboost_system -lpthread

OBJ_DIR = ../../obj

BIN_DIR = ../../bin

BASE_OBJS = $(addprefix $(OBJ_DIR)/, file_util.o anno_residue.o species.o xml_writer.o support_peak_type.o semi_align_type.o db_residue_seq.o extreme_value.o algorithm.o ion.o base_data.o fasta_reader.o proteoform.o bp_spec.o break_point.o residue_seq.o change.o acid.o ptm.o residue_freq.o residue.o xml_dom.o xml_dom_document.o ion_type.o neutral_loss.o trunc.o prot_mod.o activation.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, peak.o extend_peak.o support_peak.o prm_peak.o spectrum_set.o sp_para.o theo_peak.o extend_sp_para.o peak_tolerance.o msalign_reader.o deconv_ms.o peak_list.o ms_header.o deconv_peak.o)

PRSM_OBJS = $(addprefix $(OBJ_DIR)/, prsm_fdr.o prsm_species.o output_selector.o table_writer.o prsm_selector.o cleavage.o pair.o prsm_combine.o prsm_writer.o simple_prsm_writer.o peak_ion_pair.o prsm.o simple_prsm.o)

ZERO_PTM_OBJS = $(addprefix $(OBJ_DIR)/, zero_ptm_slow_match.o zero_ptm_fast_match.o zero_ptm_search.o)

FILTER_OBJS =  $(addprefix $(OBJ_DIR)/, comp_shift_hi_mem.o ptm_fast_filter_block.o ptm_fast_filter_hi_mem.o ptm_fast_filter_processor.o) 

PTM_OBJS = $(addprefix $(OBJ_DIR)/, ptm_processor.o ptm_slow_filter.o ptm_slow_match.o ps_align.o basic_diag_pair.o dp_pair.o diagonal.o diagonal_header.o comp_shift_low_mem.o )

TDGF_OBJS = $(addprefix $(OBJ_DIR)/, evalue_processor.o comp_prob_value.o count_test_num.o comp_pvalue_array.o) 

XPP_OBJS = $(addprefix $(OBJ_DIR)/, anno_view.o xml_generator.o transformer.o) 
