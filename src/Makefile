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

OBJ_DIR = ../obj

BIN_DIR = ../bin

ZERO_BASE_OBJS = $(addprefix $(OBJ_DIR)/, support_peak_type.o semi_align_type.o db_residue_seq.o extreme_value.o algorithm.o ion.o residue_seq.o fasta_reader.o change.o break_point.o bp_spec.o proteoform.o base_data.o activation.o ion_type.o neutral_loss.o trunc.o prot_mod.o acid.o ptm.o residue.o xml_dom.o xml_dom_document.o)

ZERO_SPEC_OBJS = $(addprefix $(OBJ_DIR)/, theo_peak.o  support_peak.o extend_peak.o prm_peak.o spectrum_set.o sp_para.o deconv_ms.o ms_header.o deconv_peak.o msalign_reader.o extend_sp_para.o peak_tolerance.o)

ZERO_PRSM_OBJS = $(addprefix $(OBJ_DIR)/, peak_ion_pair.o prsm.o prsm_writer.o)

ZERO_PTM_OBJS = $(addprefix $(OBJ_DIR)/, zero_ptm_slow_match.o zero_ptm_fast_match.o zero_ptm_search.o zero_ptm_console.o)

ZERO_OBJS = $(ZERO_BASE_OBJS) $(ZERO_SPEC_OBJS) $(ZERO_PRSM_OBJS) $(ZERO_PTM_OBJS)


#
# define the executable file 
ZERO = $(BIN_DIR)/zero_ptm
#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
zero:    $(ZERO)

$(ZERO): $(ZERO_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(ZERO) $(ZERO_OBJS) $(LFLAGS) $(LIBS)


all:
	cd base; make
	cd spec; make
	cd prsm; make
	cd zeroptmsearch; make
	cd tdgf; make
	cd filterdiagonal; make
	cd ptmsearch; make
	#cd console; make
	cd xpp; make	
	cd test; make

clean:
	cd base; make clean
	cd spec; make clean
	cd prsm; make clean
	cd zeroptmsearch; make clean
	cd tdgf; make clean
	cd filterdiagonal; make clean
	cd ptmsearch; make clean
	#cd console; make clean
	cd xpp; make clean
	cd test; make clean

