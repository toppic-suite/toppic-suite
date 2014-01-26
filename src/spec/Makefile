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
SRCS = theo_peak.cpp extend_sp_para.cpp peak_tolerance.cpp msalign_reader.cpp deconv_ms.cpp  peak_list.cpp ms_header.cpp deconv_peak.cpp
#
# define the C object files 
#

OBJ_DIR = ../../obj
BASE_OBJS = $(addprefix $(OBJ_DIR)/, support_peak_type.o neutral_loss.o semi_align_type.o prot_mod.o db_residue_seq.o trunc.o change.o break_point.o bp_spec.o acid.o ptm.o residue.o residue_seq.o proteoform.o ion.o ion_type.o xml_dom.o xml_dom_document.o activation.o)

SPEC_OBJS = $(addprefix $(OBJ_DIR)/, extend_peak.o support_peak.o prm_peak.o spectrum_set.o sp_para.o theo_peak.o extend_sp_para.o peak_tolerance.o msalign_reader.o deconv_ms.o peak_list.o ms_header.o deconv_peak.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/spec
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
	$(RM) *~ $(MAIN) $(SPEC_OBJS)

depend: $(SRCS)
