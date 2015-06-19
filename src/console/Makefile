include ../Makefile.in

TOPPIC_OBJS = $(addprefix $(OBJ_DIR)/, toppic.o)

TOPPIC_STAT_OBJS = $(addprefix $(OBJ_DIR)/, toppic_stat.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_OBJS) $(DIAG_FILTER_OBJS) $(ONE_PTM_FILTER_OBJS) $(PTM_OBJS) $(TDGF_OBJS) $(PRSM_VIEW_OBJS) $(ARGUMENT_OBJS) 


BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/toppic  $(BIN_DIR)/toppic_stat
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

$(MAIN): $(OBJS) $(TOPPIC_OBJS) $(TOPPIC_STAT_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(BIN_DIR)/toppic $(OBJS) $(TOPPIC_OBJS) $(LFLAGS) $(LIBS)
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(BIN_DIR)/top_stat $(OBJS) $(TOPPIC_STAT_OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(ARGUMENT_OBJS) $(TOPPIC_OBJS)

depend: $(SRCS)
