include ../Makefile.in

TOPPIC_OBJS = $(addprefix $(OBJ_DIR)/, toppic.o summary.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_FILTER_OBJS) $(ZERO_PTM_SEARCH_OBJS) $(ONE_PTM_FILTER_OBJS) $(ONE_PTM_SEARCH_OBJS) $(DIAG_FILTER_OBJS) $(PTM_SEARCH_OBJS) $(TDGF_OBJS) $(PRSM_VIEW_OBJS) $(LOCAL_OBJS) $(ARGUMENT_OBJS) $(TOPPIC_OBJS)


BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/toppic
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
	$(RM) *~ $(MAIN) $(ARGUMENT_OBJS) $(TOPPIC_OBJS)

depend: $(SRCS)
