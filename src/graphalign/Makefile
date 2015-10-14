include ../Makefile.in

LINK_OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(PTM_OBJS) $(TDGF_OBJS) $(ARGUMENT_OBJS) $(GRAPH_OBJS)

OBJS = $(GRAPH_ALIGN_OBJS) 

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/graph_align
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

$(MAIN): $(OBJS) $(LINK_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(BIN_DIR)/graph_align $(OBJS) $(LINK_OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(GRAPH_ALIGN_OBJS)

depend: $(SRCS)
