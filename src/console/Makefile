include ../Makefile.in

TOPPIC_OBJS = $(addprefix $(OBJ_DIR)/, toppic.o summary.o)

TOPFD_OBJS= $(addprefix $(OBJ_DIR)/, topfd.o)

TOPPIC_LINK_OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_FILTER_OBJS) $(ZERO_PTM_SEARCH_OBJS) $(ONE_PTM_FILTER_OBJS) $(ONE_PTM_SEARCH_OBJS) $(DIAG_FILTER_OBJS) $(PTM_SEARCH_OBJS) $(TDGF_OBJS) $(PRSM_VIEW_OBJS) $(LOCAL_OBJS) $(ARGUMENT_OBJS) $(TOPPIC_OBJS)

TOPFD_LINK_OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(FEATURE_OBJS) $(TOPFD_OBJS) $(addprefix $(OBJ_DIR)/, local_anno.o)

BIN_DIR = ../../bin

#
# define the executable file 
TOPPIC = $(BIN_DIR)/toppic
TOPFD = $(BIN_DIR)/topfd
#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
all:    $(TOPPIC) $(TOPFD)

$(TOPPIC): $(TOPPIC_LINK_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(TOPPIC) $(TOPPIC_LINK_OBJS) $(LFLAGS) $(LIBS)

$(TOPFD): $(TOPFD_LINK_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(TOPFD) $(TOPFD_LINK_OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(TOPPIC) $(TOPFD) $(ARGUMENT_OBJS) $(TOPPIC_OBJS)

depend: $(SRCS)
