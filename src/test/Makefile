include ../Makefile.in

TEST_ZERO_PTM_OBJS = $(addprefix $(OBJ_DIR)/, test_zero_ptm.o)

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_OBJS) $(FILTER_OBJS) $(PTM_OBJS) $(TDGF_OBJS) $(POISSON_OBJS) $(XPP_OBJS) $(CONSOLE_OBJS) 

TEST_ZERO_PTM_LINK_OBJS = $(OBJS) $(TEST_ZERO_PTM_OBJS)

BIN_DIR = ../../bin

#
# define the executable file 
TEST_ZERO_PTM = $(BIN_DIR)/test_zero_ptm
#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
all:    $(TEST_ZERO_PTM)

$(TEST_ZERO_PTM): $(TEST_ZERO_PTM_LINK_OBJS) 
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(TEST_ZERO_PTM) $(TEST_ZERO_PTM_LINK_OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(TEST_ZERO_PTM) $(TEST_ZERO_PTM_OBJS)

depend: $(SRCS)
