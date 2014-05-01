include ../var.mk

TEST_OBJS = $(addprefix $(OBJ_DIR)/, argument.o test_console.o)

OBJS = $(BASE_OBJS) $(FILE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(ZERO_PTM_OBJS) $(FILTER_OBJS) $(PTM_OBJS) $(TDGF_OBJS) $(XPP_OBJS) $(TEST_OBJS)

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
	$(CC)  $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(TEST_OBJS)

depend: $(SRCS)
