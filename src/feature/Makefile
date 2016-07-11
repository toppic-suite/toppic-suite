include ../Makefile.in

OBJS = $(FEATURE_OBJS) 

LINK_OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(addprefix $(OBJ_DIR)/, local_anno.o)

#
# define the executable file 
MAIN = $(BIN_DIR)/test_feature
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
	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LINK_OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@


clean:
	$(RM) *~ $(MAIN) $(OBJS)

depend: $(SRCS)
