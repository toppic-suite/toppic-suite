include ../Makefile.in

OBJS = $(FEATURE_OBJS) 
#
# define the executable file 
MAIN = $(BIN_DIR)/feature
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

$(MAIN): $(FILE_OBJS) $(OBJS)
# $(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)
#
$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@


clean:
	$(RM) *~ $(MAIN) $(OBJS)

depend: $(SRCS)
