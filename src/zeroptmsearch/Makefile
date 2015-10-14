include ../Makefile.in

OBJS = $(ZERO_PTM_SEARCH_OBJS)

BIN_DIR = ../../bin

#
# define the executable file 
MAIN = $(BIN_DIR)/zero_ptm
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
#	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(ZERO_PTM_OBJS)

depend: $(SRCS)
