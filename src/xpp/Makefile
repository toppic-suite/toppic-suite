include ../Makefile.in

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(XPP_OBJS)

MAIN = $(BIN_DIR)/xpp

#
# define the executable file 
#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean
#
#
#
all: $(MAIN)

$(MAIN): $(OBJS) 
#	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(MAIN) $(XPP_OBJS)

depend: $(SRCS)
