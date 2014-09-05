include ../Makefile.in

OBJS = $(BASE_OBJS) $(SPEC_OBJS) $(PRSM_OBJS) $(PRSM_VIEW_OBJS)


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
all: $(OBJS)

#	$(CC) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) *~ $(PRSM_VIEW_OBJS)

depend: $(SRCS)
