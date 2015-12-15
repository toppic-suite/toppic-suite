#
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
#
.PHONY: depend clean

all:
	cd base; make
	cd spec; make
	cd prsm; make
	cd zeroptmfilter; make
	cd zeroptmsearch; make
#	cd oneptmfilter; make
#	cd diagfilter; make
#	cd ptmsearch; make
	cd tdgf; make
#	cd prsmview; make	
	cd console; make
#	cd graph; make
#	cd graphalign; make
	cd test; make

clean:
	cd base; make clean
	cd spec; make clean
	cd prsm; make clean
	cd zeroptmfilter; make clean
	cd zeroptmsearch; make clean
	cd oneptmfilter; make clean
	cd diagfilter; make clean
	cd ptmsearch; make clean
	cd tdgf; make clean
	cd prsmview; make clean
	cd console; make clean
	cd test; make clean
