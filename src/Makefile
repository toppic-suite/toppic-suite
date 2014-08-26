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
	cd zeroptmsearch; make
#	cd poisson; make
	cd tdgf; make
	cd filterdiagonal; make
	cd ptmsearch; make
#	cd xpp; make	
#	cd console; make
	cd test; make

clean:
	cd base; make clean
	cd spec; make clean
	cd prsm; make clean
	cd poisson; make clean
	cd zeroptmsearch; make clean
	cd tdgf; make clean
	cd filterdiagonal; make clean
	cd ptmsearch; make clean
	cd xpp; make clean
	cd console; make clean
	cd test; make clean
