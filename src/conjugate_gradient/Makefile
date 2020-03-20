include ../make.inc

include makefile.common

FFLAGS := -I../xc -I../bc -I../kinetic -I../hartree -I../mixing $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.L *.LL *.optlog *.i90
