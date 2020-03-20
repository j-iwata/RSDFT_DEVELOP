include ../make.inc

include makefile.common

FFLAGS := -I../xc -I../kinetic -I../mixing -I../bc -I../hartree $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.L *.LL *.optlog *.i90
