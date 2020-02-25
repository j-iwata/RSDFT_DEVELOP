include ../make.inc

include makefile.common

FFLAGS := -I../xc -I../kinetic -I../hartree -I../bc -I../mixing $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst
