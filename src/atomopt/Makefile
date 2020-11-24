include ../make.inc

include makefile.common

FFLAGS := -I../kinetic -I../xc -I../hartree $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst
