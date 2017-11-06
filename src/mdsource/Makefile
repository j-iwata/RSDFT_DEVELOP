include ../make.inc

include makefile.common

FFLAGS := -I../hartree -I../xc -I../kinetic $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.optlog *.i90

