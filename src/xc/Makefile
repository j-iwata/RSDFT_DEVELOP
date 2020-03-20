include ../make.inc

include makefile.common

FFLAGS := -I../hartree -I../bc -I../kinetic $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.L *.LL *.optlog *.i90
