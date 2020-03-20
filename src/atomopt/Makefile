include ../make.inc

include makefile.common

FFLAGS := -I../kinetic -I../xc -I../hartree -I../bc -I../mixing -I../subspace_diag -I../gram_schmidt -I../conjugate_gradient $(FFLAGS)

all : $(OBJS)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.L *.LL *.optlog *.i90
