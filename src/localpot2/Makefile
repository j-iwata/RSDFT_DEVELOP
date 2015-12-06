include ../make.inc

OBJ = localpot2_variables.o\
      localpot2_ion_module.o\
      localpot2_density_module.o\
      localpot2_vh_module.o\
      localpot2_xc_module.o\
      localpot2_te_module.o\
      localpot2_module.o\
      localpot2_Smatrix_module.o\
      gram_schmidt_u_module.o\
      cg_u_module.o\

all : $(OBJ)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.optlog *.i90

