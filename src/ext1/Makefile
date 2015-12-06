include ../make.inc

#OBJ = ext_sub_blas.o ext_sub_lapack_aux.o ext_sub_lapack_r8.o\
      ext_sub_blas_tmp.o ext_sub_lapack_c16_2.o ext_sub_lapack_u.o\
      ext_sub_eispack.o ext_sub_lapack_c16.o ext_sub_minpac.o
OBJ = ext_sub_minpac.o

all : $(OBJ)

#ext_sub_lapack_aux.o : ext_sub_lapack_aux.f
#	$(FC) -c -CB ext_sub_lapack_aux.f

clean :
	rm -r -f *.o mpif.h *.lst
