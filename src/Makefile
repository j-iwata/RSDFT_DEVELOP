include make.inc

########################################################################
########################################################################

include makefile.common

DIR1 = ext1
EXTOBJ1 = $(DIR1)/ext_sub_minpac.o

DIR2 = ext2
EXTOBJ2 = $(DIR2)/dotp.o \
          $(DIR2)/convert_capital.o \
          $(DIR2)/write_info.o \
          $(DIR2)/stop_program.o \

DIR4 = FFTE
FFTOBJ = $(DIR4)/pzfft3dv.o \
         $(DIR4)/fft235.o \
         $(DIR4)/kernel.o \
         $(DIR4)/zfft3d.o \
         $(DIR4)/pzfft3d.o \

########################################################################
########################################################################

.PHONY: all clean re test runtest cleartest

all :
	@$(MAKE) lda0
	cd $(DIR1) ; $(MAKE)
	cd $(DIR2) ; $(MAKE)
	cd $(DIR4) ; $(MAKE)
	@$(MAKE) realspace.o
	$(FC) $(LFLAGS) $(EXTOBJ1) $(EXTOBJ2) $(MINPACOBJ) $(MDOBJ) $(FFTOBJ) $(MODS1) realspace.o $(FFTW_L) $(LAPACK_L) $(Eigen_L) -o rsdft.x

lda0 : $(MODS1)

include makefile.common.program
#include makefile.common.dep
#include makefile.test

clean :
	rm -f *.o *.mod mpif.h *.lst *.optlog *.i90
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR4) ; $(MAKE) clean

