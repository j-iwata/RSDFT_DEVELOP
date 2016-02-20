include make.inc

########################################################################
########################################################################

include makefile.common

DIR0 = ext0
EXTOBJ0 = $(DIR0)/dotp.o \
          $(DIR0)/convert_capital.o \
          $(DIR0)/write_info.o \
          $(DIR0)/stop_program.o \
          $(DIR0)/get_atomic_number.o \
          $(DIR0)/exit_program.o \

DIR1 = ext1
EXTOBJ1 = $(DIR1)/ext_sub_minpac.o

DIR2 = ext2
EXTOBJ2 = $(DIR2)/pzfft3dv.o \
          $(DIR2)/factor.o \
          $(DIR2)/fft235.o \
          $(DIR2)/kernel.o \
          $(DIR2)/zfft3d.o \
          $(DIR2)/pzfft3d.o \

DIR3 = mdsource
MDOBJ = $(DIR3)/cpmdio2_module.o \
        $(DIR3)/calc_kine_temp_module.o\
        $(DIR3)/berendsen_module.o\
        $(DIR3)/active_band.o \
        $(DIR3)/init_ion.o \
        $(DIR3)/alloc_cpmd.o \
        $(DIR3)/blue_moon_module.o \
        $(DIR3)/bomd.o \
        $(DIR3)/getforce.o \
        $(DIR3)/getforce_cpmd.o \
        $(DIR3)/wf_force.o \
        $(DIR3)/calfke.o \
        $(DIR3)/calc_overlap_module.o \
        $(DIR3)/overlap_cpmd_module.o \
        $(DIR3)/rotorb_module.o \
        $(DIR3)/rotorb2.o \
        $(DIR3)/setv.o \
        $(DIR3)/vcom.o \
        $(DIR3)/mdio.o \
        $(DIR3)/cpmd_variables.o \

########################################################################
########################################################################

.PHONY: all clean re test runtest cleartest

all :
	@$(MAKE) lda0
	cd $(DIR0) ; $(MAKE)
	cd $(DIR1) ; $(MAKE)
	cd $(DIR2) ; $(MAKE)
	cd $(DIR3) ; $(MAKE)
	@$(MAKE) main.o
	$(FC) $(LFLAGS) $(EXTOBJ0) $(EXTOBJ1) $(EXTOBJ2) $(MDOBJ) $(MODS1) main.o $(FFTW_L) $(LAPACK_L) -o rsdft.x

lda0 : $(MODS1)

clean :
	rm -f *.o *.mod mpif.h *.lst *.optlog *.i90
	cd $(DIR0) ; $(MAKE) clean
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR3) ; $(MAKE) clean
