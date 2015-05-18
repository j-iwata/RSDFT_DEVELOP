include make.inc

########################################################################
########################################################################
include makefile.common

DIR1 = ext1
MINPACOBJ  = $(DIR1)/ext_sub_minpac.o

DIR2 = ext2
EXTOBJ2 = $(DIR2)/polint.o \
          $(DIR2)/indexx.o \
          $(DIR2)/bberf.o \
          $(DIR2)/p4sn.o \
          $(DIR2)/ylm.o \
          $(DIR2)/dotp.o \
          $(DIR2)/spline.o \
          $(DIR2)/convert_capital.o \
          $(DIR2)/gaussj.o \
          $(DIR2)/write_info.o \

DIR3 = mdsource
MDOBJ = $(DIR3)/cpmdio2_module.o \
        $(DIR3)/active_band.o \
        $(DIR3)/asign_atom.o \
        $(DIR3)/init_ion.o \
        $(DIR3)/alloc_cpmd.o \
        $(DIR3)/blue_moon_module.o \
        $(DIR3)/bomd.o \
        $(DIR3)/getforce.o \
        $(DIR3)/getforce_cpmd.o \
        $(DIR3)/wf_force.o \
        $(DIR3)/calkin.o \
        $(DIR3)/calfke.o \
        $(DIR3)/calc_overlap_module.o \
        $(DIR3)/overlap_cpmd_module.o \
        $(DIR3)/rotorb_module.o \
        $(DIR3)/rotorb2.o \
        $(DIR3)/setv.o \
        $(DIR3)/vcom.o \
        $(DIR3)/mdio.o \
        $(DIR3)/mnhc.o \
        $(DIR3)/cpmd_variables.o \

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
	cd $(DIR3) ; $(MAKE) -j1
	$(FC) $(FFLAGS) $(EXTOBJ2) $(MINPACOBJ) $(MDOBJ) $(FFTOBJ) $(LAPACK_L) $(MODS1) realspace.f90 -o realspace.x

lda0 : $(MODS1)

re:
	$(MAKE) -f makefile.simple

include makefile.common.program
include makefile.common.dep
include makefile.test

clean :
	rm -f *.o *.mod a.out mpif.h *.lst *.x *.optlog *.i90
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR3) ; $(MAKE) clean
	cd $(DIR4) ; $(MAKE) clean

