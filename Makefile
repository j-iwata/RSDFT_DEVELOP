include make.inc

########################################################################
########################################################################

include makefile.common

DIR2 = ext2
EXTOBJ2 = $(DIR2)/p4sn.o \
          $(DIR2)/dotp.o \
          $(DIR2)/convert_capital.o \
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

DIR5 = localpot2
LPOT2OBJ = $(DIR5)/localpot2_variables.o\
           $(DIR5)/localpot2_module.o\
           $(DIR5)/localpot2_Smatrix_module.o\
           $(DIR5)/localpot2_density_module.o\
           $(DIR5)/localpot2_ion_module.o\
           $(DIR5)/localpot2_te_module.o\
           $(DIR5)/localpot2_xc_module.o\
           $(DIR5)/localpot2_vh_module.o\
           $(DIR5)/gram_schmidt_u_module.o\
           $(DIR5)/cg_u_module.o\

DIR6 = esm
ESMOBJ = $(ESMOBJ)/esm_cylindrical_testl.o\
         $(ESMOBJ)/esm_genpot_module.o\
         $(ESMOBJ)/esm_hartree_module.o\
         $(ESMOBJ)/esm_kinetic_module.o\
         $(ESMOBJ)/esm_rgrid_module.o\
         $(ESMOBJ)/esm_rshell_module.o\
         $(ESMOBJ)/ps_local_rs_module.o\

########################################################################
########################################################################

.PHONY: all clean re test runtest cleartest

all :
	@$(MAKE) lda0
	cd $(DIR2) ; $(MAKE)
	cd $(DIR4) ; $(MAKE)
	cd $(DIR3) ; $(MAKE) -j1
	cd $(DIR5) ; $(MAKE)
	cd $(DIR6) ; $(MAKE)
	@$(MAKE) realspace.o
	$(FC) $(LFLAGS) $(EXTOBJ2) $(MINPACOBJ) $(MDOBJ) $(FFTOBJ) $(LPOT2OBJ) $(LAPACK_L) $(MODS1) realspace.o $(LIBS) -o realspace.x

lda0 : $(MODS1)

re:
	$(MAKE) -f makefile.simple

include makefile.common.program
#include makefile.common.dep
#include makefile.test

clean :
	rm -f *.o *.mod a.out mpif.h *.lst *.x *.optlog *.i90
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR3) ; $(MAKE) clean
	cd $(DIR4) ; $(MAKE) clean
	cd $(DIR5) ; $(MAKE) clean
	cd $(DIR6) ; $(MAKE) clean

