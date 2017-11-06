include make.inc

########################################################################
########################################################################

include makefile.common

DIR0 = ext0
include $(DIR0)/makefile.common
OBJ0 := $(addprefix $(DIR0)/,$(OBJS))

DIR1 = ext1
include $(DIR1)/makefile.common
OBJ1 := $(addprefix $(DIR1)/,$(OBJS))

DIR2 = ext2
include $(DIR2)/makefile.common
OBJ2 := $(addprefix $(DIR2)/,$(OBJS))

DIR_GS = gram_schmidt
include $(DIR_GS)/makefile.common
OBJ_GS := $(addprefix $(DIR_GS)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_GS) $(INCLUDE_DIR)

DIR_SD = subspace_diag
include $(DIR_SD)/makefile.common
OBJ_SD := $(addprefix $(DIR_SD)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_SD) $(INCLUDE_DIR)

DIR_CG = conjugate_gradient
include $(DIR_CG)/makefile.common
OBJ_CG := $(addprefix $(DIR_CG)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_CG) $(INCLUDE_DIR)

DIR_HT = hartree
include $(DIR_HT)/makefile.common
OBJ_HT := $(addprefix $(DIR_HT)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_HT) $(INCLUDE_DIR)

DIR_XC = xc
include $(DIR_XC)/makefile.common
OBJ_XC := $(addprefix $(DIR_XC)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_XC) $(INCLUDE_DIR)

DIR_KI = kinetic
include $(DIR_KI)/makefile.common
OBJ_KI := $(addprefix $(DIR_KI)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_KI) $(INCLUDE_DIR)

DIR3 = mdsource
MDOBJ = $(DIR3)/cpmdio_module.o \
        $(DIR3)/calc_kine_temp_module.o\
        $(DIR3)/nose_hoover_chain_ele_module.o\
        $(DIR3)/nose_hoover_chain_module.o\
        $(DIR3)/berendsen_module.o\
        $(DIR3)/velocity_scaling_module.o\
        $(DIR3)/velocity_scaling_ele_module.o\
        $(DIR3)/active_band.o \
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

.PHONY: all clean

all :
	@$(MAKE) lda0
	cd $(DIR0) ; $(MAKE)
	cd $(DIR1) ; $(MAKE)
	cd $(DIR2) ; $(MAKE)
	cd $(DIR_KI) ; $(MAKE)
	cd $(DIR_HT) ; $(MAKE)
	cd $(DIR_XC) ; $(MAKE)
	@$(MAKE) lda2
	cd $(DIR_GS) ; $(MAKE)
	cd $(DIR_SD) ; $(MAKE)
	cd $(DIR_CG) ; $(MAKE)
	@$(MAKE) lda1
	cd $(DIR3) ; $(MAKE)
	@$(MAKE) main.o
	$(FC) $(LFLAGS) $(OBJ0) $(OBJ1) $(OBJ2)\
                        $(MODS1) $(OBJ_XC) $(OBJ_HT)\
                        $(OBJ_KI)\
                        $(MODS2)\
                        $(OBJ_GS) $(OBJ_SD) $(OBJ_CG)\
                        $(MODS0) $(MDOBJ) main.o $(LIBXC_L) $(FFTW_L) $(LAPACK_L) -o rsdft.x

lda0 : $(MODS1)
lda1 : $(MODS0)
lda2 : $(MODS2)

clean :
	rm -f *.o *.mod mpif.h *.lst *.optlog *.i90
	cd $(DIR0) ; $(MAKE) clean
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR3) ; $(MAKE) clean
	cd $(DIR_GS) ; $(MAKE) clean
	cd $(DIR_SD) ; $(MAKE) clean
	cd $(DIR_CG) ; $(MAKE) clean
	cd $(DIR_HT) ; $(MAKE) clean
	cd $(DIR_XC) ; $(MAKE) clean
	cd $(DIR_KI) ; $(MAKE) clean
