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

DIR_BC = bc
include $(DIR_BC)/makefile.common
OBJ_BC := $(addprefix $(DIR_BC)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_BC) $(INCLUDE_DIR)

DIR_KI = kinetic
include $(DIR_KI)/makefile.common
OBJ_KI := $(addprefix $(DIR_KI)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_KI) $(INCLUDE_DIR)

DIR_MIX = mixing
include $(DIR_MIX)/makefile.common
OBJ_MIX := $(addprefix $(DIR_MIX)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_MIX) $(INCLUDE_DIR)

DIR_AO = atomopt
include $(DIR_AO)/makefile.common
OBJ_AO := $(addprefix $(DIR_AO)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_AO) $(INCLUDE_DIR)

DIR_MD = mdsource
include $(DIR_MD)/makefile.common
OBJ_MD := $(addprefix $(DIR_MD)/,$(OBJS))
INCLUDE_DIR := -I$(DIR_MD) $(INCLUDE_DIR)

########################################################################
########################################################################

.PHONY: all clean

all :
	./git_info.sh
	@$(MAKE) lda0
	cd $(DIR0) ; $(MAKE)
	cd $(DIR1) ; $(MAKE)
	cd $(DIR2) ; $(MAKE)
	cd $(DIR_BC) ; $(MAKE)
	cd $(DIR_KI) ; $(MAKE)
	cd $(DIR_HT) ; $(MAKE)
	cd $(DIR_XC) ; $(MAKE)
	cd $(DIR_MIX) ; $(MAKE)
	@$(MAKE) lda2
	cd $(DIR_GS) ; $(MAKE)
	cd $(DIR_SD) ; $(MAKE)
	cd $(DIR_CG) ; $(MAKE)
	@$(MAKE) lda1
	cd $(DIR_AO) ; $(MAKE)
	@$(MAKE) lda3
	cd $(DIR_MD) ; $(MAKE)
	@$(MAKE) main.o
	$(FC) $(LFLAGS) $(OBJ0) $(OBJ1) $(OBJ2)\
                        $(MODS1) $(OBJ_XC) $(OBJ_HT) $(OBJ_BC)\
                        $(OBJ_KI) $(OBJ_MIX) $(OBJ_IO)\
                        $(MODS2)\
                        $(OBJ_GS) $(OBJ_SD) $(OBJ_CG)\
                        $(MODS0)\
                        $(OBJ_AO)\
                        $(MODS3)\
                        $(OBJ_MD) main.o $(LIBXC_L) $(FFTW_L) $(LAPACK_L) -o rsdft.x

lda0 : $(MODS1)
lda1 : $(MODS0)
lda2 : $(MODS2)
lda3 : $(MODS3)

clean :
	rm -f *.o *.mod mpif.h *.lst *.optlog *.i90
	cd $(DIR0) ; $(MAKE) clean
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR_MD) ; $(MAKE) clean
	cd $(DIR_GS) ; $(MAKE) clean
	cd $(DIR_SD) ; $(MAKE) clean
	cd $(DIR_CG) ; $(MAKE) clean
	cd $(DIR_HT) ; $(MAKE) clean
	cd $(DIR_XC) ; $(MAKE) clean
	cd $(DIR_KI) ; $(MAKE) clean
	cd $(DIR_AO) ; $(MAKE) clean
	cd $(DIR_MIX) ; $(MAKE) clean
	cd $(DIR_BC) ; $(MAKE) clean
