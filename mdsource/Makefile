include ../make.inc

MOBJ = cpmd_variables.o

OBJ = mnhc.o\
      cpmdio2_module.o\
      setv.o\
      vcom.o\
      mdio.o\
      asign_atom.o\
      init_ion.o\
      active_band.o\
      alloc_cpmd.o\
      calc_overlap_module.o \
      overlap_cpmd_module.o \
      rotorb_module.o \
      blue_moon_module.o \
      bomd.o\
      calkin.o\
      getforce_cpmd.o\
      getforce.o\
      rotorb2.o\
      mdio.o\
      calfke.o\
      wf_force.o\

all : $(MOBJ) $(OBJ)

clean :
	rm -r -f *.o *.mod mpif.h *.lst

