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
      bomd.o\
      calkin.o\
      getforce_cpmd.o\
      getforce.o\
      rotorb.o\
      rotorb2.o\
      overlap.o \
      mdio.o\
      calfke.o\
      wf_force.o\

#      cpmdio.o\

all : $(MOBJ) $(OBJ)

clean :
	rm -r -f *.o *.mod mpif.h *.lst

