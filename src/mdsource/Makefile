include ../make.inc

FFLAGS := -I../hartree -I../xc -I../kinetic $(FFLAGS)

MOBJ = cpmd_variables.o

OBJ = cpmdio_module.o\
      calc_kine_temp_module.o\
      nose_hoover_chain_ele_module.o\
      nose_hoover_chain_module.o\
      berendsen_module.o\
      velocity_scaling_module.o\
      velocity_scaling_ele_module.o\
      setv.o\
      vcom.o\
      mdio.o\
      active_band.o\
      alloc_cpmd.o\
      calc_overlap_module.o \
      overlap_cpmd_module.o \
      rotorb_module.o \
      blue_moon_module.o \
      bomd.o\
      getforce_cpmd.o\
      getforce.o\
      rotorb2.o\
      mdio.o\
      calfke.o\
      wf_force.o\

all : $(MOBJ) $(OBJ)

clean :
	rm -r -f *.o *.mod mpif.h *.lst *.optlog *.i90

