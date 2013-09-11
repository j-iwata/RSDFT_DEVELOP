include ../make.inc

OBJ = init_mesh.o init_parallel_mol.o init_fdparallel_mol.o\
      init_ps_mol.o prep_ps_mol.o prep_ps_sub_mol.o calc_Eion.o init_rho_mol.o\
      hpsi_mol.o hartree_mol.o prep_hartree.o calc_vloc_rho.o\
      force_ion_mol.o force_local_mol.o force_IonIon.o \
      init_mesh_sym.o init_parallel_mol_sym.o prep_ps_mol_sym.o prep_hartree_sym.o\
      hartree_mol_sym.o calc_vloc_rho_sym.o init_wf_mol_sym.o d_gsmol_sym.o\
      hpsi_mol_sym.o d_diag_mol_sym.o d_dtcg_mol_sym.o precond_cg_mol_sym.o e_tot_mol_sym.o\
      re_sym.o fock.o

all : $(OBJ)
	cp ../global_variables.mod ./

clean :
	rm -r -f *.o mpif.h *.mod *.lst

