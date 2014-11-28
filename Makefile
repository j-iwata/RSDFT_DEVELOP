include make.inc

MODS1 = simc_module.o\
        watch_module.o \
        info_module.o \
        modified_bessel_module.o \
	bberf_module.o \
        expint_module.o \
        parallel_module.o \
        rsdft_mpi_module.o \
        bcast_module.o \
        array_bound_module.o\
        ffte_sub_module.o \
        aa_module.o \
        bb_module.o \
        atom_module.o \
        rgrid_variables.o \
        rgrid_sol_module.o \
        rgrid_mol_module.o \
        esm_rgrid_module.o \
        rgrid_module.o \
        ggrid_module.o \
        bz_module.o \
        fd_module.o \
        bc_variables.o \
        bc_mol_module.o \
        bc_module.o \
        kinetic_variables.o  \
        kinetic_mol_module.o \
        esm_kinetic_module.o \
        tmp3_omp_variables.o \
        kinetic_fft_module.o \
        kinetic_sol_0_module.o \
        kinetic_sol_1_module.o \
        kinetic_module.o \
        strfac_module.o \
        ps_read_tm_module.o \
        ps_read_upf_module.o \
        ps_read_yb_module.o \
        ps_gth_module.o \
        pseudopot_module.o \
        ps_local_gth_module.o \
        electron_module.o\
        init_occ_electron_module.o \
        wf_sub_module.o \
        wf_module.o \
        symmetry_module.o \
        density_module.o\
        esm_rshell_module.o \
        ps_local_module.o ps_local_rs_module.o \
        esm_genpot_module.o \
        hartree_variables.o \
        hartree_sol_ffte_module.o \
        hartree_sol_module.o \
        hartree_mol_module.o \
        esm_hartree_module.o \
        hartree_module.o \
        ps_nloc_gth_module.o \
        ps_nloc_hgh_module.o \
        ps_pcc_module.o \
        ps_initrho_module.o \
        ps_nloc1_module.o \
        maskf_module.o \
        minimal_box_module.o \
        ps_nloc2_init_module.o \
        ps_nloc2_variables.o \
        ps_nloc_mr_module.o \
        ps_nloc3_module.o \
        ps_nloc2_module.o \
        construct_matrix_ps_nloc2_module.o \
        localpot2_variables.o \
        localpot2_density_module.o \
        localpot2_ion_module.o \
        localpot2_vh_module.o \
        localpot2_xc_module.o \
        localpot2_module.o \
        localpot2_Smatrix_module.o \
        localpot_module.o \
        nonlocal_module.o \
        ps_pcc_mol_module.o\
        ps_initrho_mol_module.o\
        ps_nloc2_mol_module.o\
        ps_local_mol_module.o\
        vloc_rho_mol_module.o \
        force_nloc2_mol_module.o \
        force_ion_mol_module.o \
        force_local_mol_module.o \
        eion_mol_module.o\
        ewald_variables.o \
        ewald_module.o\
        eion_module.o \
        force_ewald_module.o \
        gram_schmidt_u_module.o \
        gram_schmidt_m_module.o \
        gram_schmidt_t_module.o \
        gram_schmidt_module.o \
        vdw_grimme_module.o \
        xc_hybrid_module.o \
        fock_ffte_module.o \
        fock_fft_module.o \
        fock_cg_module.o \
        fock_module.o \
        xc_hf_module.o \
        xc_hse_module.o \
        xc_pw92_gth_module.o\
        xc_ldapz81_module.o \
        xc_ggapbe96_mol_module.o \
        xc_ggapbe96_module.o \
        xc_module.o \
        ps_pcc_force_module.o \
        hamiltonian_module.o\
        cgpc_gausseidel_module.o \
        cgpc_2_module.o \
        cgpc_diag_module.o \
        cgpc_seitsonen_module.o \
        cgpc_module.o \
        cg_u_module.o \
        cg_lobpcg_module.o \
        cg_module.o \
        fermi_mp_module.o\
        total_energy_module.o \
        mixing_broyden_module.o \
        mixing_module.o \
        esp_gather_module.o \
        io2_module.o \
        io_module.o \
        force_module.o \
        force_mol_module.o \
        esm_cylindrical_test.o \
        func2gp_module.o \
        localpot2_te_module.o \
        test_hpsi2_module.o \
        test_force_module.o \
        esp_calc_module.o \
        subspace_diag_variables.o \
        psv_initrho_module.o \
        random_initrho_module.o \
        hamiltonian_matrix_module.o \
        ChebyshevFilter_module.o \

#        momentum_mate_module.o \
#        esm_ion_module.o \
#        band_analysis_module.o \



DIAGLA_MOD = subspace_diag_la_module.o
DIAGSL_MOD = scalapack_module.o\
             subspace_mate_sl_module.o subspace_mate_sl_0_module.o\
             subspace_solv_sl_module.o \
             subspace_rotv_sl_module.o subspace_rotv_sl_0_module.o\
             subspace_diag_sl_module.o 
DIAG_MOD   = subspace_diag_module.o

MODS3 = timer_module.o iter_lin_solvers.o sseig.o prepare_sseig.o apply_sseig.o
MODS2 = momentum_module.o band_variables.o band_sseig_module.o band_module.o
MODS4 = scf_module.o scf_chefsi_module.o sweep_module.o atomopt_module.o \
        global_variables.o parameters_module.o

########################################################################
########################################################################

OBJ1 = realspace.o

DIR1 = ext1
MINPACOBJ  = $(DIR1)/ext_sub_minpac.o

DIR2 = ext2
EXTOBJ2 = $(DIR2)/polint.o \
          $(DIR2)/indexx.o \
          $(DIR2)/bberf.o \
          $(DIR2)/p4sn.o \
          $(DIR2)/ylm.o \
          $(DIR2)/dotp.o \
          $(DIR2)/spline.o\
          $(DIR2)/convert_capital.o\
          $(DIR2)/gaussj.o\

DIR3 = mdsource
MDOBJ = $(DIR3)/cpmdio2_module.o\
        $(DIR3)/active_band.o\
        $(DIR3)/asign_atom.o\
        $(DIR3)/init_ion.o\
        $(DIR3)/alloc_cpmd.o\
        $(DIR3)/blue_moon_module.o \
        $(DIR3)/bomd.o\
        $(DIR3)/getforce.o\
        $(DIR3)/getforce_cpmd.o\
        $(DIR3)/wf_force.o\
        $(DIR3)/calkin.o\
        $(DIR3)/calfke.o\
        $(DIR3)/calc_overlap_module.o \
        $(DIR3)/overlap_cpmd_module.o\
        $(DIR3)/rotorb_module.o\
        $(DIR3)/rotorb2.o\
        $(DIR3)/setv.o\
        $(DIR3)/vcom.o\
        $(DIR3)/mdio.o\
        $(DIR3)/mnhc.o\
        $(DIR3)/cpmd_variables.o

DIR4 = FFTE
FFTOBJ = $(DIR4)/pzfft3dv.o \
         $(DIR4)/fft235.o \
         $(DIR4)/kernel.o \
         $(DIR4)/zfft3d.o \
         $(DIR4)/pzfft3d.o

########################################################################
########################################################################

all : lda0 lda1 lda2 lda3 lda4 lda
	$(FC) $(LFLAGS) $(MODS1) $(DIAGLA_MOD) $(DIAGSL_MOD) $(DIAG_MOD) \
                        $(MODS2) $(MODS3) $(MODS4) $(OBJ1) \
	                $(EXTOBJ2) $(MINPACOBJ) \
                        $(MDOBJ) $(FFTOBJ) $(LAPACK_L)

lda0 : $(MODS1) $(DIAGLA_MOD) $(DIAGSL_MOD) $(DIAG_MOD) $(MODS3) $(MODS4) $(MODS2)

lda1 :
	cd $(DIR1) ; $(MAKE)
lda2 :
	cd $(DIR2) ; $(MAKE)
lda3 :
	cd $(DIR3) ; $(MAKE)
lda4 :
	cd $(DIR4) ; $(MAKE)

lda : $(OBJ1)

clean :
	rm -f *.o *.mod a.out mpif.h *.lst
	cd $(DIR1) ; $(MAKE) clean
	cd $(DIR2) ; $(MAKE) clean
	cd $(DIR3) ; $(MAKE) clean
	cd $(DIR4) ; $(MAKE) clean

