MODULE global_variables

  use aa_module
  use bb_module
  use atom_module
  use rgrid_module
  use ggrid_module
  use bz_module
  use fd_module
  use strfac_module
  use pseudopot_module
  use ps_local_module
  use ps_pcc_module
  use ps_initrho_module
  use ps_nloc1_module
  use ps_nloc2_init_module
  use ps_nloc2_module
  use bc_module
  use electron_module
  use density_module
  use parallel_module
  use wf_module
  use localpot_module
  use nonlocal_module
  use ewald_module
  use gram_schmidt_module
  use gram_schmidt_t_module
  use hartree_module
  use xc_module
  use scalapack_module
  use subspace_diag_module
  use subspace_diag_la_module
  use subspace_diag_sl_module
  use subspace_mate_sl_module
  use subspace_solv_sl_module
  use subspace_rotv_sl_module
  use kinetic_module
  use hamiltonian_module
  use cgpc_module
  use cg_module
  use total_energy_module
  use mixing_module
  use esp_gather_module
  use fermi_module
  use watch_module
  use io_module
  use array_bound_module
  use atomopt_module

  implicit none

  integer :: Diter, Nsweep, Ndiag
  integer :: iswitch_scf,iswitch_opt,iswitch_band
  real(8) :: etime_limit
  logical :: disp_switch

END MODULE global_variables
