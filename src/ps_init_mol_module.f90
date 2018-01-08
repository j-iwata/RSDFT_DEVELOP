MODULE ps_init_mol_module

  use ps_local_mol_module, only: init_ps_local_mol, construct_ps_local_mol
  use ps_pcc_mol_module, only: init_ps_pcc_mol, construct_ps_pcc_mol
  use ps_initrho_mol_module, only: init_ps_initrho_mol,construct_ps_initrho_mol
  use ps_nloc2_init_module, only: ps_nloc2_init
  use ps_nloc2_mol_module, only: prep_ps_nloc2_mol

  implicit none

  PRIVATE
  PUBLIC ps_init_mol

CONTAINS

  SUBROUTINE ps_init_mol( Gcut, rho )
    implicit none
    real(8),intent(IN) :: Gcut
    real(8),intent(INOUT) :: rho(:,:)

    call init_ps_local_mol(Gcut)
    call init_ps_pcc_mol
    call init_ps_initrho_mol

    !call Construct_RgridMol(Igrid)

    call construct_ps_local_mol
    call construct_ps_pcc_mol
    call construct_ps_initrho_mol( rho )
    !call normalize_density

    call ps_nloc2_init(Gcut)
    call prep_ps_nloc2_mol

    !call ConstructBoundary_RgridMol(Md,Igrid)

  END SUBROUTINE ps_init_mol

END MODULE ps_init_mol_module
