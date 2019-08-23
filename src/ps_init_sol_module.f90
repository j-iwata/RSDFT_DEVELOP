MODULE ps_init_sol_module

  use ps_local_module, only: init_ps_local, construct_ps_local
  use ps_pcc_module, only: init_ps_pcc, construct_ps_pcc
  use strfac_module, only: construct_strfac, destruct_strfac
  use ps_initrho_module, only: construct_ps_initrho
  use ps_nloc_initiate_module, only: ps_nloc_initiate
  use allel_module, only: flag_allel, init_ae_local_allel

  implicit none

  PRIVATE
  PUBLIC :: ps_init_sol

CONTAINS

  SUBROUTINE ps_init_sol( Gcut, rho )
    implicit none
    real(8),intent(IN) :: Gcut
    real(8),intent(INOUT) :: rho(:,:)

    if ( flag_allel ) then
       call init_ae_local_allel
    else
       call init_ps_local
       call init_ps_pcc
    end if

    call construct_strfac  !----- structure factor

    call construct_ps_local
    call construct_ps_pcc
    call construct_ps_initrho( rho )

    call destruct_strfac   !----- structure factor

    if ( flag_allel ) then
       return
    else
       call ps_nloc_initiate( Gcut )
    end if

  END SUBROUTINE ps_init_sol

END MODULE ps_init_sol_module
