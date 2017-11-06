MODULE hartree_ene_module

  use rgrid_module, only: dV
  use parallel_module, only: comm_grid
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_ene

  include 'mpif.h'

CONTAINS

  SUBROUTINE calc_hartree_ene( rho, Vh, Eh )
    implicit none
    real(8),intent(IN)  :: rho(:), Vh(:)
    real(8),intent(OUT) :: Eh
    integer :: i
    Eh = 0.5d0*dV*sum( rho*Vh )
    call rsdft_allreduce_sum( Eh, comm_grid )
  END SUBROUTINE calc_hartree_ene

END MODULE hartree_ene_module
