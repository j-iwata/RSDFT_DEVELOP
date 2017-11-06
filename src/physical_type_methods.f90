MODULE physical_type_methods

  use parallel_module
  use rgrid_variables, only: dV
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: dSpatialIntegral

CONTAINS


  SUBROUTINE dSpatialIntegral( d )
    implicit none
    real(8) :: d(:)
    call rsdft_allreduce_sum( d, comm_grid )
    d=d*dV
  END SUBROUTINE dSpatialIntegral


END MODULE physical_type_methods
