MODULE hamiltonian_ncol_module

  use ps_spinorbit_module, only: op_ps_spinorbit
  use noncollinear_module, only: op_xc_noncollinear

  implicit none

  PRIVATE
  PUBLIC :: hamiltonian_ncol

CONTAINS

  SUBROUTINE hamiltonian_ncol( k, n1,n2, tpsi, hpsi )
    implicit none
    integer,intent(IN) :: k,n1,n2
    complex(8),intent(IN) :: tpsi(:,:)
    complex(8),intent(INOUT) :: hpsi(:,:)
    call op_xc_noncollinear( tpsi, hpsi )
    call op_ps_spinorbit( k, n1,n2, tpsi, hpsi )
  END SUBROUTINE hamiltonian_ncol

END MODULE hamiltonian_ncol_module
