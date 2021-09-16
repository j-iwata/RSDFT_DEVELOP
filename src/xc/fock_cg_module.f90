module fock_cg_module

  implicit none

  private
  public :: Fock_CG_Double
  public :: Fock_CG_DoubleComplex

contains

  subroutine Fock_CG_Double( trho, tVh )
    use hartree_mol_module, only: calc_hartree_mol
    implicit none
    real(8),intent(in) :: trho(:)
    real(8),intent(inout) :: tVh(:)
    real(8) :: Edummy
    integer :: n1,n2
    n1 = 1
    n2 = size(trho)
    call calc_hartree_mol( n1, n2, 1, trho, tVh, Edummy, tol=1.d-8 )
    return
  end subroutine Fock_CG_Double

  subroutine Fock_CG_DoubleComplex( trho, tVh )
    implicit none
    complex(8),intent(in) :: trho(:)
    complex(8),intent(inout) :: tVh(:)
    call stop_program( 'Fock_CG_DoubleComplex is not implemented...' )
    return
  end subroutine Fock_CG_DoubleComplex

end module fock_cg_module
