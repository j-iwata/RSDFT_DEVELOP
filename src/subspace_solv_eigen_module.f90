MODULE subspace_solv_eigen_module

#ifdef _EIGEN_
  use eigen_libs
#endif

  implicit none

  PRIVATE
  PUBLIC :: subspace_solv_eigen

CONTAINS

  SUBROUTINE subspace_solv_eigen( N, a, w, z )
    implicit none
    integer,intent(IN) :: N
    real(8),intent(INOUT) :: a(:,:)
    real(8),intent(OUT) :: w(:),z(:,:)
#ifdef _EIGEN_
    call eigen_sx( N,N,a,size(a,1),w,z,size(z,1) )
#else
    w=0.0d0
    z=0.0d0
#endif
  END SUBROUTINE subspace_solv_eigen

END MODULE subspace_solv_eigen_module
