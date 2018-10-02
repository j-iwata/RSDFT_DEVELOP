MODULE subspace_diag_variables

  implicit none

  PRIVATE
  PUBLIC :: MB_diag, Hsub, Vsub, mat_block ,nLB&
           ,NBLK1, NBLK2, zero, one, TYPE_MAIN
!
!  PUBLIC :: esp
  PUBLIC :: Hsub_e, Vsub_e
!

  include 'mpif.h'

  integer,allocatable :: mat_block(:,:)
  integer,allocatable :: nLB(:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN = MPI_REAL8
  real(8),allocatable :: Hsub(:,:), Vsub(:,:)
  real(8),parameter :: zero=0.d0,one=1.d0
#else
#ifdef _NO_MPI_COMPLEX16_
  integer,parameter :: TYPE_MAIN = MPI_DOUBLE_COMPLEX
#else
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
#endif
  complex(8),allocatable :: Hsub(:,:), Vsub(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
#endif
  integer :: MB_diag,NBLK1,NBLK2

! Solver
!  real(8),allocatable ::  esp(:,:,:)
! Eigen
    real(8), allocatable :: Hsub_e(:,:)
    real(8), allocatable :: Vsub_e(:,:)

END MODULE subspace_diag_variables
