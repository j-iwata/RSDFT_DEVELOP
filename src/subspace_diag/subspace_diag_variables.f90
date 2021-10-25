module subspace_diag_variables

  implicit none

  private
  public :: algo_sd
  public :: MB_diag, Hsub, Vsub, mat_block &
           ,NBLK1, NBLK2, zero, one, TYPE_MAIN

  include 'mpif.h'

  integer,allocatable :: mat_block(:,:)

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

  integer :: ialgo_sd = 1

contains

  integer function algo_sd()
    use io_tools_module, only: IOTools_readIntegerKeyword
    implicit none
    call IOTools_readIntegerKeyword( 'IALGO_SD', ialgo_sd )
    algo_sd = ialgo_sd
  end function algo_sd

end module subspace_diag_variables
