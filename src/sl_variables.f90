module sl_variables

  implicit none

  private
  public :: slinfo2
  public :: sl0
  public :: sl1
  public :: sl2
  public :: Hsub

  type slinfo2
    integer :: nband
    integer :: myrow,  mycol
    integer :: nprow,  npcol
    integer :: mbsize, nbsize
    integer,allocatable :: map_1to2(:,:)
    integer,allocatable :: map_2to1(:,:)
    integer,allocatable :: usermap(:,:)
    integer :: distribution_method
    integer :: icontxt_sys
    integer :: icontxt_a, icontxt_b
    integer :: desca(9), descb(9), descz(9)
    character(8) :: idiag
    character(1) :: uplo
    integer :: lwork, lrwork, liwork
    integer :: ldr, ldc
  end type slinfo2

  type(slinfo2) :: sl0
  type(slinfo2) :: sl1
  type(slinfo2) :: sl2

  complex(8),allocatable :: Hsub(:,:)

end module sl_variables
