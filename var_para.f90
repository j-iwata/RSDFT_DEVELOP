MODULE VarPara
  implicit none
  
  include 'mpif.h'

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN = MPI_REAL8
#else
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
#endif

!  integer :: node_partition(1:6)
!  integer :: nprocs,nprocs_g
!  integer :: myrank,myrank_g

END MODULE VarPara
