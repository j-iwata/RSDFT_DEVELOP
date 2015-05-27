MODULE VarPara
  implicit none
  PRIVATE
  include 'mpif.h'

#ifdef _DRSDFT_
  integer,parameter,PUBLIC :: TYPE_MAIN = MPI_REAL8
#else
  integer,parameter,PUBLIC :: TYPE_MAIN = MPI_COMPLEX16
#endif

END MODULE VarPara
