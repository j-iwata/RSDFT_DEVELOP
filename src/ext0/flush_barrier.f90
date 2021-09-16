subroutine flush_barrier( u )
  implicit none
  integer,intent(in) :: u
  integer :: ierr
  include 'mpif.h'
  call flush(u)
  call MPI_Barrier( MPI_COMM_WORLD, ierr )
end subroutine flush_barrier
