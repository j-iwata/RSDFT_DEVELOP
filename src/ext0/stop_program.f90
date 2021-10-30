subroutine stop_program( indx )
  implicit none
  character(*),intent(in) :: indx
  integer :: ierr, myrank
  include 'mpif.h'
  call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
  if ( myrank == 0 ) then
    write(*,'(/,a25,"  stop_program is called !!!  ",a25)') repeat("-",25),repeat("-",25)
    write(*,*) indx
  end if
  call MPI_Abort( MPI_COMM_WORLD, 7, ierr )
  ! call MPI_Finalize( ierr )
  ! stop "stop@stop_program"
end subroutine stop_program


subroutine stop_program_f( indx )
  implicit none
  character(*),intent(in) :: indx
  integer :: ierr, myrank
  include 'mpif.h'
  call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
  write(*,'(/,a25," stop_program_f is called !!! ",a25)') repeat("-",25),repeat("-",25)
  write(*,*) indx, myrank
  call MPI_Abort( MPI_COMM_WORLD, 7, ierr )
end subroutine stop_program_f
