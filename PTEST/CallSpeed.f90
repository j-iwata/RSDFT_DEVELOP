PROGRAM main
  implicit none
  include 'mpif.h'
  integer :: count
  integer :: i
  integer :: n
  integer :: ierr
  integer :: myrank,nprocs
  logical :: disp_switch
  real(8) :: ct0,ct1,ct_call,ct_if
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  disp_switch=(myrank==0)
  read(1,*) n
  
  call cpu_time(ct0)
  count=0
  do i=1,n
    call call_time(count,disp_switch)
  enddo
  call cpu_time(ct1)
  ct_call=ct1-ct0
  
  write(*,*) repeat('-',40)
  call cpu_time(ct0)
  count=0
  if (myrank==0) then
    do i=1,n
      count=count+1
      if (disp_switch) write(2,'(I7,A30)') count,'if'
    enddo
  endif
  call cpu_time(ct1)
  ct_if=ct1-ct0

  write(*,'(a6,g20.7)') 'call=',ct_call
  write(*,'(a6,g20.7)') 'if  =',ct_if

  call MPI_FINALIZE(ierr)
CONTAINS
  SUBROUTINE call_time(count_,disp_switch_)
    implicit none
    integer,intent(INOUT) :: count_
    logical,intent(IN) :: disp_switch_
    count_=count_+1
    if (disp_switch_) write(3,'(I7,A30)') count_,'inside subroutine'
  END SUBROUTINE call_time
END PROGRAM main
