MODULE timer_module

  implicit none

  PRIVATE
  PUBLIC :: timer,tic,toc,clear_timer

  include 'mpif.h'

  type timer
     real(8) :: ct
     real(8) :: et
     real(8) :: tmpct
     real(8) :: tmpet
  end type timer

CONTAINS

  subroutine tic(arg,comm)
    type(timer), intent(out) :: arg
    integer, intent(in) :: comm
    integer :: ierr
    call mpi_barrier(comm,ierr)
    call cpu_time(arg%tmpct)
    arg%tmpet = mpi_wtime()
  end subroutine tic

  subroutine toc(arg,comm)
    type(timer), intent(inout) :: arg
    integer, intent(in) :: comm
    real(8) :: ct, et
    integer :: ierr
    call mpi_barrier(comm,ierr)
    call cpu_time(ct)
    et = mpi_wtime()
    arg%ct = arg%ct + ct - arg%tmpct
    arg%et = arg%et + et - arg%tmpet
  end subroutine toc

  subroutine clear_timer(arg)
    type(timer), intent(out) :: arg
    arg%ct = 0d0
    arg%et = 0d0
    arg%tmpct = 0d0
    arg%tmpet = 0d0
  end subroutine clear_timer

END MODULE timer_module
