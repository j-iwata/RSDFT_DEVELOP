MODULE TestAssert
  implicit none
  PRIVATE
  PUBLIC :: assert
  PUBLIC :: assert_init
  PUBLIC :: assert_finalize
  integer,parameter :: unit=7
  character(8),parameter :: filename='TEST.log'
  integer :: count, count_T, count_F
CONTAINS
  SUBROUTINE assert_init
    implicit none
    count   = 0
    count_T = 0
    count_F = 0
  END SUBROUTINE assert_init

  SUBROUTINE assert( bool, msg )
    implicit none
    logical,intent(IN) :: bool
    character(*),intent(IN) :: msg
    if ( bool ) then
      count_T = count_T + 1
      write(*,100) '[True]',msg
    elseif ( .not. bool ) then
      count_F = count_F + 1
      write(*,100) '[False]',msg
    endif
100 FORMAT(1X,A8,1X,A)
  END SUBROUTINE assert

  SUBROUTINE assert_finalize( msg )
    implicit none
    character(*),intent(IN),optional :: msg
    count = count_T + count_F
    if ( present(msg) ) write(*,101) msg
    write(*,102) count_T, count
101 FORMAT(1X,A)
102 FORMAT(1X,'RESULT: ',I4,' / ',I4,' : PASSED')
  END SUBROUTINE assert_finalize
END MODULE TestAssert
