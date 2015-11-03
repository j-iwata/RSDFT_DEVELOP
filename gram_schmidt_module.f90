MODULE gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_t_module
 !use gram_schmidt_u_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt, read_gram_schmidt

  integer :: iswitch_algorithm = 0
  include 'mpif.h'

CONTAINS

  SUBROUTINE read_gram_schmidt(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(2) :: cbuf,ckey
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "GS" ) then
             backspace(unit)
             read(unit,*) cbuf, iswitch_algorithm
             exit
          end if
       end do
999    continue
       write(*,*) "GS:iswitch_algorithm=",iswitch_algorithm
    end if
    call mpi_bcast(iswitch_algorithm,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE read_gram_schmidt

  SUBROUTINE gram_schmidt(n0,n1,k,s)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    select case( iswitch_algorithm )
    case default
       call gram_schmidt_t(n0,n1,k,s)
    case( 1 )
       call gram_schmidt_m(n0,n1,k,s)
    !case( 2 )
    !   call gram_schmidt_u(n0,n1,k,s)
    end select
  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module
