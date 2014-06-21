MODULE gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_t_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt

  integer :: iswitch_algorithm = 0 

CONTAINS

  SUBROUTINE gram_schmidt(n0,n1,k,s)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    select case( iswitch_algorithm )
    case default
       call gram_schmidt_t(n0,n1,k,s)
    case( 1 )
       call gram_schmidt_m(n0,n1,k,s)
    end select
  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module
