MODULE nonlocal_module

  use pseudopot_module, only: pselect
  use ps_nloc1_module
  use ps_nloc2_module
! use ps_nloc3_module
  use ps_nloc_mr_module

  implicit none

  PRIVATE
  PUBLIC :: op_nonlocal

CONTAINS

  SUBROUTINE op_nonlocal(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    select case( pselect )
    case(2,4)
       call op_ps_nloc2(k,tpsi,htpsi,n1,n2,ib1,ib2)
    case(5)
       call op_ps_nloc_mr(k,tpsi,htpsi,n1,n2,ib1,ib2)
    case default
       stop "pselect/=2,4,5 are not implemented"
    end select

  END SUBROUTINE op_nonlocal

END MODULE nonlocal_module
