MODULE nonlocal_module

  use pseudopot_module, only: pselect
  use ps_nloc1_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module

  use PSnonLocOpG2

  use parallel_module, only: myrank

  implicit none

  PRIVATE
  PUBLIC :: op_nonlocal

CONTAINS

  SUBROUTINE op_nonlocal(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
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
    case(3)
       call op_ps_nloc3(k,tpsi,htpsi,n1,n2,ib1,ib2)
    case(5)
       call op_ps_nloc_mr(k,tpsi,htpsi,n1,n2,ib1,ib2)
    case(102)
#ifdef _SHOWALL_OP_
write(200+myrank,*) '----------------------------op_nonlocal'
write(200+myrank,'(A17,6I7)') 'k,s,n1,n2,ib1,ib2',k,s,n1,n2,ib1,ib2
write(200+myrank,*) '----------------------------op_nonlocal'
#endif
      call op_ps_nloc2_uspp(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    case default
       stop "pselect/=2,4,5,102 are not implemented"
    end select

  END SUBROUTINE op_nonlocal

END MODULE nonlocal_module
