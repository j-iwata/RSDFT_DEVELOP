MODULE hamiltonian_module
use parallel_module, only: myrank

  use kinetic_module
  use localpot_module
  use nonlocal_module
  use fock_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: hamiltonian,op_kinetic,op_localpot,op_nonlocal &
           ,ctt_hamil,ett_hamil

  real(8) :: ctt_hamil(4),ett_hamil(4)

CONTAINS

  SUBROUTINE hamiltonian(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    integer :: i,ib,i1,i2,i3,j,lma,m,ML0,n,nb
    integer :: a1,a2,a3,b1,b2,b3,ierr,nreq
    real(8) :: ct0,ct1,et0,et1

    htpsi=(0.d0,0.d0)

! --- Kinetic energy ---

    call watch(ct0,et0)

    call op_kinetic(k,tpsi,htpsi,n1,n2,ib1,ib2)

    call watch(ct1,et1)
    ctt_hamil(1)=ctt_hamil(1)+ct1-ct0 ; ett_hamil(1)=ett_hamil(1)+et1-et0

! --- local potential ---

    call op_localpot(s,n2-n1+1,ib2-ib1+1,tpsi,htpsi)
do ib=ib1,ib2
do i=n1,n2
write(150,'(2I5,4g20.7)') ib,i,tpsi(i,ib),htpsi(i,ib)
enddo
enddo
stop

    call watch(ct0,et0)
    ctt_hamil(2)=ctt_hamil(2)+ct0-ct1 ; ett_hamil(2)=ett_hamil(2)+et0-et1

! --- nonlocal potential ---
#ifdef _SHOWALL_OP_
write(200+myrank,*) '----------------------------hamiltonian'
write(200+myrank,'(A17,6I7)') 'k,s,n1,n2,ib1,ib2',k,s,n1,n2,ib1,ib2
write(200+myrank,*) '----------------------------hamiltonian'
#endif
    call op_nonlocal(k,s,tpsi,htpsi,n1,n2,ib1,ib2,.false.)
!write(170,*) 'ib,i,tpsi(i,ib),htpsi(i,ib)'
do ib=ib1,ib2
do i=n1,n2
!write(170,'(2I5,4g20.7)') ib,i,tpsi(i,ib),htpsi(i,ib)
enddo
enddo
!stop

    call watch(ct1,et1)
    ctt_hamil(3)=ctt_hamil(3)+ct1-ct0 ; ett_hamil(3)=ett_hamil(3)+et1-et0

    call op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)

    call watch(ct0,et0)
    ctt_hamil(4)=ctt_hamil(4)+ct0-ct1 ; ett_hamil(4)=ett_hamil(4)+et0-et1

  END SUBROUTINE hamiltonian

END MODULE hamiltonian_module
