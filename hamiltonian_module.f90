MODULE hamiltonian_module

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
!$  use omp_lib
    implicit none
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

!!$  et0 = omp_get_wtime()

    call init_op_nonlocal

!$OMP parallel private( et0, et1 )

!$OMP workshare
    htpsi=(0.d0,0.d0)
!$OMP end workshare

! --- Kinetic energy ---

!$OMP barrier
    et0 = omp_get_wtime()

    call op_kinetic(k,tpsi,htpsi,n1,n2,ib1,ib2)

!$OMP barrier
    et1 = omp_get_wtime()

!$OMP single
    ett_hamil(1)=ett_hamil(1)+et1-et0
!$OMP end single

! --- local potential ---

!$OMP barrier
    et0 = omp_get_wtime()

    call op_localpot(s,n2-n1+1,ib2-ib1+1,tpsi,htpsi)

!$OMP barrier
    et1 = omp_get_wtime()

!$OMP single
    ett_hamil(2)=ett_hamil(2)+et1-et0
!$OMP end single

! --- nonlocal potential ---

!$OMP barrier
    et0 = omp_get_wtime()

    call op_nonlocal(k,tpsi,htpsi,n1,n2,ib1,ib2)

!$OMP barrier
    et1 = omp_get_wtime()

!$OMP single
    ett_hamil(3)=ett_hamil(3)+et1-et0
!$OMP end single

!$OMP end parallel

    call op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)

  END SUBROUTINE hamiltonian

END MODULE hamiltonian_module
