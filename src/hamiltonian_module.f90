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
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:,ib1:)
    real(8),intent(OUT) :: htpsi(n1:,ib1:)
#else
    complex(8),intent(IN)  :: tpsi(n1:,ib1:)
    complex(8),intent(OUT) :: htpsi(n1:,ib1:)
#endif
    real(8) :: ttmp(2)

!$OMP parallel

!$OMP workshare
    htpsi=(0.0d0,0.0d0)
!$OMP end workshare

    call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic( tpsi, htpsi, k )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    call op_localpot( tpsi, htpsi, s )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal( tpsi, htpsi, k, s, ib1, ib2 )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,3) )

!$OMP end parallel

    call watchb( ttmp )

    call op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)

    call watchb( ttmp, time_hmlt(1,4) )

  END SUBROUTINE hamiltonian


  SUBROUTINE hamiltonian_test(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:,ib1:)
    real(8),intent(OUT) :: htpsi(n1:,ib1:)
#else
    complex(8),intent(IN)  :: tpsi(n1:,ib1:)
    complex(8),intent(OUT) :: htpsi(n1:,ib1:)
#endif
    real(8) :: ttmp(2)
    integer :: ib

!$OMP parallel

!$OMP workshare
    htpsi=(0.0d0,0.0d0)
!$OMP end workshare

    call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic( tpsi(:,ib1:ib2), htpsi(:,ib1:ib2), k )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    call op_localpot( tpsi, htpsi, s )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal( tpsi(:,ib1:ib2), htpsi(:,ib1:ib2), k, s )

!$OMP barrier

    call watchb_omp( ttmp, time_hmlt(1,3) )

!$OMP end parallel

    call watchb( ttmp )

    call op_fock(k,s,n1,n2,ib1,ib2,tpsi,htpsi)

    call watchb( ttmp, time_hmlt(1,4) )

  END SUBROUTINE hamiltonian_test


END MODULE hamiltonian_module
