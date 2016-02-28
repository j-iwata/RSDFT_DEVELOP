MODULE nose_hoover_chain_module

  use parallel_module, only: myrank
  use calc_kine_temp_module, only: calc_kine, get_Ndof
  use cpmd_variables, only: inivel, linitnose

  implicit none

  PRIVATE
  PUBLIC :: nose_hoover_chain
  PUBLIC :: init_nose_hoover_chain
  PUBLIC :: noseene
  PUBLIC :: write_nose_data

  integer,parameter :: ncall=9
  real(8) :: dt_ST(ncall)
  integer,parameter :: nchain=4
  real(8) :: Qeta(nchain)
  real(8) :: eta(nchain)
  real(8) :: etadot(nchain)
  real(8) :: Feta(nchain)

  real(8),parameter :: TOJOUL=4.35975d-18 ! (J/hartree)
  real(8),parameter :: kB_J=1.380658d-23  ! (J/K)
  real(8),parameter :: kB=kB_J/TOJOUL     ! (hartree/K)

  real(8) :: kT,gkT

CONTAINS


  SUBROUTINE init_nose_hoover_chain( dt, temp, omega )

    implicit none
    real(8),intent(IN) :: dt, temp, omega
    real(8),parameter :: to_au=7.26d-7
    real(8) :: c,w,w2,sigma,pi,rnr(3)
    integer :: i, Ndof

    w  = omega * to_au
    w2 = w*w
!    pi = 0.0d0
    pi = acos(-1.0d0)

    call make_nose_time( dt )

    if ( .not.(inivel.or.linitnose) ) then
       call read_nose_data
       return
    end if

    call get_Ndof( Ndof )

    kT  = temp*kB
    gkT = Ndof*kT

    Qeta(1) = gkT/w2
    do i=2,nchain
       Qeta(i) = kT/w2
    end do

    do i=1,nchain
       eta(i) = 0.0d0
       sigma=sqrt( kT/Qeta(i) )
       call mprand(3,rnr)
       c=2.0d0*pi*rnr(1)
       etadot(i)=sqrt( log(xranf())*(-2.0d0) )*cos(c)*sigma
    end do

    return
  END SUBROUTINE init_nose_hoover_chain

  subroutine mprand(n,a)
    implicit none
    integer :: n
    real(8) :: a(n)
    integer i
    do i=1,n
       a(i)=xranf()
    end do
    return
  end subroutine mprand

  function xranf()
    implicit none
    real(8) :: xranf
    integer :: konst=125
    integer,save :: m=100001
    m=m*konst
    m=m-2796203*(m/2796203)
    xranf=dble(m)/2796203d0
    return
  end function xranf

  SUBROUTINE make_nose_time( dt )
    implicit none
    real(8),intent(IN) :: dt
    real(8) :: w_1,w_2,w_3,w_4,w_0
    w_1 = 0.192d0
    w_2 = 0.554910818409783619692725006662999d0
    w_3 = 0.124659619941888644216504240951585d0
    w_4 =-0.843182063596933505315033808282941d0
    w_0 = 1.0d0 - 2.0d0*( w_1 + w_2 + w_3 + w_4 )
    dt_ST(1) = w_4*dt
    dt_ST(2) = w_3*dt
    dt_ST(3) = w_2*dt
    dt_ST(4) = w_1*dt
    dt_ST(5) = w_0*dt
    dt_ST(6) = w_1*dt
    dt_ST(7) = w_2*dt
    dt_ST(8) = w_3*dt
    dt_ST(9) = w_4*dt
    return
  END SUBROUTINE make_nose_time

!----------------------------------------------------------


  SUBROUTINE nose_hoover_chain( Velocity )
    implicit none
    real(8),intent(INOUT) :: Velocity(:,:)
    integer :: i
    do i=1,ncall
       call nose( dt_ST(i), Velocity )
    end do
  END SUBROUTINE nose_hoover_chain


  SUBROUTINE nose( step, Velocity )

    implicit none
    real(8),intent(IN) :: step
    real(8),intent(INOUT) :: Velocity(:,:)
    real(8) :: ekinp,aaex
    integer :: i

    call calc_kine( Velocity, ekinp )

    Feta(1) = 2.0d0*ekinp - gkT
    do i=2,nchain
       ekinp = 0.5d0*Qeta(i-1)*etadot(i-1)**2
       Feta(i) = 2.0d0*ekinp - kT
    end do

    etadot(nchain) = etadot(nchain) &
         + 0.25d0*step*Feta(nchain)/Qeta(nchain)

    do i=1,nchain-1
       aaex = exp( -0.125d0*step*etadot(nchain+1-i) )
       etadot(nchain-i) = etadot(nchain-i)*aaex*aaex &
            + 0.25d0*step*Feta(nchain-i)*aaex/Qeta(nchain-i)
    end do

    aaex = exp( -0.5d0*step*etadot(1) )

    Velocity(:,:) = Velocity(:,:)*aaex

    do i=1,nchain
       eta(i) = eta(i) + 0.5d0*step*etadot(i)
    end do

    call calc_kine( Velocity, ekinp )

    Feta(1) = 2.0d0*ekinp - gkT

    do i=1,nchain-1
       aaex = exp( -0.125d0*step*etadot(i+1) )

       etadot(i) = etadot(i)*aaex*aaex &
            + 0.25d0*step*Feta(i)*aaex/Qeta(i)

       ekinp = 0.5d0*Qeta(i)*etadot(i)**2
       Feta(i+1) = 2.0d0*ekinp - kT
    end do

    etadot(nchain) = etadot(nchain) &
         + 0.25d0*step*Feta(nchain)/Qeta(nchain)

    return
  END SUBROUTINE nose

!---------------------------------------------------------------------------
!--- Bath energy -----------------------------------------------------------

  SUBROUTINE noseene( bathe )
    implicit none
    real(8),intent(OUT) :: bathe
    integer :: i
    bathe = 0.5d0*Qeta(1)*etadot(1)**2 + gkT*eta(1)
    do i=2,nchain
       bathe = bathe + 0.5d0*Qeta(i)*etadot(i)**2 + kT*eta(i)
    end do
    return
  END SUBROUTINE noseene

!-------------------------------------------------------------------------
!----------------Restart Nose'--------------------------------------------
!-------------------------------------------------------------------------

  SUBROUTINE write_nose_data
    implicit none
    integer :: i
    if ( myrank == 0 ) then
       open(999,file="Bathdata.dat1")
       write(999,*)(eta(i),etadot(i),i=1,nchain) 
       write(999,*) kT,gkT
       write(999,*)(Qeta(i),i=1,nchain) 
       write(999,*)(Feta(i),i=1,nchain)
       close(999)
    end if
    return
  END SUBROUTINE write_nose_data


  SUBROUTINE read_nose_data
    implicit none
    integer :: i,ierr
    include 'mpif.h'
    if ( myrank == 0 ) then
       open(999,file="Bathdata.dat")
       read(999,*)(eta(i),etadot(i),i=1,nchain) 
       read(999,*) kT,gkT
       read(999,*)(Qeta(i),i=1,nchain)
       read(999,*)(Feta(i),i=1,nchain)
       close(999)
    end if
    call mpi_bcast(eta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(etadot,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(gkT,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Qeta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Feta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    return
  END SUBROUTINE read_nose_data


END MODULE nose_hoover_chain_module
