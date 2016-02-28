MODULE nose_hoover_chain_ele_module

  use parallel_module, only: myrank
  use cpmd_variables, only: mstocck, inivel, linitnosee

  implicit none

  PRIVATE
  PUBLIC :: nose_hoover_chain_ele
  PUBLIC :: init_nose_hoover_chain_ele
  PUBLIC :: noseeneele
  PUBLIC :: write_nosee_data

  integer,parameter :: ncall=9
  real(8) :: dt_ST(ncall)
  integer,parameter :: nchain=4
  real(8) :: Qeta(nchain)
  real(8) :: eta(nchain)
  real(8) :: etadot(nchain)
  real(8) :: Feta(nchain)

  integer :: nedof
  real(8) :: ekinw

CONTAINS


  SUBROUTINE make_nosee_time( dt )
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
  END SUBROUTINE make_nosee_time


  SUBROUTINE init_nose_hoover_chain_ele( dt, ekinw_in, omega )

    implicit none
    real(8),intent(IN) :: dt, ekinw_in, omega
    real(8),parameter :: to_au=7.26d-7
    real(8) :: w,w2
    integer :: i,n
    integer,parameter :: nedof0=6

    call make_nosee_time( dt )

    if ( .not.(inivel.or.linitnosee) ) then
       call read_nosee_data
       return
    end if

    w  = omega * to_au
    w2 = w*w

    n = mstocck(1,1)
    nedof = int(n*nedof0)

    ekinw = ekinw_in

    Qeta(1)        = 2.0d0*ekinw/w2
    Qeta(2:nchain) = 2.0d0*ekinw/dble(nedof)/w2

    etadot(1) = sqrt( 2.0d0*ekinw/Qeta(1) )
    do i=2,nchain
       etadot(i) = sqrt( 2.0d0*ekinw/dble(nedof)/Qeta(i) )
    end do

    return
  END SUBROUTINE init_nose_hoover_chain_ele


  SUBROUTINE nose_hoover_chain_ele( fke, psi_v, n1,n2 )
    implicit none
    real(8),intent(INOUT) :: fke
    real(8),intent(INOUT) :: psi_v(:,:,:,:)
    integer,intent(IN) :: n1,n2
    integer :: i
    real(8) :: scale
    scale=1.0d0
    do i=1,ncall
       call enosmove( fke, dt_ST(i), scale )
    end do
    psi_v(:,n1:n2,:,:) = scale*psi_v(:,n1:n2,:,:)
  END SUBROUTINE nose_hoover_chain_ele


  SUBROUTINE enosmove( ekinc, step, sctot )

    implicit none
    real(8),intent(IN)    :: step
    real(8),intent(INOUT) :: ekinc,sctot
    real(8) :: ckewant,ckine,aae,f1,f2
    integer :: i

    ckewant = ekinw/dble(nedof)
    Feta(1) = 2.d0*(ekinc - ekinw)
    do i=2,nchain+1
       ckine = 0.5d0*Qeta(i-1)*etadot(i-1)*etadot(i-1)
       Feta(i) = 2.d0*(ckine - ckewant)
    end do

    aae = exp( -0.125d0*step*etadot(nchain-1) )
    etadot(nchain) = etadot(nchain)*aae*aae &
         + 0.25d0*step*Feta(nchain)*aae/Qeta(nchain)
    ckine = 0.5d0*Qeta(nchain)*etadot(nchain)*etadot(nchain)
    Feta(nchain+1) = 2.d0*(ckine - ckewant)
    f1 = Feta(nchain-1)
    f2 = Feta(nchain+1)
    Feta(nchain-1) = f1 + f2

    do i=1,nchain-1
       aae = exp(-0.125d0*step*etadot(nchain+1-i))
       etadot(nchain-i) = etadot(nchain-i)*aae*aae &
            + 0.25d0*step*Feta(nchain-i)*aae/Qeta(nchain-i)
    end do

    aae = exp(-0.5d0*step*etadot(1))
    sctot = sctot*aae
    ekinc = ekinc*aae*aae

    do i=1,nchain
       eta(i) = eta(i) + 0.5d0*step*etadot(i)
       Feta(i) = 0.0d0
    end do

    Feta(1) = 2.d0*(ekinc - ekinw)
    ckine = 0.5d0*Qeta(nchain)*etadot(nchain)*etadot(nchain)
    Feta(nchain-1) = 2.d0*(ckine - ckewant)

    do i=1,nchain-1
       aae = exp(-0.125d0*step*etadot(i+1))
       etadot(i) = etadot(i)*aae*aae &
            + 0.25d0*step*Feta(i)*aae/Qeta(i)
       ckine = 0.5d0*Qeta(i)*etadot(i)*etadot(i)
       Feta(i+1) = Feta(i+1) + 2.d0*(ckine - ckewant)
    end do
    aae = exp(-0.125d0*step*etadot(nchain-1))
    etadot(nchain) = etadot(nchain)*aae*aae &
         + 0.25d0*step*Feta(nchain)*aae/Qeta(nchain)

    return
  END SUBROUTINE enosmove


  SUBROUTINE noseeneele( enose )

    implicit none
    real(8) :: enose
    integer :: i

    enose = 0.5d0*Qeta(1)*etadot(1)**2 + 2.0d0*ekinw*eta(1)
    do i=2,nchain
       enose = enose + 0.5d0*Qeta(i)*etadot(i)**2 &
            + 2.0d0*ekinw/dble(nedof)*eta(i)
    end do
    enose = enose + 2.d0*ekinw/dble(nedof)*eta(nchain-1)
    return
  END SUBROUTINE noseeneele


  SUBROUTINE write_nosee_data

    implicit none
    integer :: i

    if ( myrank == 0 ) then
       open(999,file="Bathdata_ele.dat1")
       write(999,*)(eta(i),etadot(i),i=1,nchain) 
       write(999,*)(Feta(i),i=1,nchain+1) 
       write(999,*)(Qeta(i),i=1,nchain) 
       write(999,*) nedof
       close(999)
    end if
    return

  END SUBROUTINE write_nosee_data


  SUBROUTINE read_nosee_data

    implicit none
    integer :: i, ierr
    include 'mpif.h'

    if ( myrank == 0 ) then
       open(999,file="Bathdata_ele.dat")
       read(999,*)(eta(i),etadot(i),i=1,nchain) 
       read(999,*)(Feta(i),i=1,nchain+1)
       read(999,*)(Qeta(i),i=1,nchain) 
       read(999,*) nedof
       close(999)
    end if

    call mpi_bcast(eta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(etadot,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Feta,nchain+1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Qeta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nedof,1,mpi_integer,0,mpi_comm_world,ierr)

    return

  END SUBROUTINE read_nosee_data


END MODULE nose_hoover_chain_ele_module
