MODULE nose_hoover_chain_ele_module

  use parallel_module, only: myrank
  use cpmd_variables, only: mstocck, inivel, linitnosee

  implicit none

  PRIVATE
  PUBLIC :: nose_hoover_chain_ele
  PUBLIC :: init_nose_hoover_chain_ele
  PUBLIC :: write_nosee_data

  integer,parameter :: ncall=9
  real(8) :: dt_ST(ncall)
  integer,parameter :: nchain=4
  real(8) :: Qeta(nchain)
  real(8) :: eta(nchain)
  real(8) :: etadot(nchain)
  real(8) :: Feta(nchain+1)

  integer :: Ndof_e
  real(8) :: ekinw
  real(8) :: beta_e

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


  SUBROUTINE init_nose_hoover_chain_ele( dt, ekinw_in, omega, Ndof_e_in, ene )

    implicit none
    real(8),intent(IN) :: dt, ekinw_in, omega
    integer,intent(IN) :: Ndof_e_in
    real(8),optional,intent(OUT) :: ene
    real(8),parameter :: to_au=7.26d-7
    real(8) :: w,w2
    integer :: i

    call make_nosee_time( dt )

    if ( .not.(inivel.or.linitnosee) ) then
       call read_nosee_data
       return
    end if

    w  = omega * to_au
    w2 = w*w

!    Ndof_e = Ndof_e_in*3.0d0
    Ndof_e = Ndof_e_in

    ekinw = ekinw_in

    beta_e = dble(Ndof_e)/(2.0d0*ekinw)

    Qeta(1)        = 2.0d0*ekinw/w2
    Qeta(2:nchain) = 1.0d0/(beta_e*w2)

    etadot(1) = sqrt( 2.0d0*ekinw/Qeta(1) )
    do i=2,nchain
       etadot(i) = sqrt( (1.0d0/beta_e)/Qeta(i) )
    end do

    if ( present(ene) ) call noseeneele( ene )

    return
  END SUBROUTINE init_nose_hoover_chain_ele


  SUBROUTINE nose_hoover_chain_ele( fke, psi_v, n1,n2, ene )
    implicit none
    real(8),intent(INOUT) :: fke
    real(8),intent(INOUT) :: psi_v(:,:,:,:)
    integer,intent(IN) :: n1,n2
    real(8),optional,intent(OUT) :: ene
    integer :: i
    real(8) :: scale
    scale=1.0d0
    do i=1,ncall
       call enosmove( fke, dt_ST(i), scale )
    end do
    psi_v(:,n1:n2,:,:) = scale*psi_v(:,n1:n2,:,:)
    if ( present(ene) ) call noseeneele( ene )
  END SUBROUTINE nose_hoover_chain_ele


  SUBROUTINE enosmove( fke, step, sctot )

    implicit none
    real(8),intent(IN)    :: step
    real(8),intent(INOUT) :: fke,sctot
    real(8) :: aae,f1,f2
    integer :: i

    Feta(1) = 2.0d0*( fke - ekinw )
    do i=2,nchain+1
       Feta(i) = Qeta(i-1)*etadot(i-1)**2 - 1.0d0/beta_e
    end do

    aae = exp( -0.125d0*step*etadot(nchain-1) )
    etadot(nchain) = etadot(nchain)*aae*aae &
                   + 0.25d0*step*Feta(nchain)*aae/Qeta(nchain)
    Feta(nchain+1) = Qeta(nchain)*etadot(nchain)**2 -1.0d0/beta_e
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
    fke = fke*aae*aae

    do i=1,nchain
       eta(i) = eta(i) + 0.5d0*step*etadot(i)
    end do

    Feta(:) = 0.0d0
    Feta(1) = 2.0d0*( fke - ekinw )
    Feta(nchain-1) = Qeta(nchain)*etadot(nchain)**2 - 1.0d0/beta_e

    do i=1,nchain-1
       aae = exp(-0.125d0*step*etadot(i+1))
       etadot(i) = etadot(i)*aae*aae + 0.25d0*step*Feta(i)*aae/Qeta(i)
       Feta(i+1) = Feta(i+1) + Qeta(i)*etadot(i)**2 - 1.0d0/beta_e
    end do
    aae = exp(-0.125d0*step*etadot(nchain-1))
    etadot(nchain) = etadot(nchain)*aae*aae &
                   + 0.25d0*step*Feta(nchain)*aae/Qeta(nchain)

    return
  END SUBROUTINE enosmove


  SUBROUTINE noseeneele( enose )
    implicit none
    real(8),intent(OUT) :: enose
    integer :: i
    enose = 0.5d0*Qeta(1)*etadot(1)**2 + 2.0d0*ekinw*eta(1)
    do i=2,nchain
       enose = enose + 0.5d0*Qeta(i)*etadot(i)**2 + eta(i)/beta_e
    end do
    enose = enose + eta(nchain-1)/beta_e
  END SUBROUTINE noseeneele


  SUBROUTINE write_nosee_data

    implicit none
    integer :: i

    if ( myrank == 0 ) then
       open(999,file="Bathdata_ele.dat1")
       write(999,*)(eta(i),etadot(i),i=1,nchain) 
       write(999,*)(Feta(i),i=1,nchain+1) 
       write(999,*)(Qeta(i),i=1,nchain) 
       write(999,*) Ndof_e, ekinw
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
       read(999,*) Ndof_e, ekinw
       close(999)
    end if

    call mpi_bcast(eta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(etadot,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Feta,nchain+1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Qeta,nchain,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(Ndof_e,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ekinw,1,mpi_real8,0,mpi_comm_world,ierr)

    beta_e = dble(Ndof_e)/(2.0d0*ekinw)

    return

  END SUBROUTINE read_nosee_data


END MODULE nose_hoover_chain_ele_module
