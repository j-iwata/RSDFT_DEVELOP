MODULE nose_hoover_chain_ele_module

  implicit none

  PRIVATE
  PUBLIC :: nose_hoover_chain_ele
  PUBLIC :: nosepae
  PUBLIC :: noseeneele
  PUBLIC :: write_nosee_data
  PUBLIC :: read_nosee_data
  PUBLIC :: make_nosee_time

  integer,parameter :: ncall=9
  real(8) :: dtsuz(ncall)

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
    dtsuz(1) = w_4*dt
    dtsuz(2) = w_3*dt
    dtsuz(3) = w_2*dt
    dtsuz(4) = w_1*dt
    dtsuz(5) = w_0*dt
    dtsuz(6) = w_1*dt
    dtsuz(7) = w_2*dt
    dtsuz(8) = w_3*dt
    dtsuz(9) = w_4*dt
    return
  END SUBROUTINE make_nosee_time


  SUBROUTINE nosepae( dt )

    use cpmd_variables, only: etadot, nche, ekinw, qnosee, nedof, mstocck, wnosee, wnose0

    implicit none
    real(8),intent(IN) :: dt
    real(8),parameter :: wntau=7.26d-7
    real(8) :: wnosee2
    real(8) :: w_yosh7_1,w_yosh7_2,w_yosh7_3,w_yosh7_4,w_yosh7_0
    integer :: i, n
    integer,parameter :: ip=1
    integer,parameter :: nedof0=6

    wnosee = wnose0 * wntau
    n = mstocck(1,1)

    nedof = int(n*nedof0)

    wnosee2 = wnosee * wnosee
    qnosee(1)      = 2.0d0*ekinw/wnosee2
    qnosee(2:nche) = 2.0d0*ekinw/wnosee2/dble(nedof)

    etadot(:,:)=0.0d0

    etadot(1,ip) = sqrt( 2.0d0*ekinw/qnosee(1) )
    do i=2,nche
       etadot(i,ip) = sqrt( 2.0d0*ekinw/qnosee(i)/dble(nedof) )
    end do

    call make_nosee_time( dt )

    return
  END SUBROUTINE nosepae


  SUBROUTINE nose_hoover_chain_ele( fke, psi_v, n1,n2 )
    implicit none
    real(8),intent(INOUT) :: fke
    real(8),intent(INOUT) :: psi_v(:,:,:,:)
    integer,intent(IN) :: n1,n2
    integer :: i
    real(8) :: scale
    scale=1.0d0
    do i=1,ncall
       call enosmove( fke, dtsuz(i), scale )
    end do
    psi_v(:,n1:n2,:,:) = scale*psi_v(:,n1:n2,:,:)
  END SUBROUTINE nose_hoover_chain_ele


  SUBROUTINE enosmove( ekinc, step, sctot )

    use cpmd_variables

    implicit none
    real(8),intent(IN)    :: step
    real(8),intent(INOUT) :: ekinc,sctot
    real(8) ::  ckewant,ckine,aae,f1,f2
    integer ::  l
    integer,parameter :: ip=1

    ckewant = ekinw/dble(nedof)
    feta(1) = 2.d0*(ekinc - ekinw)
    do l=2,nche+1
       ckine = 0.5d0*qnosee(l-1)*etadot(l-1,ip)*etadot(l-1,ip)
       feta(l) = 2.d0*(ckine - ckewant)
    end do

    aae = exp( -0.125d0*step*etadot(nche-1,ip) )
    etadot(nche,ip) = etadot(nche,ip)*aae*aae &
         + 0.25d0*step*feta(nche)*aae/qnosee(nche)
    ckine = 0.5d0*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
    feta(nche+1) = 2.d0*(ckine - ckewant)
    f1 = feta(nche-1)
    f2 = feta(nche+1)
    feta(nche-1) = f1 + f2

    do l=1,nche-1
       aae = exp(-0.125d0*step*etadot(nche+1-l,ip))
       etadot(nche-l,ip) = etadot(nche-l,ip)*aae*aae &
            + 0.25d0*step*feta(nche-l)*aae/qnosee(nche-l)
    end do

    aae = exp(-0.5d0*step*etadot(1,ip))
    sctot = sctot*aae
    ekinc = ekinc*aae*aae

    do l=1,nche
       etap(l,ip) = etap(l,ip) + 0.5d0*step*etadot(l,ip)
       feta(l) = 0.d0
    end do

    feta(1) = 2.d0*(ekinc - ekinw)
    ckine = 0.5d0*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
    feta(nche-1) = 2.d0*(ckine - ckewant)

    do l=1,nche-1
       aae = exp(-0.125d0*step*etadot(l+1,ip))
       etadot(l,ip) = etadot(l,ip)*aae*aae &
            + 0.25d0*step*feta(l)*aae/qnosee(l)
       ckine = 0.5d0*qnosee(l)*etadot(l,ip)*etadot(l,ip)
       feta(l+1) = feta(l+1) + 2.d0*(ckine - ckewant)
    end do
    aae = exp(-0.125d0*step*etadot(nche-1,ip))
    etadot(nche,ip) = etadot(nche,ip)*aae*aae &
         + 0.25d0*step*feta(nche)*aae/qnosee(nche)

    return
  END SUBROUTINE enosmove


  SUBROUTINE noseeneele( enose )

    use cpmd_variables

    implicit none
    real(8) :: enose
    integer :: i
    integer,parameter :: ip=1

    enose=0.0d0
    enose = enose + 0.5d0*qnosee(1)*etadot(1,ip)*etadot(1,ip) &
         + 2.0d0*ekinw*etap(1,ip)
    do i=2,nche
       enose = enose + 0.5d0*qnosee(i)*etadot(i,ip)*etadot(i,ip) &
            + 2.0d0*ekinw*etap(i,ip)/dble(nedof)
    end do
    enose = enose + 2.d0*ekinw*etap(nche-1,ip)/dble(nedof)
    return
  END SUBROUTINE noseeneele


  SUBROUTINE write_nosee_data

    use cpmd_variables

    implicit none
    integer :: i
    integer,parameter :: ip=1

    if ( myrank == 0 ) then
       open(999,file="Bathdata_ele.dat1")
       write(999,*)(etap(i,ip),etadot(i,ip),i=1,nche) 
       write(999,*)(feta(i),i=1,nche+1) 
       write(999,*)(qnosee(i),i=1,nche) 
       write(999,*) nedof
       close(999)
    end if
    return

  END SUBROUTINE write_nosee_data


  SUBROUTINE read_nosee_data

    use cpmd_variables

    implicit none
    integer :: i, ierr
    integer,parameter :: ip=1

    if ( myrank == 0 ) then
       open(999,file="Bathdata_ele.dat")
       read(999,*)(etap(i,ip),etadot(i,ip),i=1,nche) 
       read(999,*)(feta(i),i=1,nche+1)
       read(999,*)(qnosee(i),i=1,nche) 
       read(999,*) nedof
       close(999)
    end if

    call mpi_bcast(etap(1,ip),nche,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(etadot(1,ip),nche,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(feta(1),nche+1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(qnosee(1),nche,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nedof,1,mpi_integer,0,mpi_comm_world,ierr)

    return

  END SUBROUTINE read_nosee_data


END MODULE nose_hoover_chain_ele_module
