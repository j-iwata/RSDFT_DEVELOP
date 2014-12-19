!--------------------------------------------
! TM format
! J.-L. Martins' code ( latest version is atom-5.7.0 )
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
MODULE ps_read_TM_module

  implicit none

  PRIVATE
  PUBLIC :: ps_read_TM, ps_tm

  integer,parameter :: nrmax = 5000
  integer,parameter :: lmax  = 3

  TYPE tm
     integer :: nrr
     integer :: norb
     integer :: NRc(lmax+1)
     integer :: inorm(lmax+1)
     integer :: lo(lmax+1)
     integer :: no(lmax+1)
     real(8) :: znuc
     real(8) :: vps(nrmax,lmax+1)
     real(8) :: ups(nrmax,lmax+1)
     real(8) :: vql(nrmax)
     real(8) :: cdc(nrmax)
     real(8) :: cdd(nrmax)
     real(8) :: rr(nrmax)
     real(8) :: rx(nrmax)
     real(8) :: Rc(lmax+1)
     real(8) :: anorm(lmax+1)
     real(8) :: Dij(lmax+1,lmax+1)
  END TYPE tm

  type(tm) :: ps_tm

CONTAINS

  SUBROUTINE ps_read_TM(g)
    implicit none
    integer,intent(IN) :: g
    character(2)  :: icorr,nameat
    character(3)  :: irel
    character(4)  :: nicore
    character(10) :: iray(6),ititle(7)
    integer :: i,j,norb,numu,nrr0,nrr,Lref
    integer,allocatable :: tlo(:),tinorm(:)
    real(8) :: a0,b0,viou(lmax,nrmax),pi4
    real(8),parameter :: eps=1.d-7
    real(8),allocatable :: tviod(:,:),tanorm(:)

    write(*,'(a40," ps_read_tm")') repeat("-",40)

    ps_tm%nrr     =0
    ps_tm%norb    =0
    ps_tm%znuc    =0.0d0
    ps_tm%vps(:,:)=0.0d0
    ps_tm%ups(:,:)=0.0d0
    ps_tm%vql(:)  =0.0d0
    ps_tm%cdc(:)  =0.0d0
    ps_tm%cdd(:)  =0.0d0
    ps_tm%rr(:)   =0.0d0
    ps_tm%rx(:)   =0.0d0
    ps_tm%Rc(:)   =0.0d0
    ps_tm%NRc(:)  =0
    ps_tm%anorm(:)=0.0d0
    ps_tm%inorm(:)=0
    ps_tm%lo(:)   =0
    ps_tm%no(:)   =0
    ps_tm%Dij(:,:)=0.0d0

    pi4 = 4.0d0*acos(-1.0d0)

    allocate( tlo(2*lmax) )
    allocate( tviod(lmax,nrmax) )
    allocate( tanorm(lmax) )
    allocate( tinorm(lmax) )

    read(g) nameat,icorr,irel,nicore,iray,ititle,norb,numu,nrr,a0,b0,ps_tm%znuc
    nrr=nrr+1
    read(g) ( ps_tm%rr(i) ,i=2,nrr )
    do i=1,norb
       read(g) tlo(i),(tviod(i,j),j=2,nrr)
    end do
    do i=1,numu
       read(g) tlo(i+norb),(viou(tlo(i+norb)+1,j),j=2,nrr)
    end do
    read(g) ( ps_tm%cdc(i),i=2,nrr)
    read(g) ( ps_tm%cdd(i),i=2,nrr)
    read(g) ( ps_tm%vql(i),i=2,nrr)
    read(g) norb
    do i=1,norb
       read(g) tinorm(tlo(i)+1),tanorm(tlo(i)+1)
    end do

    j=0
    do i=1,norb
       if ( tinorm(tlo(i)+1) /= 0 ) then
          j=j+1
          ps_tm%lo(j)=tlo(i)
          ps_tm%inorm(j)=tinorm(tlo(i)+1)
          ps_tm%anorm(j)=tanorm(tlo(i)+1)
          ps_tm%vps(2:nrr,j)=tviod(i,2:nrr)
       else
          Lref=tlo(i)
       end if
    end do
    norb=j

    deallocate( tinorm, tanorm, tviod, tlo )

! Values at the origin
!
    ps_tm%rr(1)    = 0.0d0
    ps_tm%vql(1)   = ps_tm%vql(2)
    ps_tm%vps(1,:) = ps_tm%vps(2,:)

    do j=1,norb
       ps_tm%vps(:,j) = ps_tm%vps(:,j)*ps_tm%rr(:)
    end do
    do i=2,nrr
       ps_tm%cdc(i) = ps_tm%cdc(i)*0.5d0/(pi4*ps_tm%rr(i)**2)
    end do

    ps_tm%cdc(1) = ps_tm%cdc(2)
    ps_tm%cdd(1) = ps_tm%cdd(2)

! dr/dx
!
    ps_tm%rx(:) = b0*( ps_tm%rr(:) + a0 )

! Convert unit from Rydberg to Hartree
!
    ps_tm%vps(:,:) = ps_tm%vps(:,:)*sqrt(0.5d0)
    ps_tm%anorm(:) = ps_tm%anorm(:)*sqrt(0.5d0)
    ps_tm%vql(:)   = ps_tm%vql(:)*0.5d0

! Cut-off radius
!
    do j=1,norb
       nrr0=1
       do i=nrr,1,-1
          if ( abs( ps_tm%vps(i,j) ) > 1.d-2 ) then
             nrr0=i
             exit
          end if
       end do ! i
       do i=nrr0,nrr
          if ( abs( ps_tm%vps(i,j) ) < eps ) then
             ps_tm%Rc(j)  = ps_tm%rr(i)
             ps_tm%NRc(j) = i
             exit
          end if
       end do ! i
    end do ! j

    ps_tm%nrr  = nrr
    ps_tm%norb = norb

    write(*,*) "*** TM Format ***"
    write(*,*) "Znuc=",ps_tm%znuc
    write(*,*) "# of radial mesh points =",nrr
    write(*,*) "# of orbitals =",norb
    write(*,*) "angular momentum =",ps_tm%lo(1:norb)
    write(*,*) "cut off radius =",ps_tm%Rc(1:norb)
    write(*,*) "# of grid points within cut off radius",ps_tm%NRc(1:norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( ps_tm%inorm(i)*ps_tm%anorm(i),i=1,norb )
    write(*,*) "Dij ="
    write(*,'(1x,9f10.5)') (( ps_tm%Dij(i,j),i=1,norb ),j=1,norb)
    write(*,*) "sum(rhov)=",sum(ps_tm%cdd*ps_tm%rx)
    write(*,*) "sum(rhoc)=",sum(ps_tm%cdc*ps_tm%rx*(ps_tm%rr)**2)*pi4

    write(*,'(a40," ps_read_tm(end)")') repeat("-",40)

  END SUBROUTINE ps_read_TM

END MODULE ps_read_TM_module
