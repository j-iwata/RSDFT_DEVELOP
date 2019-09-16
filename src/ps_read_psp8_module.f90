module ps_read_psp8_module

  use var_ps_member, only: ps1d, ps_allocate_ps1d

  implicit none

  private
  public :: ps_read_psp8

contains

  subroutine ps_read_psp8( g, psp )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: psp
    integer :: lloc,lmaxabinit,i,j,l,j0,j1
    integer,allocatable :: nproj(:)
    real(8) :: psdat,pspcod,pspxc,r2well,zatm,zion
    real(8) :: rchrg,fchrg,qchrg,dr,pi4,r,tmp
    character(1) :: dummy_text

    write(*,'(a40," ps_read_psp8(start)")') repeat("-",40)

    read(g,*) dummy_text
    read(g,*) zatm, zion, psdat
    read(g,*) pspcod,pspxc,lmaxabinit,lloc,psp%mr,r2well
    read(g,*) rchrg,fchrg,qchrg

    allocate( nproj(0:lmaxabinit) ); nproj=0
    read(g,*) nproj(:)

    psp%norb = sum( nproj )
    call ps_allocate_ps1d( psp )

    psp%Zelement = nint(zatm)
    psp%Zps = nint(zion)

    read(g,*)

    do l=0,lmaxabinit
       j0=sum(nproj(0:l))-nproj(l)+1
       j1=j0+nproj(l)-1
       read(g,*) dummy_text, (psp%anorm(j),j=j0,j1)
       do i=1,psp%mr
          read(g,*) dummy_text, psp%rad(i), (psp%viod(i,j),j=j0,j1)
       end do
    end do

    read(g,*)

    do i=1,psp%mr
       read(g,*) dummy_text, psp%rad(i), psp%vql(i)
    end do

    dr=psp%rad(2)-psp%rad(1)
    pi4=4.0d0*acos(-1.0d0)

    if ( fchrg > 0.0d0 ) then
       do i=1,psp%mr
          read(g,*) dummy_text,r,psp%cdc(i)
       end do
       psp%cdc=psp%cdc/pi4
    end if

    do i=1,psp%mr
       read(g,*) dummy_text,r,psp%cdd(i)
       psp%cdd(i)=psp%cdd(i)*psp%rad(i)*psp%rad(i)
    end do

    psp%rab = dr

    i=0
    do l=0,lmaxabinit
    do j=1,nproj(l)
       i=i+1
       psp%lo(i)=l
       psp%no(i)=j
    end do
    end do

    do j=1,psp%norb
       do i=psp%Mr,1,-1
          if ( abs(psp%viod(i,j)) >=1.d-13 ) then
             psp%Rps(j) = psp%rad(i+1)
             psp%NRps(j)= i+1
             exit
          end if
       end do
    end do

    if ( any(psp%anorm/=0.0d0) ) then
       do j=1,psp%norb
          psp%inorm(j)=1
          if ( psp%anorm(j) < 0.0d0 ) psp%inorm(j)=-1
          tmp = sqrt( abs( psp%anorm(j) ) )
          psp%viod(:,j) = psp%viod(:,j)*tmp
       end do
    end if

    write(*,*) "*** ABINIT psp8 Format (psp8) ***"
    write(*,*) "Znuc=",psp%Zps
    write(*,*) "# of radial mesh points =",psp%Mr
    write(*,*) "# of orbitals =",psp%norb
    if ( psp%norb > 0 ) then
    write(*,*) "angular momentum =",psp%lo(1:psp%norb)
    write(*,*) "cut off radius =",psp%Rps(1:psp%norb)
    write(*,*) "# of grid points within cut off radius",psp%NRps(1:psp%norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( psp%inorm(i)*psp%anorm(i),i=1,psp%norb )
    end if
    write(*,*) "sum(rhov)=",sum(psp%cdd*psp%rab)
    write(*,*) "sum(rhoc)=",sum(psp%cdc*psp%rab*psp%rad**2)*pi4

    write(*,'(a40," ps_read_psp8(end)")') repeat("-",40)

    return
  end subroutine ps_read_psp8

end module ps_read_psp8_module
