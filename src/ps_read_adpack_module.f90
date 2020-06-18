module ps_read_adpack_module

  use var_ps_member, only: ps1d, ps_allocate_ps1d

  implicit none

  private
  public :: ps_read_adpack

contains

  subroutine ps_read_adpack( g, psp )
    implicit none
    integer,intent(in) :: g
    type(ps1d),intent(inout) :: psp
    real(8) :: zion, zatm, rdummy, vdummy, r, tmp, dx
    real(8),allocatable :: x(:)
    integer :: icount,i,idummy
    integer :: iorb,norb,l,nproj(0:9)
    character(50) :: cbuf

    write(*,'(a40," ps_read_adpack(start)")') repeat("-",40)

    rewind g

    icount=-1
    do
       read(g,*) cbuf
       if ( cbuf == "<project.energies" ) then
          read(g,*) psp%norb
          cycle
       end if
       if ( cbuf == "<Pseudo.Potentials" ) then
          icount=0
          cycle
       end if
       if ( icount >= 0 ) then
          if ( cbuf == "Pseudo.Potentials>" ) then
             exit
          else
             icount=icount+1
          end if
       end if
    end do
    psp%Mr=icount

    call ps_allocate_ps1d( psp )

    rewind g
    do
       read(g,*) cbuf
       if ( cbuf == "valence.electron" ) then
          backspace(g)
          read(g,*) cbuf, zion
          psp%Zps=nint(zion)
       end if
       if ( cbuf == "AtomSpecies" ) then
          backspace(g)
          read(g,*) cbuf, zatm
          psp%Zelement=nint(zatm)
       end if
       if ( cbuf == "<project.energies" ) then
          read(g,*) idummy
          nproj=0
          do iorb=1,psp%norb
             read(g,*) l, psp%anorm(iorb), vdummy
             psp%lo(iorb)=l
             nproj(l)=nproj(l)+1
             psp%no(iorb)=nproj(l)
          end do
          exit
       end if
    end do

    allocate( x(psp%Mr) ); x=0.0d0

    icount=-1
    do
       read(g,*) cbuf
       if ( cbuf == "<Pseudo.Potentials" ) then
          icount=0
          cycle
       end if
       if ( icount >= 0 ) then
          if ( cbuf == "Pseudo.Potentials>" ) then
             exit
          else
             icount=icount+1
             backspace(g)
             read(g,*) x(icount), psp%rad(icount), psp%vql(icount), &
                  ( psp%viod(icount,iorb), vdummy, iorb=1,psp%norb )
             do iorb=1,psp%norb
                psp%viod(icount,iorb)=psp%viod(icount,iorb)*psp%rad(icount)
             end do
          end if
       end if
    end do

!    do i=1,size(psp%viod,1)
!      write(10,'(1x,10f20.15)') psp%rad(i),(psp%viod(i,icount),icount=1,psp%norb)
!    end do

    do i=1,psp%Mr-1
       dx = x(i+1) - x(i)
       psp%rab(i) = psp%rad(i)*dx
    end do
    psp%rab(psp%Mr) = psp%rad(psp%Mr)*dx

    deallocate( x )

    icount=-1
    do
       read(g,*,END=9) cbuf
       if ( cbuf == "<density.PCC" ) then
          icount=0
          cycle
       end if
       if ( icount >= 0 ) then
          if ( cbuf == "density.PCC>" ) then
             exit
          else
             icount=icount+1
             backspace(g)
             read(g,*) rdummy, r, psp%cdc(icount)
          end if
       end if
    end do

9   continue

    do iorb=1,psp%norb
       do i=psp%Mr,1,-1
          if ( abs(psp%viod(i,iorb))**2 >= 1.d-8 ) then
             icount=min(i+1,psp%Mr)
             psp%Rps(iorb) = psp%rad(icount)
             psp%NRps(iorb)= icount
             exit
          end if
       end do
    end do

    if ( any(psp%anorm/=0.0d0) ) then
       do iorb=1,psp%norb
          psp%inorm(iorb)=1
          if ( psp%anorm(iorb) < 0.0d0 ) psp%inorm(iorb)=-1
          tmp = sqrt( abs( psp%anorm(iorb) ) )
          psp%viod(:,iorb) = psp%viod(:,iorb)*tmp
       end do
    end if

    write(*,*) "*** OpenMX ADPACK Format (.vps) ***"
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
    write(*,*) "sum(rhoc)=",sum(psp%cdc*psp%rab*psp%rad**2)*4.0d0*acos(-1.0d0)

    write(*,'(a40," ps_read_adpack(end)")') repeat("-",40)

    return
  end subroutine ps_read_adpack

end module ps_read_adpack_module
