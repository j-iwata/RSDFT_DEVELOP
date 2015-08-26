MODULE esm_rshell_module

  use hsort_module

  implicit none

  PRIVATE
  PUBLIC :: construct_lattice_rshell,Mshell,num_member_shell,lattice_shell &
            ,read_rshell

  integer :: Mshell
  integer,allocatable :: num_member_shell(:)
  integer,allocatable :: lattice_shell(:,:,:)

  integer :: Mcell
  integer :: Ncell1,Ncell2,Ncell3
  
CONTAINS

  SUBROUTINE read_rshell(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(5) :: cbuf,ckey
    include 'mpif.h'
    Ncell1=0
    Ncell2=0
    Ncell3=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "NCELL" ) then
             backspace(unit)
             read(unit,*) cbuf,Ncell1,Ncell2,Ncell3
          end if
       end do
999    continue
       write(*,*) "Ncell1,Ncell2,Ncell3=",Ncell1,Ncell2,Ncell3
    end if
    call mpi_bcast(Ncell1,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Ncell2,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Ncell3,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE read_rshell

  SUBROUTINE construct_lattice_rshell
    use aa_module
    implicit none
    integer,allocatable :: indx(:),lat_tmp(:,:)
    integer :: i,m,mm,i1,i2,i3
    real(8),parameter :: eps=1.d-5
    real(8),allocatable :: rraa(:)
    real(8) :: x,y,z,r

    Mcell  = (2*Ncell1+1)*(2*Ncell2+1)*(2*Ncell3+1)

    allocate( rraa(Mcell)      ) ; rraa=0.d0
    allocate( indx(Mcell)      ) ; indx=0
    allocate( lat_tmp(3,Mcell) ) ; lat_tmp=0

    i=0
    do i3=-Ncell3,Ncell3
    do i2=-Ncell2,Ncell2
    do i1=-Ncell1,Ncell1
       i=i+1
       x=aa(1,1)*i1+aa(1,2)*i2+aa(1,3)*i3
       y=aa(2,1)*i1+aa(2,2)*i2+aa(2,3)*i3
       z=aa(3,1)*i1+aa(3,2)*i2+aa(3,3)*i3
       rraa(i)=x*x+y*y+z*z
       lat_tmp(1,i)=i1
       lat_tmp(2,i)=i2
       lat_tmp(3,i)=i3
    end do
    end do
    end do
    call indexx(Mcell,rraa,indx)

    r=rraa(indx(1))
    Mshell=1
    m=1
    mm=m
    do i=2,Mcell
       if ( abs(rraa(indx(i))-r) > eps ) then
          Mshell=Mshell+1
          m=1
          r=rraa(indx(i))
       else
          m=m+1
       end if
       mm=max(m,mm)
    end do

    allocate( lattice_shell(3,mm,Mshell) ) ; lattice_shell=0
    allocate( num_member_shell(Mshell)   ) ; num_member_shell=1

    r=rraa(indx(1))
    Mshell=1
    m=1
    do i=2,Mcell
       if ( abs(rraa(indx(i))-r) > eps ) then
          Mshell=Mshell+1
          m=1
          r=rraa(indx(i))
          lattice_shell(1:3,m,Mshell)=lat_tmp(1:3,indx(i))
       else
          m=m+1
          num_member_shell(Mshell)=m
          lattice_shell(1:3,m,Mshell)=lat_tmp(1:3,indx(i))
       end if
    end do

    deallocate( lat_tmp )
    deallocate( indx )
    deallocate( rraa )

!    mm=0
!    do i=1,Mshell
!       do m=1,num_member_shell(i)
!          mm=mm+1
!          write(*,'(1x,9i6)') mm,i,m,lattice_shell(1:3,m,i),sum(lattice_shell(1:3,m,i)**2)
!       end do
!    end do

  END SUBROUTINE construct_lattice_rshell

END MODULE esm_rshell_module
