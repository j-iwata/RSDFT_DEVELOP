MODULE rgrid_module

  implicit none

  PRIVATE
  PUBLIC :: Ngrid,Hgrid,dV,zdV,Igrid &
           ,read_rgrid,construct_rgrid &
           ,parallel_rgrid

  integer :: Ngrid(0:3),Igrid(2,0:3)
  real(8) :: Hgrid(3)
  real(8) :: dV
#ifdef _DRSDFT_
  real(8) :: zdV
#else
  complex(8) :: zdV
#endif

CONTAINS

  SUBROUTINE read_rgrid(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(5) :: cbuf,ckey
    Ngrid(:)=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "NGRID" ) then
             backspace(unit)
             read(unit,*) cbuf,Ngrid(1:3)
          end if
       end do
999    continue
       write(*,*) "Ngrid(1:3)=",Ngrid(1:3)
    end if
    call send_rgrid(0)
  END SUBROUTINE read_rgrid

!  SUBROUTINE read_rgrid(unit)
!    integer,intent(IN) :: unit
!    Ngrid(:)=0
!    read(unit,*) Ngrid(1:3)
!    write(*,*) "Ngrid(1:3)=",Ngrid(1:3)
!  END SUBROUTINE read_rgrid


  SUBROUTINE send_rgrid(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Ngrid(1),3,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_rgrid


  SUBROUTINE construct_rgrid(aa)
    implicit none
    real(8),intent(IN) :: aa(3,3)
    real(8) :: Vaa
    Vaa = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
         +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
         -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
    Ngrid(0)=Ngrid(1)*Ngrid(2)*Ngrid(3)
    Hgrid(1)=sqrt( sum(aa(1:3,1)**2) )/Ngrid(1)
    Hgrid(2)=sqrt( sum(aa(1:3,2)**2) )/Ngrid(2)
    Hgrid(3)=sqrt( sum(aa(1:3,3)**2) )/Ngrid(3)
    dV=abs(Vaa)/Ngrid(0)
    zdV=dV
    Igrid(:,:)=0
    Igrid(2,1)=Ngrid(1)-1
    Igrid(2,2)=Ngrid(2)-1
    Igrid(2,3)=Ngrid(3)-1
    Igrid(1,0)=1
    Igrid(2,0)=Ngrid(0)
  END SUBROUTINE construct_rgrid


  SUBROUTINE parallel_rgrid(np_grid,myrank)
    implicit none
    integer,intent(IN) :: np_grid(3),myrank
    integer :: i,j,n,i1,i2,i3
    integer,allocatable :: np1(:),np2(:),np3(:)
    allocate( np1(np_grid(1)) ) ; np1=0
    allocate( np2(np_grid(2)) ) ; np2=0
    allocate( np3(np_grid(3)) ) ; np3=0
    do i=1,Ngrid(1)
       n=mod(i-1,np_grid(1))+1
       np1(n)=np1(n)+1
    end do
    do i=1,Ngrid(2)
       n=mod(i-1,np_grid(2))+1
       np2(n)=np2(n)+1
    end do
    do i=1,Ngrid(3)
       n=mod(i-1,np_grid(3))+1
       np3(n)=np3(n)+1
    end do
    n=-1
    i=0
    do i3=1,np_grid(3)
    do i2=1,np_grid(2)
    do i1=1,np_grid(1)
       n=n+1
       if ( n == myrank ) then
          Igrid(1,0)=i+1
          Igrid(2,0)=i+np1(i1)*np2(i2)*np3(i3)
          Igrid(1,1)=sum(np1(1:i1))-np1(i1)
          Igrid(2,1)=sum(np1(1:i1))-1
          Igrid(1,2)=sum(np2(1:i2))-np2(i2)
          Igrid(2,2)=sum(np2(1:i2))-1
          Igrid(1,3)=sum(np3(1:i3))-np3(i3)
          Igrid(2,3)=sum(np3(1:i3))-1
       end if
       i=i+np1(i1)*np2(i2)*np3(i3)
    end do
    end do
    end do
    deallocate( np3,np2,np1 )
  END SUBROUTINE parallel_rgrid

END MODULE rgrid_module
