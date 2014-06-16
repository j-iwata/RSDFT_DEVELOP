MODULE rgrid_module

  implicit none

  PRIVATE
  PUBLIC :: Ngrid,Hgrid,dV,zdV,Igrid &
           ,read_rgrid,read_oldformat_rgrid,construct_rgrid &
           ,ParallelInit_RgridSol

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


  SUBROUTINE read_oldformat_rgrid(rank,unit)
    integer,intent(IN) :: rank,unit
    Ngrid(:)=0
    if ( rank == 0 ) then
       read(unit,*) Ngrid(1:3)
       write(*,*) "Ngrid(1:3)=",Ngrid(1:3)
    end if
    call send_rgrid(0)
  END SUBROUTINE read_oldformat_rgrid


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
  END SUBROUTINE construct_rgrid


  SUBROUTINE ParallelInit_RgridSol( np, np_grid, pinfo_grid )
    implicit none
    integer,intent(IN)  :: np(3), np_grid
    integer,intent(OUT) :: pinfo_grid(8,0:np_grid-1)
    integer :: i,j,n,i1,i2,i3
    integer,allocatable :: ntmp(:,:)
    n=maxval(np)
    allocate( ntmp(0:n-1,3) ) ; ntmp=0
    do j=1,3
       do i=0,Ngrid(j)-1
          n=mod(i,np(j))
          ntmp(n,j)=ntmp(n,j)+1
       end do
    end do
    n=-1
    i= 0
    do i3=0,np(3)-1
    do i2=0,np(2)-1
    do i1=0,np(1)-1
       n=n+1
       pinfo_grid(1,n)=sum( ntmp(0:i1,1) )-ntmp(i1,1)
       pinfo_grid(2,n)=ntmp(i1,1)
       pinfo_grid(3,n)=sum( ntmp(0:i2,2) )-ntmp(i2,2)
       pinfo_grid(4,n)=ntmp(i2,2)
       pinfo_grid(5,n)=sum( ntmp(0:i3,3) )-ntmp(i3,3)
       pinfo_grid(6,n)=ntmp(i3,3)
       pinfo_grid(7,n)=i
       pinfo_grid(8,n)=ntmp(i1,1)*ntmp(i2,2)*ntmp(i3,3)
       i=i+pinfo_grid(8,n)
    end do
    end do
    end do
    deallocate( ntmp )
  END SUBROUTINE ParallelInit_RgridSol


END MODULE rgrid_module
