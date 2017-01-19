MODULE io2_module

  use rgrid_module, only: Igrid
  use wf_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: read_io2
  PUBLIC :: read_data_io2
  PUBLIC :: read_data2_io2

  integer :: myrank
  integer :: nprocs_new, nmax
  integer,allocatable :: nmap(:),iomap(:,:,:)
  integer :: MB_0_IO, MB_1_IO

  character(16) :: map_file_name = "input_wfiodir"

CONTAINS


  SUBROUTINE read_io2( SYStype, info )
    implicit none
    integer,intent(IN)  :: SYStype
    integer,intent(OUT) :: info
    integer,parameter  :: unit0=10
    integer :: i,j,ierr,idummy
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
    integer :: MB_tmp,MB1_tmp,MB2_tmp
    character(10) :: cbuf
    logical :: flag
    include 'mpif.h'

! ---

    if ( SYStype /= 0 ) then
       write(*,*) "SYStype/=0 is not supported yet: SYSType=",SYStype
       call stop_program( "stop@read_io2" )
    end if

! ---

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

! ---

    call check_io2( flag )
    if ( .not.flag ) then
       info=-1
       return
    end if

!---

    call write_border( 0, " read_io2(start)" )

    info = 0

    if ( myrank == 0 ) then

       open(unit0,file=map_file_name,status='old')
       read(unit0,*) nprocs_new
       allocate( nmap(nprocs_new) ) ; nmap=0
       do i=1,nprocs_new
          read(unit0,*) nmap(i)
          do j=1,nmap(i)
             read(unit0,*) idummy
          end do
       end do
       nmax=maxval( nmap )
       allocate( iomap(0:12,nmax,nprocs_new) ) ; iomap=0
       rewind unit0
       read(unit0,*) nprocs_new
       do i=1,nprocs_new
          read(unit0,*) nmap(i)
          do j=1,nmap(i)
             read(unit0,*) iomap(0,j,i)
          end do
       end do
       do i=1,nprocs_new
          read(unit0,*)
          do j=1,nmap(i)
             read(unit0,*) cbuf,idummy,iomap(1:12,j,i)
          end do
       end do
       close(unit0)

    end if

!---

    call mpi_bcast(nmax,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nprocs_new,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( myrank /= 0 ) then
       allocate( nmap(nprocs_new) ) ; nmap=0
       allocate( iomap(0:12,nmax,nprocs_new) ) ; iomap=0
    end if
    call mpi_bcast(nmap,nprocs_new,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iomap,size(iomap),mpi_integer,0,mpi_comm_world,ierr)

!---

!    if ( myrank == 0 ) then
!       open(3,file="wf.dat1",status="old",form="unformatted")
!       read(3) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
!       read(3) MB_tmp,MB1_tmp,MB2_tmp
!       close(3)
!       write(*,*) "MB_0_IO,MB_1_IO=",1,MB_tmp
!    end if
!    call mpi_bcast(MB_tmp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!    MB_0_IO=1
!    MB_1_IO=MB_tmp

!---

!    call read_data_io2

! ---

    call write_border( 0, " read_io2(end)" )

  END SUBROUTINE read_io2


  SUBROUTINE check_io2( flag )
    implicit none
    logical,intent(OUT) :: flag
    integer :: ierr
    include 'mpif.h'
    if ( myrank == 0 ) inquire( FILE=map_file_name, EXIST=flag )
    call MPI_BCAST( flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
  END SUBROUTINE check_io2


  SUBROUTINE read_data_io2
    implicit none
    integer,parameter :: u1=10
    integer :: i,n,k,s,i1,i2,i3,i0,a1a,b1a,a2a,b2a,a3a,b3a
    integer :: a1c,b1c,a2c,b2c,a3c,b3c,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: a4a,b4a,a5a,b5a,a6a,b6a,a4b,b4b,a5b,b5b,a6b,b6b
    integer,allocatable :: LLL(:,:,:)
    character(5) :: crank_old
    character(32) :: filename
    real(8) :: ct0,ct1,et0,et1
#ifdef _DRSDFT_
    real(8),allocatable :: w(:,:,:)
    real(8),parameter :: z0=0.0d0
#else
    complex(8),allocatable :: w(:,:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
#endif
    logical :: disp_sw

    call write_border( 0, " read_data_io2(start)" )
    call check_disp_switch( disp_sw, 0 )

    a1b=Igrid(1,1)
    b1b=Igrid(2,1)
    a2b=Igrid(1,2)
    b2b=Igrid(2,2)
    a3b=Igrid(1,3)
    b3b=Igrid(2,3)
    a4b=MB_0_WF
    b4b=MB_1_WF
    a5b=MK_0_WF
    b5b=MK_1_WF
    a6b=MS_0_WF
    b6b=MS_1_WF
    allocate( LLL(a1b:b1b,a2b:b2b,a3b:b3b) ) ; LLL=0
    i=ML_0_WF-1
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    do i=1,nmap(myrank+1)
       write(crank_old,'(i5.5)') iomap(0,i,myrank+1)
       filename="wf.dat1."//crank_old
       if ( disp_sw ) write(*,*) filename
       open(u1,file=filename,form='unformatted')
       a1a=iomap( 1,i,myrank+1)
       b1a=iomap( 2,i,myrank+1)
       a2a=iomap( 3,i,myrank+1)
       b2a=iomap( 4,i,myrank+1)
       a3a=iomap( 5,i,myrank+1)
       b3a=iomap( 6,i,myrank+1)
       a4a=iomap( 7,i,myrank+1)
       b4a=iomap( 8,i,myrank+1)
       a5a=iomap( 9,i,myrank+1)
       b5a=iomap(10,i,myrank+1)
       a6a=iomap(11,i,myrank+1)
       b6a=iomap(12,i,myrank+1)
       allocate( w(a1a:b1a,a2a:b2a,a3a:b3a) ) ; w=z0
       a1c=max(a1b,a1a)
       b1c=min(b1b,b1a)
       a2c=max(a2b,a2a)
       b2c=min(b2b,b2a)
       a3c=max(a3b,a3a)
       b3c=min(b3b,b3a)
       do s=a6a,b6a
       do k=a5a,b5a
       do n=a4a,b4a
          w=z0
          read(u1) w
          if ( ( a4b <= n .and. n <= b4b ) &
          .and.( a5b <= k .and. k <= b5b ) &
          .and.( a6b <= s .and. s <= b6b ) ) then
             do i3=a3c,b3c
             do i2=a2c,b2c
             do i1=a1c,b1c
                i0=LLL(i1,i2,i3)
                unk(i0,n,k,s)=w(i1,i2,i3)
             end do
             end do
             end do
          end if
       end do ! n
       end do ! k
       end do ! s
       deallocate( w )
       close(u1)
    end do ! i

    deallocate( LLL )
    deallocate( iomap, nmap )

    call write_border( 0," read_data_io2(end)" )

  END SUBROUTINE read_data_io2


  SUBROUTINE read_data2_io2( string_in_1, string_in_2, u, v, nb0, nb1 )
    implicit none
    character(*),intent(IN) :: string_in_1, string_in_2
    integer,optional,intent(IN) :: nb0,nb1
    integer,parameter :: u1=10
    integer :: i,j,n,k,s,i1,i2,i3,i0,a1a,b1a,a2a,b2a,a3a,b3a
    integer :: a1c,b1c,a2c,b2c,a3c,b3c,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: a4a,b4a,a5a,b5a,a6a,b6a,a4b,b4b,a5b,b5b,a6b,b6b
    integer :: n0,k0,s0,nb_min,nb_max
    integer,allocatable :: LLL(:,:,:)
    character(5) :: crank_old
    character(32) :: filename
    real(8) :: ct0,ct1,et0,et1
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: u(:,:,:,:)
    real(8),optional,intent(INOUT) :: v(:,:,:,:)
    real(8),allocatable :: w(:,:,:)
    real(8),parameter :: z0=0.0d0
#else
    complex(8),intent(INOUT),allocatable :: u(:,:,:)
    complex(8),optional,intent(INOUT),allocatable :: v(:,:,:)
    complex(8),allocatable :: w(:,:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
#endif
    logical :: disp_sw

    call write_border( 0, " read_data2_io2(start)" )
    call check_disp_switch( disp_sw, 0 )

    nb_min=1
    nb_max=MB_WF
    if ( present(nb0) ) nb_min=nb0
    if ( present(nb1) ) nb_max=nb1

    a1b=Igrid(1,1)
    b1b=Igrid(2,1)
    a2b=Igrid(1,2)
    b2b=Igrid(2,2)
    a3b=Igrid(1,3)
    b3b=Igrid(2,3)
    a4b=MB_0_WF
    b4b=MB_1_WF
    a5b=MK_0_WF
    b5b=MK_1_WF
    a6b=MS_0_WF
    b6b=MS_1_WF
    allocate( LLL(a1b:b1b,a2b:b2b,a3b:b3b) ) ; LLL=0
    i=ML_0_WF-1
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    do i=1,nmap(myrank+1)
       write(crank_old,'(i5.5)') iomap(0,i,myrank+1)
       filename=string_in_1//crank_old//string_in_2
       if ( disp_sw ) write(*,*) filename
       open(u1,file=filename,form='unformatted')
       a1a=iomap( 1,i,myrank+1)
       b1a=iomap( 2,i,myrank+1)
       a2a=iomap( 3,i,myrank+1)
       b2a=iomap( 4,i,myrank+1)
       a3a=iomap( 5,i,myrank+1)
       b3a=iomap( 6,i,myrank+1)
       a4a=iomap( 7,i,myrank+1)
       b4a=iomap( 8,i,myrank+1)
       a5a=iomap( 9,i,myrank+1)
       b5a=iomap(10,i,myrank+1)
       a6a=iomap(11,i,myrank+1)
       b6a=iomap(12,i,myrank+1)
       allocate( w(a1a:b1a,a2a:b2a,a3a:b3a) ) ; w=z0
       a1c=max(a1b,a1a)
       b1c=min(b1b,b1a)
       a2c=max(a2b,a2a)
       b2c=min(b2b,b2a)
       a3c=max(a3b,a3a)
       b3c=min(b3b,b3a)
       do s=a6a,b6a
       do k=a5a,b5a
       do n=a4a,b4a
          if ( n < nb_min .or. n > nb_max ) cycle
          w=z0
          do i3=a3a,b3a
          do i2=a2a,b2a
          do i1=a1a,b1a
             read(u1) w(i1,i2,i3)
          end do
          end do
          end do
          if ( ( a4b <= n .and. n <= b4b ) &
          .and.( a5b <= k .and. k <= b5b ) &
          .and.( a6b <= s .and. s <= b6b ) ) then
             do i3=a3c,b3c
             do i2=a2c,b2c
             do i1=a1c,b1c
                i0=LLL(i1,i2,i3)-ML_0_wf+1
                n0=n-MB_0_WF+1
                k0=k-MK_0_WF+1
                s0=s-MS_0_WF+1
                u(i0,n0,k0,s0)=w(i1,i2,i3)
             end do
             end do
             end do
          end if
       end do ! n
       end do ! k
       end do ! s

       if ( present(v) ) then
          rewind(u1)
          do s=a6a,b6a
          do k=a5a,b5a
          do n=a4a,b4a
             if ( n < nb_min .or. n > nb_max ) cycle
             do i3=a3a,b3a
             do i2=a2a,b2a
             do i1=a1a,b1a
                read(u1) w(i1,i2,i3)
             end do
             end do
             end do
          end do
          end do
          end do
          do s=a6a,b6a
          do k=a5a,b5a
          do n=a4a,b4a
             if ( n < nb_min .or. n > nb_max ) cycle
             w=z0
             do i3=a3a,b3a
             do i2=a2a,b2a
             do i1=a1a,b1a
                read(u1) w(i1,i2,i3)
             end do
             end do
             end do
             if ( ( a4b <= n .and. n <= b4b ) &
             .and.( a5b <= k .and. k <= b5b ) &
             .and.( a6b <= s .and. s <= b6b ) ) then
                do i3=a3c,b3c
                do i2=a2c,b2c
                do i1=a1c,b1c
                   i0=LLL(i1,i2,i3)-ML_0_wf+1
                   n0=n-MB_0_WF+1
                   k0=k-MK_0_WF+1
                   s0=s-MS_0_WF+1
                   v(i0,n0,k0,s0)=w(i1,i2,i3)
                end do
                end do
                end do
             end if
          end do ! n
          end do ! k
          end do ! s
       end if

       deallocate( w )
       close(u1)
    end do ! i

    deallocate( LLL )
    deallocate( iomap, nmap )

    call write_border( 0," read_data2_io2(end)" )

  END SUBROUTINE read_data2_io2


END MODULE io2_module
