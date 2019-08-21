module rsdft_allgather_module

  implicit none

  private
  public :: d_rsdft_allgather
  public :: test_allgather

  integer :: nblock_default=4
  integer :: n_opt, n_opt_h

contains


  subroutine d_rsdft_allgather( a, b, comm, ierr, nblock_in )
    implicit none
    real(8),intent(in)  :: a(:)
    real(8),intent(out) :: b(:)
    integer,intent(in)  :: comm
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: nblock_in
    integer :: na, nb, nprocs, nblock
    integer :: i,i0,i1,mm,p
    real(8),allocatable :: c(:)
    include 'mpif.h'

    call MPI_Comm_size( comm, nprocs, ierr )

    na=size(a)
    nb=size(b)

    if ( present(nblock_in) ) then
       nblock = nblock_in
    else
       nblock = nblock_default
    end if

    allocate( c(nblock*nprocs) ); c=0.0d0

    do i=1,na,nblock

       i0 = i
       i1 = min(i0+nblock-1,na)
       mm = i1-i0+1

       call MPI_Allgather(a(i0),mm,MPI_REAL8,c,mm,MPI_REAL8,comm,ierr)

       do p=0,nprocs-1
          b(p*na+i0:p*na+i1) = c(p*mm+1:p*mm+mm)
       end do

    end do

    deallocate( c )

  end subroutine d_rsdft_allgather


  subroutine test_allgather
    implicit none
    integer :: i,j,k,m,mt,n,i0,ierr,myrank,nprocs,npow,p
    real(8),allocatable :: a(:), b(:), c(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
    include 'mpif.h'

    !return

    call write_border( 0, " test_allgather(start)" )

    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)

    npow=16

    n=2**npow
    allocate( a(n) ); a(:)=myrank+1
    allocate( b(n*nprocs) ); b(:)=0.0d0
    allocate( c(n*nprocs) ); c(:)=0.0d0
    ctmin=1.d100
    etmin=1.d100
    do j=npow,0,-1
       m = 2**j
       mt=0
       do i=1,2**(npow-j)
          mt=mt+m
       end do
       b=0.0d0
       c=0.0d0
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct0) ; et0=mpi_wtime()
       do k=1,10
          do i=1,2**(npow-j)
             i0=(i-1)*m+1
             call MPI_Allgather(a(i0),m,MPI_REAL8,b,m,MPI_REAL8,MPI_COMM_WORLD,ierr)
             do p=0,nprocs-1
                c(p*mt+i0:p*mt+i0+m-1)=b(p*m+1:p*m+m)
             end do
          end do
       end do ! k
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct1) ; et1=mpi_wtime()
       ct=ct1-ct0
       et=et1-et0
       ctmin=min(ct,ctmin)
       if ( et < etmin ) then
          n_opt = m
          etmin = et
       end if
       if ( myrank == 0 ) then
          write(*,'(1x,3i12,4f10.5,2x,i8,g15.8)') m,i-1,mt,ct,ctmin,et,etmin,count(c/=0.0d0),sum(c)
          !do p=0,nprocs-1
          !   write(*,*) p, count(nint(c)==p+1)
          !end do
       end if
    end do ! j
    deallocate( c )
    deallocate( b )
    deallocate( a )
    call MPI_BCAST(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h
    call write_border( 0, " test_allgather(end)" )
  end subroutine test_allgather


end module rsdft_allgather_module
