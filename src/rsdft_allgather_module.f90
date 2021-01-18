module rsdft_allgather_module

  implicit none

  private
  public :: d_rsdft_allgatherv_div

  integer :: nblock_default=4
  integer :: n_opt, n_opt_h

contains

  subroutine d_rsdft_allgatherv_div( n, a, ir, id, comm, nblk_in )
    implicit none
    integer,intent(in) :: n
    real(8),intent(inout) :: a(n)
    integer,intent(in) :: ir(0:), id(0:)
    integer,intent(in) :: comm
    integer,intent(in) :: nblk_in
    integer :: nblk
    logical :: disp_sw
    integer :: i0,i1,nprc,mrnk,ierr,p
    integer :: nmax,ndat,i,j
    integer,allocatable :: irr(:),idd(:),id0(:),id1(:)
    real(8),allocatable :: tmp(:)
    include 'mpif.h'

    call write_border( 1, " d_rsdft_allgatherv_div(start)" )

    call check_disp_switch( disp_sw, 0 )

    nblk = nblk_in
    nprc = size(ir)
    call MPI_Comm_rank( comm, mrnk, ierr )

    allocate( tmp(nblk*nprc) ); tmp=0.0d0
    allocate( irr(0:nprc-1)  ); irr=0
    allocate( idd(0:nprc-1)  ); idd=0
    allocate( id0(0:nprc-1)  ); id0=0
    allocate( id1(0:nprc-1)  ); id1=0

    id0(:) = id(:)

    nmax = maxval(ir)
    do i = 1, nmax, nblk
      do p=0,nprc-1
        id1(p) = min( id0(p)+nblk, id(p)+ir(p) ) - 1
        irr(p) = id1(p) - id0(p) + 1
        idd(p) = sum(irr(0:p))-irr(p)
      end do
      call MPI_Allgatherv( a(id0(mrnk)+1),irr(mrnk),MPI_REAL8 &
                          ,tmp,irr,idd,MPI_REAL8,comm,ierr )
      do p=0,nprc-1
        if ( p /= mrnk ) then
          do j=1,irr(p)
            a(id0(p)+j)=tmp(idd(p)+j)
          end do
        end if        
        id0(p) = id1(p) + 1
      end do

    end do

    deallocate( id1 )
    deallocate( id0 )
    deallocate( idd )
    deallocate( irr )
    deallocate( tmp )

    call write_border( 1, " d_rsdft_allgatherv_div(end)" )
  end subroutine d_rsdft_allgatherv_div

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
