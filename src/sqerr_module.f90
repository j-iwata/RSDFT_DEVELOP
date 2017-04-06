MODULE sqerr_module

  use parallel_module, only: comm_grid, comm_spin

  implicit none

  PRIVATE
  PUBLIC :: calc_sqerr

  real(8),allocatable :: fold(:,:,:)

CONTAINS


  SUBROUTINE calc_sqerr( m, n, nmax, ndat, f_in, sqerr_out )
    implicit none
    integer,intent(IN)  :: m,n,nmax,ndat
    real(8),intent(IN)  :: f_in(m,n,ndat)
    real(8),intent(OUT) :: sqerr_out(nmax*ndat)
    integer :: ierr,i,n0,nn
    real(8),allocatable :: f(:,:,:)
    include 'mpif.h'

    call write_border( 1, " calc_sqerr(start)" )

    allocate( f(m,nmax,ndat) ) ; f=0.0d0

    n0=size( f_in(:,:,1) )
    nn=size( f(:,:,1) )
    do i=1,ndat
       call mpi_allgather( f_in(1,1,i),n0,MPI_REAL8 &
                          ,f(1,1,i),nn,MPI_REAL8,comm_spin,ierr )
    end do

    if ( .not.allocated(fold) ) then
       allocate( fold(m,nmax,ndat) ) ; f=0.0d0
       sqerr_out=1.d100
    else
       call calc_sqerr_sub( nmax, ndat, f, sqerr_out )
    end if

    fold=f

    deallocate( f )

    call write_border( 1, " calc_sqerr(end)" )

    return
  END SUBROUTINE calc_sqerr


  SUBROUTINE calc_sqerr_sub( nmax, ndat, f, sqerr_out )

    implicit none
    integer,intent(IN)  :: nmax, ndat
    real(8),intent(IN)  :: f(:,:,:)
    real(8),intent(OUT) :: sqerr_out(:)
    real(8) :: err0(nmax,ndat),err1(nmax,ndat)
    integer :: i,j,mm,m
    include 'mpif.h'

    m=size(f,1)
    call MPI_ALLREDUCE( m,mm,1,MPI_INTEGER,MPI_SUM,comm_grid,i )

    do j=1,ndat
    do i=1,nmax
       err0(i,j) = sum( (f(:,i,j)-fold(:,i,j))**2 )/mm
    end do
    end do
    call mpi_allreduce(err0,err1,nmax*ndat,MPI_REAL8,MPI_SUM,comm_grid,i)

    sqerr_out(:)=0.0d0
    m=0
    do j=1,ndat
    do i=1,nmax
       m=m+1
       sqerr_out(m) = err1(i,j)
    end do
    end do

  END SUBROUTINE calc_sqerr_sub


END MODULE sqerr_module
