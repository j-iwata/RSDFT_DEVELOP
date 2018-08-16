module calc_overlap_bp_module

  use parallel_module, only: nprocs_b, myrank_b, comm_band
  use rsdft_mpi_module

  implicit none

  private
  public :: calc_overlap_bp
  public :: mochikae, mochikae_back

  integer :: nblk0, nblk1

contains

  subroutine calc_overlap_bp( nb, a, b, alpha, ab )
    implicit none
    integer,intent(in)  :: nb
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    include 'mpif.h'
    integer :: m,n,nblk,blk_size
    integer :: ib,i0,i1,j0,j1,i,j,iblk,k0,k1,b0,b1
    integer :: irank, jrank, istep, nstep
    real(8),allocatable :: sendbuf(:,:), recvbuf(:,:)
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr, itags
    logical,allocatable :: ttt(:,:)

    call write_border( 1, "calc_overlap_bp(start)" )

    ab(:,:)=0.0d0

    m = size( a, 1 )
    n = size( a, 2 )

    if ( mod(nprocs_b,2) == 0 ) then
       nblk = 2
    else
       nblk = 1
    end if

    blk_size = n/nblk

    allocate( sendbuf(m,blk_size) ); sendbuf=0.0d0
    allocate( recvbuf(m,blk_size) ); recvbuf=0.0d0

    do iblk=1,nblk

       b0 = (iblk-1)*nb/nblk + 1
       b1 = b0 + nb/nblk - 1

       j0 = b0 + myrank_b*blk_size
       j1 = j0 + blk_size - 1

       i0 = j0 - blk_size
       i1 = j1

       irank = myrank_b - 1; if ( irank <  0        ) irank = irank + nprocs_b
       jrank = myrank_b + 1; if ( jrank >= nprocs_b ) jrank = jrank - nprocs_b
       itags = 10

       k0 = (iblk-1)*blk_size + 1
       k1 = k0 + blk_size - 1

       sendbuf(:,:) = a(:,k0:k1)

       do istep=1,nprocs_b

          call MPI_Irecv( recvbuf, m*blk_size, MPI_REAL8, jrank, itags, comm_band, ireq(1), ierr )
          call MPI_Isend( sendbuf, m*blk_size, MPI_REAL8, irank, itags, comm_band, ireq(2), ierr )

          i0 = i0 + blk_size
          if ( i0 > b1 ) then
             i0 = b0
             j0 = j0 + nb/nblk
             k0 = k0 + blk_size
             if ( j0 > nb ) then
                j0 = 1 + myrank_b*blk_size
                k0 = 1
             end if
             j1 = j0 + blk_size - 1
             k1 = k0 + blk_size - 1
          end if
          i1 = i0 + blk_size - 1

          do j=j0,j1
             do i=i0,i1
                ab(i,j) = sum( sendbuf(:,i-i0+1)*b(:,j-j0+k0) )
             end do
          end do

          call MPI_Waitall( 2, ireq, istatus, ierr )

          sendbuf(:,:) = recvbuf(:,:)

       end do ! istep

       if ( iblk == 2 ) then
          i0 = b0 + myrank_b*blk_size
          i1 = i0 + blk_size - 1
          k0 = blk_size + 1
          j0 = myrank_b*blk_size + 1
          j1 = j0 + blk_size - 1
          do j=j0,j1
             do i=i0,i1
                ab(i,j) = sum( a(:,i-i0+k0)*b(:,j-j0+1) )
             end do
          end do
       end if

    end do ! iblk

    call rsdft_allreduce_sum( ab, comm_band )

    do j=1,nb
       do i=j+1,nb
          if ( ab(i,j) == 0.0d0 ) ab(i,j)=ab(j,i)
       end do
    end do

    ab=ab*alpha

    deallocate( recvbuf )
    deallocate( sendbuf )

    call write_border( 1, "calc_overlap_bp(end)" )

  end subroutine calc_overlap_bp


  subroutine calc_overlap_bp_1( a, b, alpha, ab )
    implicit none
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    integer :: m,n,j,k
    call write_border( 1, "calc_overlap_bp_1(start)" )
    nblk0 = size(a,2)
    nblk1 = 4
    call calc_overlap_sub( nblk0, a, b, ab )
    ab=ab*alpha
    call write_border( 1, "calc_overlap_bp_1(end)" )
  end subroutine calc_overlap_bp_1

  recursive subroutine calc_overlap_sub( nblk, a, b, c )
    implicit none
    integer,intent(in)    :: nblk
    real(8),intent(in)    :: a(:,:),b(:,:)
    real(8),intent(inout) :: c(:,:)
    integer :: m,n,i,j,i0,i1,j0,j1,ni,nj,nblkh
    real(8),allocatable :: ctmp(:,:)
    real(8),parameter :: one=1.0d0, zero=0.0d0

    m=size(a,1)
    n=size(a,2)

    do i0=1,n,nblk
       i1=min(i0+nblk-1,n)
       ni=i1-i0+1

       do j0=1,n,nblk
          j1=min(j0+nblk-1,n)
          nj=j1-j0+1

          if ( j0 > i1 ) then

             cycle

          else if ( j1 <= i0 ) then

             call dgemm &
                  ('T','N',ni,nj,m,one,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)

          else

             if ( ni > nblk1 ) then
                allocate( ctmp(ni,ni) ); ctmp=zero
                nblkh=nblk/2
                call calc_overlap_sub( nblkh, a(:,i0:i1), b(:,j0:j1), ctmp )
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
                do j=j0,j1
                do i=j ,i1
                    c(i,j)=sum( a(:,i)*b(:,j) )
                end do
                end do
             end if

          end if

       end do ! j0

    end do ! i0

  end subroutine calc_overlap_sub


  subroutine calc_overlap_bp_0( a, b, alpha, ab )
    implicit none
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    integer :: m,n,j,k
    call write_border( 1, "calc_overlap_bp_0(start)" )
    m=size(a,1)
    n=size(a,2)
    do k=1,n
       do j=k,n
          ab(j,k) = sum( a(:,j)*b(:,k) )*alpha
       end do
    end do
    call write_border( 1, "calc_overlap_bp_0(end)" )
  end subroutine calc_overlap_bp_0


  subroutine mochikae( f, nblk )
    implicit none
    real(8),intent(inout) :: f(:,:)
    integer,intent(in)    :: nblk
    real(8),allocatable :: g(:,:)
    integer :: m,n,nskip,ib,i0,i1,j0,j1
    integer :: blk_size

    m=size(f,1)
    n=size(f,2)

    allocate( g(m,n) ); g=0.0d0

    blk_size = n/(nblk*nprocs_b)

    j0 = 1
    j1 = blk_size

    do ib=1,n,n/nblk

       i0 = ib + myrank_b*blk_size
       i1 = i0 + blk_size - 1

       g(:,j0:j1) = f(:,i0:i1)

       j0 = j0 + blk_size
       j1 = j1 + blk_size

    end do

    f(:,:) = g(:,:)

    deallocate( g )

  end subroutine mochikae


  subroutine mochikae_back( f, nblk )
    implicit none
    real(8),intent(inout) :: f(:,:)
    integer,intent(in)    :: nblk
    real(8),allocatable :: g(:,:)
    integer :: m,n,nskip,ib,i0,i1,j0,j1
    integer :: blk_size

    m=size(f,1)
    n=size(f,2)

    allocate( g(m,n) ); g=0.0d0

    blk_size = n/(nblk*nprocs_b)

    j0 = 1
    j1 = blk_size

    do ib=1,n,n/nblk

       i0 = ib + myrank_b*blk_size
       i1 = i0 + blk_size - 1

       g(:,i0:i1) = f(:,j0:j1)

       j0 = j0 + blk_size
       j1 = j1 + blk_size

    end do

    f(:,:) = g(:,:)

    deallocate( g )

  end subroutine mochikae_back


end module calc_overlap_bp_module
