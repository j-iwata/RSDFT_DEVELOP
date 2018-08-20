module calc_overlap_bp_module

  use parallel_module, only: nprocs_b, myrank_b, comm_band, myrank
  use rsdft_mpi_module
  use calc_overlap_module, only: calc_overlap_no_mpi
  use watch_module

  implicit none

  private
  public :: calc_overlap_bp
  public :: mochikae_matrix

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
    real(8),allocatable :: sendbuf(:,:), recvbuf(:,:), ab_blk(:,:)
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr, itags, nreq
    logical,allocatable :: ttt(:,:)
    real(8) :: ttmp(2),tttt(2,13)

    call write_border( 1, "calc_overlap_bp(start)" )

    !call watchb( ttmp ); tttt=0.0d0

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
    allocate( ab_blk(blk_size,blk_size) ); ab_blk=0.0d0

    !call watchb( ttmp, tttt(:,1) )

    do iblk=1,nblk

       !call watchb( ttmp )

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

       !call watchb( ttmp, tttt(:,2) )

       do istep=1,nprocs_b

          !call watchb( ttmp )

          if ( istep < nprocs_b ) then
             call MPI_Irecv( recvbuf, m*blk_size, MPI_REAL8, jrank, itags, comm_band, ireq(1), ierr )
             call MPI_Isend( sendbuf, m*blk_size, MPI_REAL8, irank, itags, comm_band, ireq(2), ierr )
             nreq=2
          else if ( istep == nprocs_b ) then
             nreq=0
          end if

          !call watchb( ttmp, tttt(:,3) )

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

          !call watchb( ttmp, tttt(:,4) )

! --- (1)
!
!          do j=j0,j1
!             do i=i0,i1
!                ab(i,j) = sum( sendbuf(:,i-i0+1)*b(:,j-j0+k0) )*alpha
!             end do
!          end do
!
! --- (2)
!
          if ( istep == 1 ) then

             !call watchb( ttmp )

             call calc_overlap_no_mpi( sendbuf, b(:,k0:k1), alpha, ab_blk )

             !call watchb( ttmp, tttt(:,5) )

             do j=j0,j1
                do i=i0,i1
                   if ( i >= j ) ab(i,j)=ab_blk(i-i0+1,j-j0+1)
                end do
             end do

             !call watchb( ttmp, tttt(:,6) )

          else

             !call watchb( ttmp )

             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab(i0,j0),nb)
!             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab_blk,blk_size)
!             ab(i0:i1,j0:j1) = ab_blk(:,:)

             !call watchb( ttmp, tttt(:,7) )

          end if
!
! -------

          !call watchb( ttmp )

          if ( nreq == 2 ) then
             call MPI_Waitall( 2, ireq, istatus, ierr )
             sendbuf(:,:) = recvbuf(:,:)
          end if

          !call watchb( ttmp, tttt(:,8) )

       end do ! istep

       !call watchb( ttmp )

       if ( iblk == 2 ) then
          i0 = b0 + myrank_b*blk_size
          i1 = i0 + blk_size - 1
          k0 = blk_size + 1
          j0 = myrank_b*blk_size + 1
          j1 = j0 + blk_size - 1
! --- (1)
!          do j=j0,j1
!             do i=i0,i1
!                ab(i,j) = sum( a(:,i-i0+k0)*b(:,j-j0+1) )*alpha
!             end do
!          end do
! --- (2)
          call DGEMM('T','N',blk_size,blk_size,m,alpha,a(1,k0),m,b,m,0.0d0,ab(i0,j0),nb)
! -------
       end if

       !call watchb( ttmp, tttt(:,9) )

    end do ! iblk

    !call watchb( ttmp )

    call rsdft_allreduce_sum( ab, comm_band )

    !call watchb( ttmp, tttt(:,10) )

    do j=1,nb
       do i=j+1,nb
          if ( ab(i,j) == 0.0d0 ) ab(i,j)=ab(j,i)
       end do
    end do

    !call watchb( ttmp, tttt(:,11) )

    do j=1,nb
       do i=1,j-1
          ab(i,j)=0.0d0
       end do
    end do

    !call watchb( ttmp, tttt(:,12) )

    !ab=ab*alpha

    deallocate( ab_blk  )
    deallocate( recvbuf )
    deallocate( sendbuf )

    !call watchb( ttmp, tttt(:,13) )

    !if ( myrank == 0 ) then
    !    do i=1,13
    !       write(*,'(1x,"time_calc_overlap_bp(",i2.2,")",2f10.5)') i,tttt(:,i)
    !    end do
    !end if

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


  subroutine mochikae_matrix( a, nblk )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,intent(in) :: nblk
    integer :: i,nb, blk_size, i0,i1,j0,j1,ib,irank_b
    real(8),allocatable :: b(:,:), c(:,:)
    real(8) :: ttmp(2),tttt(2,9)

    !call write_border( 1, "mochikae_matrix(start)" )

    !call watchb( ttmp ); tttt=0.0d0

    nb = size( a, 1 )

    allocate( b(nb,nb) ); b=0.0d0
    allocate( c(nb,nb) ); c=0.0d0

    !call watchb( ttmp, tttt(:,1) )

    call fill_matrix( a )

    !call watchb( ttmp, tttt(:,2) )

    blk_size = nb/(nblk*nprocs_b)

    j0 = 1
    j1 = blk_size
    do irank_b=0,nprocs_b-1
       do ib=1,nb,nb/nblk
          i0 = ib + irank_b*blk_size
          i1 = i0 + blk_size - 1
          do i=1,blk_size
             !b(i0+i-1,j0+i-1)=1.0d0
             b(:,j0+i-1)=a(:,i0+i-1)
          end do
          j0 = j0 + blk_size
          j1 = j1 + blk_size
       end do
    end do ! irank
    a=b
    j0 = 1
    j1 = blk_size
    do irank_b=0,nprocs_b-1
       do ib=1,nb,nb/nblk
          i0 = ib + irank_b*blk_size
          i1 = i0 + blk_size - 1
          do i=1,blk_size
             !b(i0+i-1,j0+i-1)=1.0d0
             b(j0+i-1,:)=a(i0+i-1,:)
          end do
          j0 = j0 + blk_size
          j1 = j1 + blk_size
       end do
    end do ! irank
    a=b

    !call watchb( ttmp, tttt(:,3) )

    !a = matmul( a, b )
    !c = transpose( b )
    !b = matmul( c, a )
    !a = b
    !call DGEMM('N','N',nb,nb,nb,1.0d0,a,nb,b,nb,0.0d0,c,nb)
    !call DGEMM('T','N',nb,nb,nb,1.0d0,b,nb,c,nb,0.0d0,a,nb)

    !call watchb( ttmp, tttt(:,4) )

    deallocate( c )
    deallocate( b )

    call cut_matrix( a )

    !call watchb( ttmp, tttt(:,5) )

    !if ( myrank == 0 ) then
    !   do i=1,5
    !      write(*,'(4x,"time_mochikae_matrix(",i1,")",2f10.5)') i, tttt(:,i)
    !   end do
    !end if

    !call write_border( 1, "mochikae_matrix(end)" )

  end subroutine mochikae_matrix


  subroutine fill_matrix( a )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer :: i,j
    do j=1,size(a,1)
       do i=1,j-1
          a(i,j)=a(j,i)
       end do
    end do
  end subroutine fill_matrix


  subroutine cut_matrix( a )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer :: i,j
    do j=1,size(a,1)
       do i=1,j-1
          a(i,j)=0.0d0
       end do
    end do
  end subroutine cut_matrix


end module calc_overlap_bp_module
