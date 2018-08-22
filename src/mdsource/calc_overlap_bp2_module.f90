module calc_overlap_bp2_module

  use parallel_module, only: nprocs_b, myrank_b, comm_band, myrank, comm_grid
  use rsdft_mpi_module
  use calc_overlap_module, only: calc_overlap_no_mpi
  use watch_module

  implicit none

  private
  public :: calc_overlap_bp2
  public :: mochikae_matrix2

  integer :: nblk0, nblk1

contains

  subroutine calc_overlap_bp2( nb, a, b, alpha, ab )
    implicit none
    integer,intent(in)  :: nb
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    include 'mpif.h'
    integer :: m,n,nblk,blk_size
    integer :: ib,i0,i1,j0,j1,i,j,iblk,k0,k1,b0,b1,l0,l1
    integer :: irank, jrank, istep, nstep
    real(8),allocatable :: sendbuf(:,:), recvbuf(:,:), ab_blk(:,:)
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr, itags, nreq
    logical,allocatable :: ttt(:,:)
    real(8) :: ttmp(2),tttt(2,14)

    call write_border( 1, "calc_overlap_bp2(start)" )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp ); tttt=0.0d0

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

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,1) )

    do iblk=1,nblk

       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call watchb( ttmp )

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

       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call watchb( ttmp, tttt(:,2) )

       do istep=1,nprocs_b

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp )

          if ( istep < nprocs_b ) then
             call MPI_Irecv( recvbuf, m*blk_size, MPI_REAL8, jrank, itags, comm_band, ireq(1), ierr )
             call MPI_Isend( sendbuf, m*blk_size, MPI_REAL8, irank, itags, comm_band, ireq(2), ierr )
             nreq=2
          else if ( istep == nprocs_b ) then
             nreq=0
          end if

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,3) )

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

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,4) )

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

             call MPI_Barrier( MPI_COMM_WORLD, ierr )
             call watchb( ttmp )

             !call calc_overlap_no_mpi( sendbuf, b(:,k0:k1), alpha, ab_blk )
             !call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab(i0,j0),nb)
             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab_blk,blk_size)

             call MPI_Barrier( MPI_COMM_WORLD, ierr )
             call watchb( ttmp, tttt(:,5) )

             do j=j0,j1
                do i=i0,i1
                   if ( i >= j ) ab(i,j)=ab_blk(i-i0+1,j-j0+1)
                end do
             end do

             call MPI_Barrier( MPI_COMM_WORLD, ierr )
             call watchb( ttmp, tttt(:,6) )

          else

             call MPI_Barrier( MPI_COMM_WORLD, ierr )
             call watchb( ttmp )

             !call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab(i0,j0),nb)
             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab_blk,blk_size)
             ab(i0:i1,j0:j1) = ab_blk(:,:)

             call MPI_Barrier( MPI_COMM_WORLD, ierr )
             call watchb( ttmp, tttt(:,7) )

          end if
!
! -------

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp )

!          if ( i0 < j0 ) then
!             do j=j0,j1
!                do i=i0,i1
!                   ab(j,i) = ab(i,j)
!                end do
!             end do
!          else if ( i0 > j0 ) then
             do j=j0,j1
                do i=i0,i1
                   ab(j,i) = ab(i,j)
                end do
             end do
!          end if

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,8) )

          if ( nreq == 2 ) then
             call MPI_Waitall( 2, ireq, istatus, ierr )
             sendbuf(:,:) = recvbuf(:,:)
          end if

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,9) )

       end do ! istep

       if ( iblk == 2 ) then

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp )

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

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,10) )

          do j=j0,j1
             do i=i0,i1
                ab(j,i) = ab(i,j)
             end do
          end do

          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          call watchb( ttmp, tttt(:,11) )

       end if

    end do ! iblk

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp )

! --- (1)
!
    call rsdft_allreduce_sum( ab, comm_band )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,12) )

    !do j=1,nb
    !   do i=j+1,nb
    !      if ( ab(i,j) == 0.0d0 ) ab(i,j)=ab(j,i)
    !   end do
    !end do
    !!call watchb( ttmp, tttt(:,11) )
    !do j=1,nb
    !   do i=1,j-1
    !      ab(i,j)=0.0d0
    !   end do
    !end do
    !!call watchb( ttmp, tttt(:,12) )
!
! --- (2)
!
#ifdef test
    deallocate( sendbuf )
    allocate( sendbuf(nb+blk_size,nb/nblk) ); sendbuf=0.0d0

    j0 = nb/nblk + 1 + myrank_b*blk_size
    j1 = j0 + blk_size - 1

    k0 = 1 + myrank_b*blk_size
    k1 = k0 + blk_size - 1

    i0 = 1
    i1 = myrank_b*blk_size
    if ( i0 <= i1 ) then
       do i=i0,i1,blk_size
          ab_blk(:,:) = ab(i:i+blk_size-1,j0:j1)
          ab(i:i+blk_size-1,j0:j1) = 0.0d0
          sendbuf(i:i+blk_size-1,k0:k1) = transpose( ab_blk )
       end do
    end if

    i0 = 1 + myrank_b*blk_size
    i1 = i0 + nb/nblk + blk_size - 1
    sendbuf(i0:i1,k0:k1) = ab(i0:i1,k0:k1)

    i0 = nb/nblk + (myrank_b+1)*blk_size + 1
    i1 = nb
    if ( i0 <= i1 ) then
       sendbuf(i0:i1,k0:k1) = ab(i0:i1,j0:j1)
    end if

    i0 = nb/nblk + 1 + myrank_b*blk_size
    i1 = i0 + blk_size - 1
    sendbuf(nb+1:nb+blk_size,k0:k1) = ab(i0:i1,i0:i1)

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,12) )

    call rsdft_allreduce_sum( sendbuf(:,k0:k1), comm_grid )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,13) )

    call rsdft_allgather( sendbuf(:,k0:k1), sendbuf, comm_band )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,12) )

    do irank=0,nprocs_b-1

       j0 = nb/nblk + 1 + irank*blk_size
       j1 = j0 + blk_size - 1

       k0 = 1 + irank*blk_size
       k1 = k0 + blk_size - 1

       i0 = nb/nblk + 1 + irank*blk_size
       i1 = i0 + blk_size - 1
       ab(i0:i1,i0:i1) = sendbuf(nb+1:nb+blk_size,k0:k1)

       i0 = nb/nblk + (irank+1)*blk_size + 1
       i1 = nb
       if ( i0 <= i1 ) then
          ab(i0:i1,j0:j1) = sendbuf(i0:i1,k0:k1)
       end if

       i0 = 1 + irank*blk_size
       i1 = i0 + nb/nblk + blk_size - 1
       ab(i0:i1,k0:k1) = sendbuf(i0:i1,k0:k1)

       i0 = 1
       i1 = irank*blk_size
       if ( i0 <= i1 ) then

          l0=nb/nblk+irank*blk_size+1
          l1=l0+blk_size-1
          do i=i0,i1,blk_size
             j0=i-i0+1
             j1=j0+blk_size-1
             ab(l0:l1,j0:j1) = sendbuf(i:i+blk_size-1,k0:k1)
          end do

       end if

    end do ! irank

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,13) )
#endif
!
! -------

    !ab=ab*alpha

    deallocate( ab_blk  )
    deallocate( recvbuf )
    deallocate( sendbuf )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,14) )

    !if ( myrank == 0 ) then
    !    do i=1,14
    !       write(*,'(1x,"time_calc_overlap_bp(",i2.2,")",2f10.5)') i,tttt(:,i)
    !    end do
    !end if

    call write_border( 1, "calc_overlap_bp(end)" )

  end subroutine calc_overlap_bp2


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


  subroutine mochikae_matrix2( a, nblk )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,intent(in) :: nblk
    integer :: i,j,nb, blk_size, i0,i1,j0,j1,ib,irank_b,ierr
    real(8),allocatable :: b(:,:)
    real(8) :: ttmp(2),tttt(2,8)
    include 'mpif.h'

    call write_border( 1, "mochikae_matrix2(start)" )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp ); tttt=0.0d0

    nb = size( a, 1 )

    allocate( b(nb,nb) ); b=0.0d0

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,1) )

!    call fill_matrix( a, nblk )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,2) )

    blk_size = nb/(nblk*nprocs_b)

! --- (1)
!
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
!    a=b
    j0 = 1
    j1 = blk_size
    do irank_b=0,nprocs_b-1
       do ib=1,nb,nb/nblk
          i0 = ib + irank_b*blk_size
          i1 = i0 + blk_size - 1
          do j=1,nb
          do i=1,blk_size
             !b(i0+i-1,j0+i-1)=1.0d0
             !b(j0+i-1,:)=a(i0+i-1,:)
             a(j0+i-1,j)=b(i0+i-1,j)
          end do
          end do
          j0 = j0 + blk_size
          j1 = j1 + blk_size
       end do
    end do ! irank
!    a=b
!
! --- (2)
#ifdef test
    j0 = 1 + myrank_b*nb/nprocs_b
    j1 = j0 + blk_size - 1
    do ib=1,nb,nb/nblk
       i0 = ib + myrank_b*blk_size
       i1 = i0 + blk_size - 1
       do i=1,blk_size
          b(:,j0+i-1)=a(:,i0+i-1)
       end do
       j0 = j0 + blk_size
       j1 = j1 + blk_size
    end do

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,3) )

    j0 = 1 + myrank_b*nb/nprocs_b
    j1 = j0 + nb/nprocs_b - 1
    call rsdft_allgather( b(:,j0:j1), b, comm_band )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,4) )

    a=0.0d0

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,5) )

    j0 = 1 + myrank_b*nb/nprocs_b
    j1 = j0 + blk_size - 1
    do ib=1,nb,nb/nblk
       i0 = ib + myrank_b*blk_size
       i1 = i0 + blk_size - 1
       do j=1,nb
       do i=1,blk_size
          a(j0+i-1,j)=b(i0+i-1,j)
       end do
       end do
       j0 = j0 + blk_size
       j1 = j1 + blk_size
    end do

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,6) )

    call rsdft_allreduce_sum( a, comm_band )
#endif
! -------

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,7) )

    deallocate( b )

    call cut_matrix( a )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,8) )

    !if ( myrank == 0 ) then
    !   do i=1,8
    !      write(*,'(4x,"time_mochikae_matrix(",i1,")",2f10.5)') i, tttt(:,i)
    !   end do
    !end if

    call write_border( 1, "mochikae_matrix2(end)" )

  end subroutine mochikae_matrix2


  subroutine fill_matrix( a, nblk )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,intent(in) :: nblk
    integer :: i,j
    integer :: iblk,j0,j1,nb,blk_size
!
! --- (1)
!
!    do j=1,size(a,1)
!       do i=1,j-1
!          a(i,j)=a(j,i)
!       end do
!    end do
!
! --- (2)
!
    nb = size(a,1)
    blk_size = nb/(nprocs_b*nblk)
    do iblk=1,nblk
       j0 = 1 + myrank_b*blk_size + (iblk-1)*nb/nblk
       j1 = j0 + blk_size - 1
       do j=j0,j1
          do i=1,j-1
             a(i,j)=a(j,i)
          end do
       end do
    end do
!
! -------
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


end module calc_overlap_bp2_module
