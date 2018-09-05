module calc_overlap_bp_module

  use parallel_module, only: nprocs_b, myrank_b, comm_band, myrank, comm_grid
  use rsdft_mpi_module
  use calc_overlap_module, only: calc_overlap_no_mpi, calc_overlap
  use watch_module
  use transpose_module, only: rsdft_transpose

  implicit none

  private
  public :: calc_overlap_bp

  real(8),allocatable :: sendbuf(:,:), recvbuf(:,:)
  real(8),allocatable :: ab_blk(:,:), tr_blk(:,:)
  real(8),allocatable :: ab_tmp(:,:)

  interface calc_overlap_bp
     module procedure calc_overlap_bp &
          , calc_overlap_bp_zd, calc_overlap_bp_zz
  end interface

contains

  subroutine calc_overlap_bp( nb, a, b, alpha, ab )
    implicit none
    integer,intent(in)  :: nb
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    include 'mpif.h'
    integer :: m,n,nblk,blk_size
    integer :: ib,i0,i1,j0,j1,i,j,iblk,k0,k1,b0,b1,l0,l1
    integer :: irank, jrank, istep, nstep
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr, itags, nreq
    real(8) :: ttmp(2),tttt(2,15)

    if ( nprocs_b == 1 ) then
       call calc_overlap( size(a,1), nb, a, b, alpha, ab )
       return
    end if

    !call write_border( 1, "calc_overlap_bp(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    ab(:,:)=0.0d0

    if ( mod(nprocs_b,2) == 0 ) then
       nblk = 2
    else
       nblk = 1
    end if

    m = size( a, 1 )
    n = size( a, 2 )
    blk_size = n/nblk

    if ( .not.allocated(sendbuf) ) then
       allocate( sendbuf(m,blk_size) ); sendbuf=0.0d0
       allocate( recvbuf(m,blk_size) ); recvbuf=0.0d0
       allocate( ab_blk(blk_size,blk_size) ); ab_blk=0.0d0
       allocate( tr_blk(blk_size,blk_size) ); tr_blk=0.0d0
    end if

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    do iblk=1,nblk

       !call watchb( ttmp, barrier="on" )

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

       !call watchb( ttmp, tttt(:,2), barrier="on" )

       do istep=1,nprocs_b

          !call watchb( ttmp, barrier="on" )

          if ( istep < nprocs_b ) then
             call MPI_Irecv( recvbuf, m*blk_size, MPI_REAL8, jrank, itags, comm_band, ireq(1), ierr )
             call MPI_Isend( sendbuf, m*blk_size, MPI_REAL8, irank, itags, comm_band, ireq(2), ierr )
             nreq=2
          else if ( istep == nprocs_b ) then
             nreq=0
          end if

          !call watchb( ttmp, tttt(:,3), barrier="on" )

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

          !call watchb( ttmp, tttt(:,4), barrier="on" )

          if ( istep == 1 ) then

             !call watchb( ttmp, barrier="on" )
! ---(1)
!             call calc_overlap_no_mpi( sendbuf, b(:,k0:k1), alpha, ab_blk )
!             call rsdft_transpose( ab_blk, 16 )
! ---(2)
             !call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab(i0,j0),nb)
             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab_blk,blk_size)
! ------
             !call watchb( ttmp, tttt(:,5), barrier="on" )

             call rsdft_allreduce_sum( ab_blk, comm_grid )

             ab(i0:i1,j0:j1) = ab_blk(:,:)

             !call watchb( ttmp, tttt(:,6), barrier="on" )

          else

             !call watchb( ttmp, barrier="on" )

             !call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab(i0,j0),nb)
             call DGEMM('T','N',blk_size,blk_size,m,alpha,sendbuf,m,b(1,k0),m,0.0d0,ab_blk,blk_size)

             !call watchb( ttmp, tttt(:,7), barrier="on" )

             call rsdft_allreduce_sum( ab_blk, comm_grid )

             ab(i0:i1,j0:j1) = ab_blk(:,:)

             !call watchb( ttmp, tttt(:,8), barrier="on" )

          end if

          !call watchb( ttmp, barrier="on" )

          if ( istep > 1 ) then
             tr_blk = transpose( ab_blk )
             ab(j0:j1,i0:i1) = tr_blk(:,:)
          end if

          !call watchb( ttmp, tttt(:,9), barrier="on" )

          if ( nreq == 2 ) then
             call MPI_Waitall( 2, ireq, istatus, ierr )
             sendbuf(:,:) = recvbuf(:,:)
          end if

          !call watchb( ttmp, tttt(:,10), barrier="on" )

       end do ! istep

       if ( iblk == 2 ) then

          !call watchb( ttmp, barrier="on" )

          i0 = b0 + myrank_b*blk_size
          i1 = i0 + blk_size - 1
          k0 = blk_size + 1
          j0 = myrank_b*blk_size + 1
          j1 = j0 + blk_size - 1

!          call DGEMM('T','N',blk_size,blk_size,m,alpha,a(1,k0),m,b,m,0.0d0,ab(i0,j0),nb)
          call DGEMM('T','N',blk_size,blk_size,m,alpha,a(1,k0),m,b,m,0.0d0,ab_blk,blk_size)

          !call watchb( ttmp, tttt(:,11), barrier="on" )

          call rsdft_allreduce_sum( ab_blk, comm_grid )

          ab(i0:i1,j0:j1) = ab_blk(:,:)

          !call watchb( ttmp, tttt(:,12), barrier="on" )

          tr_blk = transpose( ab_blk )
          ab(j0:j1,i0:i1) = tr_blk(:,:)

          !call watchb( ttmp, tttt(:,13), barrier="on" )

       end if

    end do ! iblk

    !call watchb( ttmp, barrier="on" )

    if ( nprocs_b > 1 ) then
! ---(1)
!       call rsdft_allreduce_sum( ab, comm_band )
! ---(2)
       if ( .not.allocated(ab_tmp) ) then
          allocate( ab_tmp(nb,nb) ); ab_tmp=0.0d0
       end if
       ab_tmp=ab
       call MPI_Allreduce( ab_tmp, ab, size(ab), MPI_REAL8, MPI_SUM, comm_band, ierr )
! ------
    end if

    !call watchb( ttmp, tttt(:,14), barrier="on" )

    if ( nprocs_b > 1 ) call matrix_permutation( ab, nblk )

    !call watchb( ttmp, tttt(:,15), barrier="on" )

!    if ( myrank == 0 ) then
!        do i=1,15
!           write(*,'(1x,"time_calc_overlap_bp(",i2.2,")",2f10.5)') i,tttt(:,i)
!        end do
!    end if

    !call write_border( 1, "calc_overlap_bp(end)" )

  end subroutine calc_overlap_bp


  subroutine matrix_permutation( a, nblk )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,intent(in) :: nblk
    integer :: i,j,nb, blk_size, i0,i1,j0,j1,ib,irank_b
    real(8),save,allocatable :: b(:,:)
    real(8) :: ttmp(2),tttt(2,2)

    !call write_border( 1, "matrix_permutation(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    nb = size( a, 1 )
    blk_size = nb/(nblk*nprocs_b)

    if ( .not.allocated(b) ) then
       allocate( b(nb,nb) ); b=0.0d0
    end if

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    j0 = 1
    j1 = blk_size
    do irank_b=0,nprocs_b-1
       do ib=1,nb,nb/nblk
          i0 = ib + irank_b*blk_size
          i1 = i0 + blk_size - 1
          do i=1,blk_size
          do j=1,nb
             b(j,j0+i-1)=a(j,i0+i-1)
          end do
          end do
          j0 = j0 + blk_size
          j1 = j1 + blk_size
       end do
    end do ! irank
    j0 = 1
    j1 = blk_size
    do irank_b=0,nprocs_b-1
       do ib=1,nb,nb/nblk
          i0 = ib + irank_b*blk_size
          i1 = i0 + blk_size - 1
          do j=1,nb
          do i=1,blk_size
             a(j0+i-1,j)=b(i0+i-1,j)
          end do
          end do
          j0 = j0 + blk_size
          j1 = j1 + blk_size
       end do
    end do ! irank

    !call watchb( ttmp, tttt(:,2), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,2
!          write(*,'(4x,"time_matrix_permutation(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "matrix_permutation(end)" )

  end subroutine matrix_permutation

! ---------- dummy routines for different argument type
!
  subroutine calc_overlap_bp_zd( nb, a, b, alpha, ab )
    implicit none
    integer,intent(in)  :: nb
    complex(8),intent(in)  :: a(:,:)
    real(8),intent(in)  :: b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    call stop_program( "calc_overlap_bp_zd is not implemented" )
  end subroutine calc_overlap_bp_zd
  subroutine calc_overlap_bp_zz( nb, a, b, alpha, ab )
    implicit none
    integer,intent(in)  :: nb
    complex(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    call stop_program( "calc_overlap_bp_zz is not implemented" )
  end subroutine calc_overlap_bp_zz
!
! --------------------------------------------------------------

end module calc_overlap_bp_module
