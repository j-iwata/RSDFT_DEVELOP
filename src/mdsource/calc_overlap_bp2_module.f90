module calc_overlap_bp2_module

  ! use parallel_module, only: nprocs_b, myrank_b, comm_band, myrank, comm_grid
  ! use rsdft_mpi_module
  use watch_module, only: watchb
  ! use transpose_module, only: rsdft_transpose

  implicit none

  private
  public :: calc_overlap_bp2

  real(8),allocatable :: sendbuf(:,:), recvbuf(:,:)
  real(8),allocatable :: ab_blk(:,:), tr_blk(:,:)
  real(8),allocatable :: ab_tmp(:,:)


contains

  !DGEMMを使ってないので要修正
  subroutine calc_overlap_bp2( nb, a, b, alpha, ab )
    use parallel_module, only: comm_band, np_band, myrank_b, myrank
    use cpmd_variables, only: id_band_cpmd, ir_band_cpmd
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    integer,intent(in)  :: nb !フルのバンドサイズ
    real(8),intent(in)  :: a(:,:), b(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: ab(:,:)
    include 'mpif.h'
    integer :: ng
    integer :: nblk0, nblk0_size, nblk1, nblk1_size, iblk1
    integer :: mb0,mb1,nb0,nb1,ib0,ib1,jb0,jb1
    integer :: m,n,nblk,blk_size
    integer :: ib,i0,i1,j0,j1,i,j,iblk,k0,k1,b0,b1,l0,l1
    integer :: irank_send, irank_recv, istep, nstep, irank
    integer,allocatable :: ndata_recv(:), ndata_send(:)
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr
    real(8) :: ttmp(2),tttt(2,0:15)
    real(8),allocatable :: f_send(:,:), f_recv(:,:)
    real(8) :: tmp1,tmp2

    ! call write_border( 1, ' calc_overlap_bp2(start)' )
    ! tttt=0.0d0; call watchb( tttt(:,0), barrier='on' )

    ab(:,:)=0.0d0

    nstep = np_band/2 !左隣にsendする回数

    nb0 = id_band_cpmd(myrank_b) + 1
    nb1 = nb0 + ir_band_cpmd(myrank_b) - 1
    
    mb0 = nb0
    mb1 = nb1

    nblk0_size = nb1 - nb0 + 1
    nblk0 = ( nb + nblk0_size - 1 )/nblk0_size 

    nblk1 = 1 !バンド並列担当分をさらにnblk1に分割
    nblk1_size = ( nblk0_size + nblk1 - 1 )/nblk1

    if ( nstep == 0 ) then
      do iblk1 = 1, nblk1
        jb0 = nb0 + (iblk1-1)*nblk1_size
        jb1 = min( jb0 + nblk1_size - 1, nb1 )
        ib0 = jb0
        ib1 = jb1
        do j = jb0, jb1
          j0 = j - nb0 + 1
          do i = j, ib1
            i0 = i - nb0 + 1
            ab(i,j) = sum( a(:,i0)*b(:,j0) )*alpha
            ab(j,i) = ab(i,j)
          end do
        end do
      end do !iblk1
      call rsdft_allreduce( ab, 'g' )
      ! call watchb( tttt(:,1), barrier='on' )
      ! if ( myrank == 0 ) write(*,'(1x,"time:",2f10.3)') (tttt(i,1)-tttt(i,0),i=1,2)
      ! call write_border( 0, ' calc_overlap_bp2(return)' )
      return
    end if

    ng = size( a, 1 ) !(並列分割された)グリッドサイズ

    n = maxval( ir_band_cpmd )
    allocate( f_send(ng,n) ); f_send=0.0d0
    allocate( f_recv(ng,n) ); f_recv=0.0d0

    irank_send = mod( myrank_b-1+np_band, np_band )
    irank_recv = mod( myrank_b+1+np_band, np_band )

    allocate( ndata_send(nstep) ); ndata_send=0
    allocate( ndata_recv(nstep) ); ndata_recv=0

    !np_bandが偶数の時は最後のstepは半分送れば良い（が、サボってそういう実装はしていない）
    do istep = 1, nstep
      irank = mod( myrank_b+istep-1+np_band, np_band )
      ndata_send(istep) = ir_band_cpmd(irank)
      irank = mod( myrank_b+1+istep-1+np_band, np_band )
      ndata_recv(istep) = ir_band_cpmd(irank)
    end do
    
    n = ndata_send(1)
    if ( n > 0 ) f_send(:,1:n) = a(:,1:n)

    do istep = 1, nstep

      call MPI_Irecv(f_recv,ng*ndata_recv(istep),MPI_REAL8,irank_recv,0,comm_band,ireq(1),ierr)
      call MPI_Isend(f_send,ng*ndata_send(istep),MPI_REAL8,irank_send,0,comm_band,ireq(2),ierr)

      if ( istep == 1 ) then
        
        do iblk1 = 1, nblk1
          jb0 = nb0 + (iblk1-1)*nblk1_size
          jb1 = min( jb0 + nblk1_size - 1, nb1 )
          ib0 = jb0
          ib1 = jb1
          do j = jb0, jb1
            j0 = j - nb0 + 1
            do i = j, ib1
              i0 = i - nb0 + 1
              ab(i,j) = sum( a(:,i0)*b(:,j0) )*alpha
            end do
          end do
        end do !iblk1

      else ![ istep > 1 ]

        mb0 = mod( mb1 + nb, nb ) + 1
        mb1 = mod( mb0 + ndata_recv(istep-1) - 2 + nb, nb ) + 1

        do j = nb0, nb1
          j0 = j - nb0 + 1
          do i = mb0, mb1
            i0 = i - mb0 + 1
            ab(i,j) = sum( f_recv(:,i0)*b(:,j0) )*alpha
          end do
        end do

      end if ![ istep == 1 ]

      call MPI_Waitall( 2, ireq, istatus, ierr )

      if ( istep < nstep ) then
        f_send(:,1:ndata_recv(istep)) = f_recv(:,1:ndata_recv(istep))
      end if

    end do !istep

    mb0 = mod( mb1 + nb, nb ) + 1
    mb1 = mod( mb0 + ndata_recv(istep-1) - 2 + nb, nb ) + 1
    
    do j = nb0, nb1
      j0 = j - nb0 + 1
      do i = mb0, mb1
        i0 = i - mb0 + 1
        ab(i,j) = sum( f_recv(:,i0)*b(:,j0) )*alpha
      end do
    end do

    deallocate( f_recv, f_send )

    call rsdft_allreduce( ab, 'b' )
    call rsdft_allreduce( ab, 'g' )
    do j = 1, nb
      do i = j, nb
        tmp1 = ab(i,j)
        tmp2 = ab(j,i)
        if ( tmp1 == 0.0d0 ) then
          if ( tmp2 == 0.0d0 ) call stop_program('calc_overlap_bp2')
          ab(i,j) = tmp2
        else
          ab(j,i) = tmp1
        end if
      end do !i
    end do !j

    ! call watchb( tttt(:,1), barrier='on' )
    ! if ( myrank == 0 ) write(*,'(1x,"time:",2f10.3)') (tttt(i,1)-tttt(i,0),i=1,2)
    ! call write_border( 1, ' calc_overlap_bp2(end)' )

  end subroutine calc_overlap_bp2


!   subroutine matrix_permutation( a, nblk )
!     implicit none
!     real(8),intent(inout) :: a(:,:)
!     integer,intent(in) :: nblk
!     integer :: i,j,nb, blk_size, i0,i1,j0,j1,ib,irank_b
!     real(8),save,allocatable :: b(:,:)
!     real(8) :: ttmp(2),tttt(2,2)

!     !call write_border( 1, "matrix_permutation(start)" )
!     !call watchb( ttmp, barrier="on" ); tttt=0.0d0

!     nb = size( a, 1 )
!     blk_size = nb/(nblk*nprocs_b)

!     if ( .not.allocated(b) ) then
!        allocate( b(nb,nb) ); b=0.0d0
!     end if

!     !call watchb( ttmp, tttt(:,1), barrier="on" )

!     j0 = 1
!     j1 = blk_size
!     do irank_b=0,nprocs_b-1
!        do ib=1,nb,nb/nblk
!           i0 = ib + irank_b*blk_size
!           i1 = i0 + blk_size - 1
!           do i=1,blk_size
!           do j=1,nb
!              b(j,j0+i-1)=a(j,i0+i-1)
!           end do
!           end do
!           j0 = j0 + blk_size
!           j1 = j1 + blk_size
!        end do
!     end do ! irank
!     j0 = 1
!     j1 = blk_size
!     do irank_b=0,nprocs_b-1
!        do ib=1,nb,nb/nblk
!           i0 = ib + irank_b*blk_size
!           i1 = i0 + blk_size - 1
!           do j=1,nb
!           do i=1,blk_size
!              a(j0+i-1,j)=b(i0+i-1,j)
!           end do
!           end do
!           j0 = j0 + blk_size
!           j1 = j1 + blk_size
!        end do
!     end do ! irank

!     !call watchb( ttmp, tttt(:,2), barrier="on" )

! !    if ( myrank == 0 ) then
! !       do i=1,2
! !          write(*,'(4x,"time_matrix_permutation(",i1,")",2f10.5)') i, tttt(:,i)
! !       end do
! !    end if

!     !call write_border( 1, "matrix_permutation(end)" )

!   end subroutine matrix_permutation

! ! ---------- dummy routines for different argument type
! !
!   subroutine calc_overlap_bp_zd( nb, a, b, alpha, ab )
!     implicit none
!     integer,intent(in)  :: nb
!     complex(8),intent(in)  :: a(:,:)
!     real(8),intent(in)  :: b(:,:)
!     real(8),intent(in)  :: alpha
!     real(8),intent(inout) :: ab(:,:)
!     call stop_program( "calc_overlap_bp_zd is not implemented" )
!   end subroutine calc_overlap_bp_zd
!   subroutine calc_overlap_bp_zz( nb, a, b, alpha, ab )
!     implicit none
!     integer,intent(in)  :: nb
!     complex(8),intent(in)  :: a(:,:), b(:,:)
!     real(8),intent(in)  :: alpha
!     real(8),intent(inout) :: ab(:,:)
!     call stop_program( "calc_overlap_bp_zz is not implemented" )
!   end subroutine calc_overlap_bp_zz
! !
! ! --------------------------------------------------------------

end module calc_overlap_bp2_module
