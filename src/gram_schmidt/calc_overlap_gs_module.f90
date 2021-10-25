module calc_overlap_gs_module

  implicit none

  private
  public :: calc_overlap_gs

  ! real(8),allocatable :: sendbuf(:,:), recvbuf(:,:)
  ! real(8),allocatable :: ab_blk(:,:), tr_blk(:,:)
  ! real(8),allocatable :: ab_tmp(:,:)

  interface calc_overlap_gs
    module procedure d_calc_overlap_gs, z_calc_overlap_gs
  end interface

contains

  ! mdsource/calc_overlap_bp2_module.f90 がベース
  ! ただし入力が１つ（ a(ng,nb））である事と、上三角行列を作る用である事が違う

  subroutine d_calc_overlap_gs( nb, a, alpha, aa )
    use parallel_module, only: comm_band,np_band,myrank_b,ir_band,id_band
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    integer,intent(in)  :: nb !フルのバンドサイズ
    real(8),intent(in)  :: a(:,:)
    real(8),intent(in)  :: alpha
    real(8),intent(inout) :: aa(:,:)
    include 'mpif.h'
    integer :: ng,nstep,mb0,mb1,nb0,nb1,istep,i,j,n
    integer :: nblk0,nblk0_size,nblk1,nblk1_size,iblk1
    integer :: ib0,ib1,jb0,jb1,i0,j0,nib,njb,mib,nbh,mbh
    integer :: irank_send,irank_recv,irank
    integer,allocatable :: ndata_recv(:), ndata_send(:)
    integer :: istatus(MPI_STATUS_SIZE,2), ireq(2), ierr
    real(8),allocatable :: f_send(:,:), f_recv(:,:)
    real(8),parameter :: z0=0.0d0, z1=1.0d0
    real(8) :: tmp1,tmp2
    ! real(8) :: ttmp(2),tt(2,0:15)

    call write_border( 1, ' d_calc_overlap_gs(start)' )
    ! tt=0.0d0; call watchb( tt(:,0), barrier='on' )

    aa(:,:)=0.0d0

    ng = size( a, 1 ) !(並列分割された)グリッドサイズ

    if ( np_band == 1 ) then !バンド並列を行わない場合
      call DSYRK( 'U', 'C', nb, ng, alpha, a, ng, z0, aa, nb )
      call rsdft_allreduce( aa, 'g' )
      call write_border( 1, ' d_calc_overlap_gs(return)' )
      return
    end if

    nstep = np_band/2 !上隣にsendする回数

    mb0 = id_band(myrank_b) + 1
    mb1 = mb0 + ir_band(myrank_b) - 1
    
    nb0 = mb0
    nb1 = mb1

    nblk0_size = mb1 - mb0 + 1
    nblk0 = ( nb + nblk0_size - 1 )/nblk0_size 

    nblk1 = 1 !バンド並列担当分をさらにnblk1に分割
    nblk1_size = ( nblk0_size + nblk1 - 1 )/nblk1

    n = maxval( ir_band )
    allocate( f_send(ng,n) ); f_send=z0
    allocate( f_recv(ng,n) ); f_recv=z0

    irank_send = mod( myrank_b-1+np_band, np_band )
    irank_recv = mod( myrank_b+1+np_band, np_band )

    allocate( ndata_send(nstep) ); ndata_send=0
    allocate( ndata_recv(nstep) ); ndata_recv=0

    !np_bandが偶数の時は最後のstepは半分送れば良い（が、サボってそういう実装はしていない）
    do istep = 1, nstep
      irank = mod( myrank_b+istep-1+np_band, np_band )
      ndata_send(istep) = ir_band(irank)
      irank = mod( myrank_b+1+istep-1+np_band, np_band )
      ndata_recv(istep) = ir_band(irank)
    end do
 
    n = ndata_send(1)
    if ( n > 0 ) f_send(:,1:n) = a(:,1:n)

    do istep = 1, nstep

      call MPI_Irecv(f_recv,ng*ndata_recv(istep),MPI_REAL8,irank_recv,0,comm_band,ireq(1),ierr)
      call MPI_Isend(f_send,ng*ndata_send(istep),MPI_REAL8,irank_send,0,comm_band,ireq(2),ierr)

      if ( istep == 1 ) then

        !対角ブロック。DGEMMとDSYRKのどちらの選択が良いかは性能依存。
        do iblk1 = 1, nblk1
          ib0 = mb0 + (iblk1-1)*nblk1_size
          ib1 = min( ib0 + nblk1_size - 1, mb1 )
          jb0 = ib0
          jb1 = ib1
          nib = ib1 - ib0 + 1
          njb = nib
          i0 = ib0 - mb0 + 1
          j0 = jb0 - nb0 + 1
          ! call DGEMM('T','N',nib,njb,ng,alpha,a(1,i0),ng,a(1,j0),ng,z0,aa(ib0,jb0),nb)
          call DSYRK( 'U', 'C', nib, ng, alpha, a(1,i0), ng, z0, aa(ib0,jb0), nb )
        end do !iblk1

      else ![ istep > 1 ]

        nb0 = mod( nb1 + nb, nb ) + 1
        nb1 = mod( nb0 + ndata_recv(istep-1) - 2 + nb, nb ) + 1
        mib = mb1 - mb0 + 1
        njb = nb1 - nb0 + 1
        i0 = 1
        j0 = 1
        if ( nb0 >= mb0 ) then
          call DGEMM('T','N',mib,njb,ng,alpha,a(1,i0),ng,f_send(1,j0),ng,z1,aa(mb0,nb0),nb)
        else
          call DGEMM('T','N',njb,mib,ng,alpha,f_send(1,j0),ng,a(1,i0),ng,z1,aa(nb0,mb0),nb)
        end if
      
      end if ![ istep == 1 ]

      call MPI_Waitall( 2, ireq, istatus, ierr )

      if ( istep < nstep ) then
        f_send(:,1:ndata_recv(istep)) = f_recv(:,1:ndata_recv(istep))
      end if

    end do !istep

    if ( mod(np_band,2) == 0 ) then

      nb0 = mod( nb1 + nb, nb ) + 1
      nb1 = mod( nb0 + ndata_recv(istep-1) - 2 + nb, nb ) + 1
      mib = mb1 - mb0 + 1
      njb = nb1 - nb0 + 1
      i0 = 1
      j0 = 1
      if ( nb0 >= mb0 ) then
        !
        ! call DGEMM('T','N',mib,njb,ng,alpha,a(1,i0),ng,f_recv(1,j0),ng,z1,aa(mb0,nb0),nb)
        !
        nbh = nb0 + (njb+1)/2 - 1
        njb = nbh - nb0 + 1
        call DGEMM('T','N',mib,njb,ng,alpha,a(1,i0),ng,f_recv(1,j0),ng,z1,aa(mb0,nb0),nb)
        !
      else
        !
        ! call DGEMM('T','N',njb,mib,ng,alpha,f_recv(1,j0),ng,a(1,i0),ng,z1,aa(nb0,mb0),nb)
        !
        i0 = i0 + (mib+1)/2
        mbh = mb0 + (mib+1)/2
        mib = mb1 - mbh + 1
        call DGEMM('T','N',njb,mib,ng,alpha,f_recv(1,j0),ng,a(1,i0),ng,z1,aa(nb0,mbh),nb)
        !
      end if

    else

      nb0 = mod( nb1 + nb, nb ) + 1
      nb1 = mod( nb0 + ndata_recv(istep-1) - 2 + nb, nb ) + 1
      mib = mb1 - mb0 + 1
      njb = nb1 - nb0 + 1
      i0 = 1
      j0 = 1
      if ( nb0 >= mb0 ) then
        call DGEMM('T','N',mib,njb,ng,alpha,a(1,i0),ng,f_recv(1,j0),ng,z1,aa(mb0,nb0),nb)
      else
        call DGEMM('T','N',njb,mib,ng,alpha,f_recv(1,j0),ng,a(1,i0),ng,z1,aa(nb0,mb0),nb)
      end if

    end if

    deallocate( f_recv, f_send )

    call rsdft_allreduce( aa, 'b' )
    call rsdft_allreduce( aa, 'g' )

    ! do j = 1, nb
    !   do i = 1, j-1
    !     tmp1 = aa(i,j)
    !     if ( tmp1 == 0.0d0 ) then
    !       tmp2 = aa(j,i)
    !       if ( tmp2 /= 0.0d0 ) aa(i,j) = tmp2
    !     end if
    !     aa(j,i) = 0.0d0 !この処理は必要ない（キレイな上三角にするためだけ）
    !   end do !i
    ! end do !j

    ! call watchb( tt(:,1), barrier='on' )
    ! if ( myrank == 0 ) write(*,'(1x,"time:",2f10.3)') (tt(i,1)-tt(i,0),i=1,2)

    call write_border( 1, ' d_calc_overlap_gs(end)' )

  end subroutine d_calc_overlap_gs


  subroutine z_calc_overlap_gs( nb, a, alpha, aa )
    integer,intent(in) :: nb !フルのバンドサイズ
    complex(8),intent(in) :: a(:,:)
    real(8),intent(in) :: alpha
    real(8),intent(inout) :: aa(:,:)
    call stop_program('z_calc_overlap_gs is not implemented yet')
  end subroutine z_calc_overlap_gs


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

end module calc_overlap_gs_module
