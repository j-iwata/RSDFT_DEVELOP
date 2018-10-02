module rotorb_sub_module

  use parallel_module, only: nprocs_b, myrank_b, comm_band, myrank
  use cpmd_variables, only: id_band_cpmd, ir_band_cpmd
  use watch_module

  implicit none

  private
  public :: rotorb_sub

contains

  subroutine rotorb_sub( a, b, c, alpha_in, beta_in )

    implicit none
    real(8),intent(in) :: a(:,:), b(:,:)
    real(8),intent(inout) :: c(:,:)
    real(8),optional,intent(in) :: alpha_in,beta_in
    include 'mpif.h'
    integer :: m,n,nb,istep,irank,jrank,itags,ireq(2),ierr
    integer :: istatus(MPI_STATUS_SIZE,2)
    integer :: m1,m2,mm,m1_0,m2_0,n1,n2,n1_0,n2_0,k1_0,k2_0
    integer :: nblk,iblk,blk_size,i
    real(8),save,allocatable :: workbuf(:,:), recvbuf(:,:)
    real(8) :: alpha, beta
    real(8) :: ttmp(2),tttt(2,9)

    !call write_border( 1, "rotorb_sub(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m  = size( a, 1 )
    n  = size( a, 2 )
    nb = size( b, 1 )

    alpha=1.0d0; if ( present(alpha_in) ) alpha=alpha_in
    beta =0.0d0; if ( present(beta_in ) ) beta =beta_in

    m1 = id_band_cpmd(myrank_b) + 1
    m2 = id_band_cpmd(myrank_b) + ir_band_cpmd(myrank_b)
    mm = m2 - m1 + 1

    nblk = 1
    blk_size = mm/nblk
    itags = 10

    if ( .not.allocated(workbuf) ) then
       allocate( workbuf(m,blk_size) ); workbuf=0.0d0
       allocate( recvbuf(m,blk_size) ); recvbuf=0.0d0
    end if

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    do iblk=1,nblk

       m1_0 = m1 + (iblk-1)*blk_size
       m2_0 = m1_0 + blk_size - 1

       k1_0 = 1 + (iblk-1)*blk_size
       k2_0 = k1_0 + blk_size - 1

       workbuf(:,:) = a(:,k1_0:k2_0)

       n1 = id_band_cpmd(myrank_b) + 1
       n2 = id_band_cpmd(myrank_b) + ir_band_cpmd(myrank_b)
       n1_0 = n1 + (iblk-1)*blk_size
       n2_0 = n1_0 + blk_size - 1

       do istep=1,nprocs_b

          !call watchb( ttmp, barrier="on" )

          if ( istep < nprocs_b ) then
             irank = mod( myrank_b + istep, nprocs_b )
             jrank = mod( nprocs_b - istep + myrank_b, nprocs_b )
             call MPI_Irecv( recvbuf  , m*blk_size, MPI_REAL8, jrank, itags, comm_band, ireq(1), ierr )
             call MPI_Isend( a(1,k1_0), m*blk_size, MPI_REAL8, irank, itags, comm_band, ireq(2), ierr )
          end if

          !call watchb( ttmp, tttt(:,2), barrier="on" )

          call DGEMM( 'N','N',m,mm,blk_size,alpha,workbuf,m,b(n1_0,m1),nb,beta,c,m)

          !call watchb( ttmp, tttt(:,3), barrier="on" )

          if ( istep == nprocs_b ) exit

          n1 = id_band_cpmd(jrank) + 1
          n2 = id_band_cpmd(jrank) + ir_band_cpmd(jrank)
          n1_0 = n1 + (iblk-1)*blk_size
          n2_0 = n1_0 + blk_size - 1

          call MPI_Waitall( 2, ireq, istatus, ierr )
          workbuf = recvbuf

          !call watchb( ttmp, tttt(:,4), barrier="on" )

       end do ! istep

    end do ! iblk

    call watchb( ttmp, barrier="on" )

    !call watchb( ttmp, tttt(:,5), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,5
!          write(*,'("time_rotorb_sub(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    call write_border( 1, "rotorb_sub(end)" )

  end subroutine rotorb_sub

end module rotorb_sub_module
