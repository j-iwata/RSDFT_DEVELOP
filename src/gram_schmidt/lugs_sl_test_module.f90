MODULE lugs_sl_test_module

  use io_tools_module
  use sl_tools_module
  use vector_tools_module
  use watch_module, only: timer => watchb, write_timer => write_watchb
  use calc_overlap_module

  implicit none

  PRIVATE
  PUBLIC :: LUGS_sl_test

  integer :: myrank, nprocs
  integer :: nprow, npcol
  integer :: mbsize, nbsize
  integer :: icontxt
  integer :: icontxt_sys
  integer :: myrow, mycol
  integer :: desca(9), descz(9)
  integer :: LDR,LDC
  character(1) :: UPLO='U'

  logical :: flag_init = .false.

  INTERFACE LUGS_sl_test
     MODULE PROCEDURE d_LUGS_sl_test, z_LUGS_sl_test
  END INTERFACE

  real(8) :: time(2,10),ttmp(2)
  character(10) :: tlabel(10)

CONTAINS

  SUBROUTINE init_sl_test( n )

    implicit none
    integer,intent(IN) :: n
    integer :: i,j,k,l,ierr
    integer,allocatable :: usermap(:,:)
    integer,external :: NUMROC
    include 'mpif.h'

    if ( flag_init ) return

    call write_border( 1, " init_sl_test(start)" )

    call IOTools_readIntegerKeyword( "NPROW" , nprow )
    call IOTools_readIntegerKeyword( "NPCOL" , npcol )
    call IOTools_readIntegerKeyword( "MBSIZE", mbsize )
    call IOTools_readIntegerKeyword( "NBSIZE", nbsize )
    call IOTools_readStringKeyword( "UPLO", UPLO )

    call blacs_get( 0, 0, icontxt_sys )

    icontxt = icontxt_sys
!    call blacs_gridinit( icontxt, 'R', nprow, npcol )

    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

!#ifdef TEST
    allocate( usermap(0:nprow-1,0:npcol-1) ) ; usermap=0
    i=nprocs/(nprow*npcol)
    do j=0,i-1
       l=j*nprow*npcol
       if ( myrank >= l ) k=l
    end do
    k=k-1
    do i=0,nprow-1
    do j=0,npcol-1
       k=k+1   
       usermap(i,j) = k
    end do
    end do
    call blacs_gridmap( icontxt, usermap, nprow, nprow, npcol )
    deallocate( usermap )
!#endif

    call blacs_gridinfo( icontxt, nprow, npcol, myrow, mycol )

    LDR = NUMROC( n, mbsize, myrow, 0, nprow )
    LDC = NUMROC( n, nbsize, mycol, 0, npcol )

    call descinit( desca, n,n,mbsize,nbsize,0,0,icontxt,LDR,ierr )
    call descinit( descz, n,n,mbsize,nbsize,0,0,icontxt,LDR,ierr )

    flag_init = .true.

    call write_border( 1, " init_sl_test(end)" )

  END SUBROUTINE init_sl_test


  SUBROUTINE d_LUGS_sl_test( u, v )
    implicit none
    real(8),intent(INOUT)  :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    real(8),allocatable :: S(:,:),Ssub(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0
    real(8) :: dV
    integer :: m,n,ierr,i,j,it
    integer,allocatable :: ipiv(:)
    type(slinfo) :: sl
    include 'mpif.h'

    call write_border( 1, " d_LUGS_sl_test(start)" )

    it=0
    time(:,:)=0.0d0
    call timer( ttmp )

    m  = size( u, 1 )
    n  = sum( v(2)%pinfo%ir )
    dV = v(1)%factor

    call init_sl_test( n )

    sl%myrow  = myrow
    sl%mycol  = mycol
    sl%nprow  = nprow
    sl%npcol  = npcol
    sl%mbsize = mbsize
    sl%nbsize = nbsize

    allocate( S(n,n)        ) ; S=zero
    allocate( Ssub(LDR,LDC) ) ; Ssub=zero

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="init"

    call dsyrk( UPLO, 'C', n, m, dV, u, m, zero, S, n )
!    call calc_overlap( u, u, dV, S )
    call mpi_allreduce( mpi_in_place, S, size(S), mpi_real8, &
         mpi_sum, v(1)%pinfo%comm, ierr )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dsyrnk"

    call distribute_matrix( sl, S, Ssub )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dist"

    call pdpotrf( UPLO, n, Ssub, 1, 1, desca, ierr ) ! Cholesky factorization

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="pdpotrf"

    call pdtrtri( UPLO, 'N', n, Ssub, 1, 1, desca, ierr )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="pdtrtri"

    call gather_matrix( sl, Ssub, S, icontxt )

    deallocate( Ssub )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="gather"

!    allocate( utmp(m,n) ) ; utmp=zero
!    call dgemm( 'N', 'N', m, n, n, one, u, m, S, n, zero, utmp, m )
!    u = utmp
!    deallocate( utmp )
!    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dgemm"

    if ( UPLO == 'U' .or. UPLO == 'u' ) then
       call dtrmm( 'R', 'U', 'N', 'N', m, n, 1.0d0, S, n, u, m )
    else
       call dtrmm( 'R', 'L', 'T', 'N', m, n, 1.0d0, S, n, u, m )
    end if
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dtrmm"

    deallocate( S )

!    call blacs_gridexit( icontxt )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="finalize"

!    call write_timer( time(:,1:it), it, tlabel(1:it) )

    call write_border( 1, " d_LUGS_sl_test(end)" )

  END SUBROUTINE d_LUGS_sl_test


  SUBROUTINE z_LUGS_sl_test( u, v )
    implicit none
    complex(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN)   :: v(2)
    complex(8),allocatable :: S(:,:),Ssub(:,:),Bsub(:,:),utmp(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    integer :: m,n,ierr,i,j
    integer,allocatable :: ipiv(:)
    type(slinfo) :: sl
    include 'mpif.h'

    call write_border( 1, " z_LUGS_sl_test(start)" )

    m = size( u, 1 )
    n = sum( v(2)%pinfo%ir )

    call init_sl_test( n )

    sl%myrow  = myrow
    sl%mycol  = mycol
    sl%nprow  = nprow
    sl%npcol  = npcol
    sl%mbsize = mbsize
    sl%nbsize = nbsize

    allocate( S(n,n)        ) ; S=zero
    allocate( Ssub(LDR,LDC) ) ; Ssub=zero
    allocate( Bsub(LDR,LDC) ) ; Bsub=zero

    do i=1,n
       S(i,i) = one
    end do
    call distribute_matrix( sl, S, Bsub )

    S = zero
    call zherk( 'U', 'C', n, m, one, u, m, zero, S, n )
    call mpi_allreduce( mpi_in_place, S, size(S), mpi_complex16, &
         mpi_sum, v(1)%pinfo%comm, ierr )

    call fill_matrix( S )
    call distribute_matrix( sl, S, Ssub )

    call pzpotrf( 'U', n, Ssub, 1, 1, desca, ierr )

    call gather_matrix( sl, Ssub, S )
    do j=1,n
       do i=j+1,n
          S(i,j)=zero
       end do
    end do
    call distribute_matrix( sl, S, Ssub )

    allocate( ipiv(n) ) ; ipiv=0
    do i=1,n
       ipiv(i)=i
    end do
!
! If pivot information is necessary, following routine should be called
!    call pzgetrf( n, n, Ssub, 1, 1, desca, ipiv, ierr )

    call pzgetrs( 'N', n, n, Ssub, 1, 1, desca, ipiv, Bsub, 1, 1, descz, ierr )

    deallocate( ipiv )

    call gather_matrix( sl, Bsub, S )

    deallocate( Bsub )
    deallocate( Ssub )

    allocate( utmp(m,n) ) ; utmp=zero
    call zgemm( 'N', 'N', m, n, n, one, u, m, S, n, zero, utmp, m )
    u = utmp
    deallocate( utmp )

    deallocate( S )

    call blacs_gridexit( icontxt )

    call write_border( 1, " z_LUGS_sl_test(end)" )

  END SUBROUTINE z_LUGS_sl_test


END MODULE lugs_sl_test_module
