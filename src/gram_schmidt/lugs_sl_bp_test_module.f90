MODULE lugs_sl_bp_test_module

  use io_tools_module
  use sl_tools_module
  use vector_tools_module
  use watch_module, only: timer => watchb, write_timer => write_watchb
  use calc_overlap_module
  use overlap_bp_module
  use reconst_module, only: reconstruct_allgatherv
  use dtrmm_bp_module

  implicit none

  PRIVATE
  PUBLIC :: LUGS_sl_bp_test

  integer :: myrank, nprocs
  integer :: nprow, npcol
  integer :: mbsize, nbsize
  integer :: icontxt, icontxt2
  integer :: icontxt_sys
  integer :: myrow, mycol
  integer :: desca(9), descz(9)
  integer :: descb(9)
  integer :: LDR,LDC
  integer,allocatable :: usermap_backup(:,:)
  character(1) :: UPLO='U'

  logical :: flag_init  = .false.
  logical :: flag_init2 = .false.

  INTERFACE LUGS_sl_bp_test
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

    call write_border( 1, " init_sl_bp_test(start)" )

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

    allocate( usermap_backup(2,0:nprocs-1) ) ; usermap_backup=0
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
       if ( k == myrank ) usermap_backup(1:2,k)=(/i,j/)
    end do
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, usermap_backup, 2*nprocs &
         , MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

    call blacs_gridmap( icontxt, usermap, nprow, nprow, npcol )
    deallocate( usermap )

    call blacs_gridinfo( icontxt, nprow, npcol, myrow, mycol )

    LDR = NUMROC( n, mbsize, myrow, 0, nprow )
    LDC = NUMROC( n, nbsize, mycol, 0, npcol )

    call descinit( desca, n,n,mbsize,nbsize,0,0,icontxt,LDR,ierr )
    call descinit( descz, n,n,mbsize,nbsize,0,0,icontxt,LDR,ierr )

    flag_init = .true.

    call write_border( 1, " init_sl_bp_test(end)" )

  END SUBROUTINE init_sl_test


  SUBROUTINE init_sl_test2( n, v )
    implicit none
    integer,intent(IN) :: n
    type(vinfo) :: v(2)
    integer,allocatable :: usermap(:,:)
    integer :: i,j,k,np_grid,np_band
    integer :: mgs,mbs,ierr,myrank
    include 'mpif.h'
    integer :: npr,npc,myr,myc
    integer :: NUMROC

    if ( flag_init2 ) return

    call write_border( 1, " init_sl_bp_test2(start)" )

    call mpi_comm_rank( mpi_comm_world, myrank, ierr )

    np_grid = v(1)%pinfo%np
    np_band = v(2)%pinfo%np

    allocate( usermap(np_grid,np_band) ) ; usermap=0

    k=-1
    do j=1,np_band
    do i=1,np_grid
       k=k+1
       usermap(i,j)=k
    end do
    end do

    icontxt2 = icontxt_sys
    call blacs_gridmap( icontxt2, usermap, np_grid, np_grid, np_band )
    call blacs_gridinfo( icontxt2, npr,npc,myr,myc )
    !write(*,*) npr,npc,myr,myc,myrank
    !call stop_program("")

    mgs = n/npr
    mbs = n/npc

!    write(*,*) mgs,mbs,NUMROC(n,mgs,myr,0,npr),NUMROC(n,mbs,myc,0,npc)
!    call stop_program("")

    call descinit( descb, n,n,mgs,mbs,0,0,icontxt2,mgs,ierr )

    deallocate( usermap )

    flag_init2 = .true.

    call write_border( 1, " init_sl_bp_test2(end)" )

  END SUBROUTINE init_sl_test2


  SUBROUTINE d_LUGS_sl_test( u, v )
    implicit none
    real(8),intent(INOUT)  :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    real(8),allocatable :: S(:,:),Ssub(:,:),utmp(:,:)
    real(8),allocatable :: T1(:,:),T2(:,:),T3(:,:),TT(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8) :: check(3,3)
    integer :: m,n,ierr,i,j,k,it,n0,n1
    integer,allocatable :: ipiv(:)
    type(slinfo) :: sl
    include 'mpif.h'

    call write_border( 1, " d_LUGS_sl_bp_test(start)" )

    it=0
    time(:,:)=0.0d0
    call timer( ttmp )

    m = size( u, 1 )
    n = sum( v(2)%pinfo%ir )
    n0 = v(2)%pinfo%id(v(2)%pinfo%me) + 1
    n1 = n0 + size(u,2) - 1

    call init_sl_test( n )
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="ini1"
    call init_sl_test2( n, v )
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="ini2"

    sl%myrow     = myrow
    sl%mycol     = mycol
    sl%nprow     = nprow
    sl%npcol     = npcol
    sl%mbsize    = mbsize
    sl%nbsize    = nbsize
    sl%icontxt_a = icontxt
    sl%icontxt_b = icontxt2
    sl%desca(:)  = desca(:)
    sl%descb(:)  = descb(:)

    allocate( sl%map_1to2(2,0:nprocs-1) ) ; sl%map_1to2=0
    sl%map_1to2(:,:)=usermap_backup(:,:)

    allocate( sl%map_2to1(0:sl%nprow-1,0:sl%npcol-1) ) ; sl%map_2to1=0
    do k=0,sl%nprow*sl%npcol-1
       i=sl%map_1to2(1,k)
       j=sl%map_1to2(2,k)
       sl%map_2to1(i,j)=k
    end do

    call sl_block_map( n, n, sl )

!    allocate( S(n,n)        ) ; S=zero     !(a1)

    allocate( Ssub(LDR,LDC) ) ; Ssub=zero

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="ini3"

!    call dsyrk( UPLO, 'C', n, m, one, u, m, zero, S, n )
!    call calc_overlap( u, u, 1.0d0, S )
!    call mpi_allreduce( mpi_in_place, S, size(S), mpi_real8, &
!                        mpi_sum, v(1)%pinfo%comm, ierr )

    call calc_overlap_bp( Ssub, u, v, sl, UPLO )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="ovlp"

!    call distribute_matrix( sl, S, Ssub )
!    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dist"

! ---
!    allocate( T1(n,n) ) ; T1=0.0d0
!    allocate( T2(n,n) ) ; T2=0.0d0
!    allocate( T3(n,n) ) ; T3=0.0d0
!    call gather_matrix( sl, Ssub, T1, icontxt )
! ---

    call pdpotrf( UPLO, n, Ssub, 1, 1, desca, ierr ) ! Cholesky factorization

!    call gather_matrix( sl, Ssub, T2, icontxt )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="pdpotrf"

    call pdtrtri( UPLO, 'N', n, Ssub, 1, 1, desca, ierr )

!    call gather_matrix( sl, Ssub, T3, icontxt )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="pdtrtri"

! ---
!    call clear_lu( "U", T2 )
!    call clear_lu( "U", T3 )
!    call check_lu( check(:,1), T1 )
!    call check_lu( check(:,2), T2 )
!    call check_lu( check(:,3), T3 )
!    write(*,*) "S",check(:,1)
!    write(*,*) "C",check(:,2)
!    write(*,*) "B",check(:,3)
!    allocate( TT(n,n) ) ; TT=0.0d0
!    TT = matmul( T2, T3 )
!    write(*,*) sum(TT)
!    TT = transpose( T2 )
!    T3 = matmul( TT, T2 )
!    write(*,*) sum((T1-T3)**2)
! ---

#ifdef _TEST_
    allocate( S(n,n)        ) ; S=zero     !(a1)

    call gather_matrix( sl, Ssub, S, icontxt )

    !deallocate( Ssub )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="gather"

#ifdef _TEST_OLD_

    allocate( utmp(m,size(u,2)) ) ; utmp=0.0d0
    utmp=u
    call reconstruct_allgatherv( utmp, v )
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="gathe2"

    if ( UPLO == 'U' .or. UPLO == 'u' ) then
       call dtrmm( 'R', 'U', 'N', 'N', m, n, 1.0d0, S, n, utmp, m )
    else
       call dtrmm( 'R', 'L', 'T', 'N', m, n, 1.0d0, S, n, utmp, m )
    end if
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dtrmm"

    u(:,1:n1-n0+1) = utmp(:,n0:n1)

    deallocate( utmp )

#else

    !call dtrmm_bp( UPLO, u, v, S )
    !it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dtrmmbp"
    call dtrmm_bp( UPLO, sl, u, v, Ssub, S )
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="adtrmmbp"

#endif

    deallocate( S )

# else

    allocate( S(1,1) )

    call dtrmm_bp( UPLO, sl, u, v, Ssub, S )
    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="bdtrmmbp"

    deallocate( S )

#endif

! ---

    deallocate( sl%map_1to2 )
    deallocate( sl%map_2to1 )
    deallocate( Ssub )

!    call blacs_gridexit( icontxt )

    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="finalize"

!    call write_timer( time(:,1:it), it, tlabel(1:it) )

    call write_border( 1, " d_LUGS_sl_bp_test(end)" )

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


  SUBROUTINE check_lu( check, T )
    implicit none
    real(8),intent(OUT) :: check(3)
    real(8),intent(IN)  :: T(:,:)
    integer :: i,j
    do j=1,size(T,2)
    do i=1,size(T,1)
       if ( i <  j ) check(1)=check(1)+abs( T(i,j) )
       if ( j <  i ) check(2)=check(2)+abs( T(i,j) )
       if ( i == j ) check(3)=check(3)+abs( T(i,j) )
    end do
    end do
  END SUBROUTINE check_lu


  SUBROUTINE clear_lu( indx, T )
    implicit none
    character(1),intent(IN) :: indx
    real(8),intent(INOUT)  :: T(:,:)
    integer :: i,j
    select case( indx )
    case( "u", "U" )
       do j=1,size(T,2)
       do i=1,size(T,1)
          if ( j < i ) T(i,j)=0.0d0
       end do
       end do
    case( "l", "L" )
       do j=1,size(T,2)
       do i=1,size(T,1)
          if ( i < j ) T(i,j)=0.0d0
       end do
       end do
    end select
  END SUBROUTINE clear_lu


END MODULE lugs_sl_bp_test_module
