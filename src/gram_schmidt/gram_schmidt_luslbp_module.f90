module gram_schmidt_luslbp_module

  use io_tools_module
  use sl_tools_module
  use watch_module, only: timer => watchb, write_timer => write_watchb
  use overlap_bp_module, only: init_overlap_bp, calc_overlap_bp
  use dtrmm_bp_module, only: init_dtrmm_bp, dtrmm_bp
  use rgrid_variables, only: dV
  use parallel_module, only: comm_grid, comm_band

  implicit none

  private
  public :: gram_schmidt_luslbp

  integer :: nband
  integer :: myrank, nprocs
  integer :: icontxt_sys
  integer :: LDR,LDC
  character(1) :: UPLO='U'

  logical :: has_init_done  = .false.
  logical :: has_init1_done = .false.
  logical :: has_init2_done = .false.

  interface gram_schmidt_luslbp
     module procedure d_gram_schmidt_luslbp
  end interface

  real(8) :: time(2,0:9),ttmp(2)
  character(10) :: tlabel(0:9)

  type(slinfo) :: sl

  integer,external :: NUMROC

contains

  subroutine init1_lusl( n )

    implicit none
    integer,intent(in) :: n
    integer :: i,j,k,l,ierr
    integer,allocatable :: usermap(:,:)
    include 'mpif.h'
    integer,parameter :: unit=90

    if ( has_init1_done ) return

!    call write_border( 1, " init1_lusl(start)" )

    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    sl%nprow=0
    sl%npcol=0
    sl%mbsize=0
    sl%nbsize=0

!    open( unit, file='lusl.in', status='old' )
!    call IOTools_readIntegerKeyword( "NPROW" , sl%nprow, unit )
!    call IOTools_readIntegerKeyword( "NPCOL" , sl%npcol, unit )
!    call IOTools_readIntegerKeyword( "MBSIZE", sl%mbsize, unit )
!    call IOTools_readIntegerKeyword( "NBSIZE", sl%nbsize, unit )
!    call IOTools_readStringKeyword( "UPLO", UPLO, unit )
!    close( unit )

    if( sl%nprow==0 .or. sl%npcol==0 .or. sl%mbsize==0 .or. sl%nbsize==0 )then
       call restore_param_sl( sl )
    end if

    if( myrank == 0 )then
       write(*,*) "sl%nprow=",sl%nprow
       write(*,*) "sl%npcol=",sl%npcol
       write(*,*) "sl%mbsize=",sl%mbsize
       write(*,*) "sl%nbsize=",sl%nbsize
    end if

! ---

    call blacs_get( 0, 0, icontxt_sys )

    sl%icontxt_a = icontxt_sys

    allocate( sl%map_1to2(2,0:nprocs-1) ); sl%map_1to2=0

    allocate( usermap(0:sl%nprow-1,0:sl%npcol-1) ); usermap=0
    i=nprocs/(sl%nprow*sl%npcol)
    do j=0,i-1
       l=j*sl%nprow*sl%npcol
       if ( myrank >= l ) k=l
    end do
    k=k-1
    do i=0,sl%nprow-1
    do j=0,sl%npcol-1
       k=k+1   
       usermap(i,j) = k
       if ( k == myrank ) sl%map_1to2(1:2,k)=(/i,j/)
    end do
    end do
    call blacs_gridmap( sl%icontxt_a, usermap, sl%nprow, sl%nprow, sl%npcol )
    deallocate( usermap )

    call MPI_Allreduce( MPI_IN_PLACE, sl%map_1to2, size(sl%map_1to2) &
                      , MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

    allocate( sl%map_2to1(0:sl%nprow-1,0:sl%npcol-1) ); sl%map_2to1=0
    do k=0,sl%nprow*sl%npcol-1
       i=sl%map_1to2(1,k)
       j=sl%map_1to2(2,k)
       sl%map_2to1(i,j)=k
    end do

    call blacs_gridinfo( sl%icontxt_a, sl%nprow, sl%npcol, sl%myrow, sl%mycol )

    LDR = NUMROC( n, sl%mbsize, sl%myrow, 0, sl%nprow )
    LDC = NUMROC( n, sl%nbsize, sl%mycol, 0, sl%npcol )

    call descinit( sl%desca,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )
    call descinit( sl%descz,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )

    has_init1_done = .true.

!    call write_border( 1, " init1_lusl(end)" )

  end subroutine init1_lusl


  subroutine init2_lusl( n )
    implicit none
    integer,intent(in) :: n
    integer,allocatable :: usermap(:,:)
    integer :: i,j,k,mgs,mbs,ierr,npr,npc,myr,myc
    integer :: np_grid, np_band
    include 'mpif.h'

    if ( has_init2_done ) return

!    call write_border( 1, " init2_lusl(start)" )

    call MPI_Comm_size( comm_grid, np_grid, ierr )
    call MPI_Comm_size( comm_band, np_band, ierr )

    allocate( usermap(np_grid,np_band) ); usermap=0

    k=-1
    do j=1,np_band
    do i=1,np_grid
       k=k+1
       usermap(i,j)=k
    end do
    end do

    sl%icontxt_b = icontxt_sys
    call blacs_gridmap( sl%icontxt_b, usermap, np_grid, np_grid, np_band )
    call blacs_gridinfo( sl%icontxt_b, npr,npc,myr,myc )

    mgs = n/npr
    mbs = n/npc

    call descinit( sl%descb,n,n,mgs,mbs,0,0,sl%icontxt_b,mgs,ierr )

    deallocate( usermap )

    has_init2_done = .true.

!    call write_border( 1, " init2_lusl(end)" )

  end subroutine init2_lusl


  subroutine get_nband( n, nb )
    implicit none
    integer,intent(in) :: n
    integer,intent(out) :: nb
    integer :: ierr
    include 'mpif.h'
    call MPI_Allreduce(n,nb,1,MPI_INTEGER,MPI_SUM,comm_band,ierr)
  end subroutine get_nband


  subroutine d_gram_schmidt_luslbp( u )
    implicit none
    real(8),intent(inout)  :: u(:,:)
    real(8),allocatable :: Ssub(:,:)
    real(8),allocatable :: S(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    integer :: m,n,ierr,i,j,k,it,n0,n1
    integer,allocatable :: ipiv(:)
!    real(8),allocatable :: T1(:,:),T2(:,:),T3(:,:),TT(:,:)
!    real(8) :: check(3,3)

!    call write_border( 1, " d_gram_schmidt_luslbp(start)" )

!    call timer( ttmp ); time(:,:)=0.0d0; it=0

    m = size( u, 1 )
    n = size( u, 2 )

    if ( .not.has_init_done ) then
       call get_nband( n, nband )
       call init1_lusl( nband )
       call init2_lusl( nband )
       call init_overlap_bp( m, n, dV, comm_grid, comm_band )
       call init_dtrmm_bp( m, n, comm_grid, comm_band )
       has_init_done = .true.
    end if

    call sl_block_map( n, n, sl )

    allocate( S(nband,nband) ); S=zero

    allocate( Ssub(LDR,LDC) ); Ssub=zero

!    call timer( ttmp,time(:,it) ); tlabel(it)="ini3"; it=it+1

!    call dsyrk( UPLO, 'C', n, m, one, u, m, zero, S, n )
!    call calc_overlap( u, u, 1.0d0, S )
!    call mpi_allreduce( mpi_in_place, S, size(S), mpi_real8, &
!                        mpi_sum, v(1)%pinfo%comm, ierr )

    call calc_overlap_bp( Ssub, u, sl, UPLO )
!write(*,'("Ssub",3i4,3f10.5)') myrank,LDR,LDC,Ssub
!call gather_matrix( sl, Ssub, S, sl%icontxt_a )
!write(*,*) "S",sum(S**2),count(S/=zero),size(S,2),nband
!do i=1,size(Ssub,1)
!write(*,'("Ssub",2i3,2x,8f10.5)') myrank,i, (Ssub(i,j),j=1,size(Ssub,2))
!end do
!do i=1,size(S,1)
!write(*,'("S   ",2i3,2x,8f10.5)') myrank,i, (S(i,j),j=1,size(S,2))
!end do
!call stop_program("aaa")

!    call timer( ttmp,time(:,it) ); tlabel(it)="ovlp"; it=it+1

!    call distribute_matrix( sl, S, Ssub )
!    it=it+1 ; call timer( ttmp, time(:,it) ) ; tlabel(it)="dist"

! ---
!    allocate( T1(n,n) ) ; T1=0.0d0
!    allocate( T2(n,n) ) ; T2=0.0d0
!    allocate( T3(n,n) ) ; T3=0.0d0
!    call gather_matrix( sl, Ssub, T1, sl%icontxt_a )
! ---

    call pdpotrf( UPLO, n, Ssub, 1, 1, sl%desca, ierr ) ! Cholesky factorization

!    call gather_matrix( sl, Ssub, T2, sl%icontxt_a )

!    call timer( ttmp,time(:,it) ); tlabel(it)="pdpotrf"; it=it+1

    call pdtrtri( UPLO, 'N', n, Ssub, 1, 1, sl%desca, ierr )

!    call gather_matrix( sl, Ssub, T3, sl%icontxt_a )

!    call timer( ttmp,time(:,it) ); tlabel(it)="pdtrtri"; it=it+1

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
    allocate( S(n,n) ) ; S=zero     !(a1)
    call gather_matrix( sl, Ssub, S, sl%icontxt_a )
    call timer( ttmp,time(:,it) ); tlabel(it)="gather"; it=it+1
#ifdef _TEST_OLD_
    allocate( utmp(m,size(u,2)) ); utmp=u
    call reconstruct_allgatherv( utmp, v )
    call timer( ttmp,time(:,it) ); tlabel(it)="gathe2"; it=it+1
    if ( UPLO == 'U' .or. UPLO == 'u' ) then
       call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, utmp, m )
    else
       call dtrmm( 'R', 'L', 'T', 'N', m, n, one, S, n, utmp, m )
    end if
    call timer( ttmp,time(:,it) ); tlabel(it)="dtrmm"; it=it+1
    u(:,1:n1-n0+1) = utmp(:,n0:n1)
    deallocate( utmp )
#else
    !call dtrmm_bp( UPLO, u, v, S )
    call dtrmm_bp( UPLO, sl, u, v, Ssub, S )
    call timer( ttmp, time(:,it) ); tlabel(it)="adtrmmbp"; it=it+1
#endif
    deallocate( S )
# else

!write(*,*) "before dtrmm_bp",sum(Ssub**2),count(Ssub/=zero),sum(u**2)
    call dtrmm_bp( UPLO, sl, u, Ssub )
!write(*,*) "after dtrmm_bp",sum(Ssub**2),count(Ssub/=zero),sum(u**2)
!call stop_program("bbb")
!    call timer( ttmp,time(:,it) ); tlabel(it)="dtrmmbp"; it=it+1

#endif

! ---

    !deallocate( sl%map_1to2 )
    !deallocate( sl%map_2to1 )
    deallocate( Ssub )

!    call blacs_gridexit( sl%icontxt_a )

!    call timer( ttmp, time(:,it) ); tlabel(it)="finalize"; it=it+1

!    call write_timer( time(:,0:it-1), it-1, tlabel(0:it-1) )

!    call write_border( 1, " d_gram_schmidt_luslbp(end)" )

  end subroutine d_gram_schmidt_luslbp


  subroutine check_lu( check, T )
    implicit none
    real(8),intent(out) :: check(3)
    real(8),intent(in)  :: T(:,:)
    integer :: i,j
    do j=1,size(T,2)
    do i=1,size(T,1)
       if ( i <  j ) check(1)=check(1)+abs( T(i,j) )
       if ( j <  i ) check(2)=check(2)+abs( T(i,j) )
       if ( i == j ) check(3)=check(3)+abs( T(i,j) )
    end do
    end do
  end subroutine check_lu


  subroutine clear_lu( indx, T )
    implicit none
    character(1),intent(in) :: indx
    real(8),intent(inout)  :: T(:,:)
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
  end subroutine clear_lu


end module gram_schmidt_luslbp_module
