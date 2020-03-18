module gram_schmidt_lusl_module

  use io_tools_module, only: IOTools_readIntegerKeyword, IOTools_readStringKeyword
  use sl_tools_module, only: slinfo, distribute_matrix, gather_matrix, restore_param_sl
!  use watch_module, only: timer => watchb, write_timer => write_watchb
  use rgrid_variables, only: dV  
  use parallel_module, only: comm_grid

  implicit none

  private
  public :: gram_schmidt_lusl

  integer :: myrank, nprocs
  integer :: icontxt
  integer :: icontxt_sys
  integer :: desca(9), descz(9)
  integer :: LDR,LDC
  character(1) :: UPLO='U'

  logical :: flag_init = .false.

  interface gram_schmidt_lusl
     module procedure d_gram_schmidt_lusl
  end interface

!  real(8) :: time(2,0:9),ttmp(2)
!  character(10) :: tlabel(0:9)

  type(slinfo) :: sl

contains

  subroutine init_lusl( n )

    implicit none
    integer,intent(in) :: n
    integer :: i,j,k,l,ierr
    integer,allocatable :: usermap(:,:)
    integer,external :: NUMROC
    include 'mpif.h'
    integer,parameter :: unit=90

    if ( flag_init ) return

!    call write_border( 1, " init_lusl(start)" )

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

!    write(*,*) "sl%nprow=",sl%nprow
!    write(*,*) "sl%npcol=",sl%npcol
!    write(*,*) "sl%mbsize=",sl%mbsize
!    write(*,*) "sl%nbsize=",sl%nbsize

! ---

    call blacs_get( 0, 0, icontxt_sys )

    icontxt = icontxt_sys

    allocate( usermap(0:sl%nprow-1,0:sl%npcol-1) ); usermap=0
    i=nprocs/(sl%nprow*sl%npcol)
    k=0
    do j=0,i-1
       l=j*sl%nprow*sl%npcol
       if ( myrank >= l ) k=l
    end do
    k=k-1
    do i=0,sl%nprow-1
    do j=0,sl%npcol-1
       k=k+1   
       usermap(i,j) = k
    end do
    end do
    call blacs_gridmap( icontxt, usermap, sl%nprow, sl%nprow, sl%npcol )
    deallocate( usermap )

    call blacs_gridinfo( icontxt, sl%nprow, sl%npcol, sl%myrow, sl%mycol )

    LDR = NUMROC( n, sl%mbsize, sl%myrow, 0, sl%nprow )
    LDC = NUMROC( n, sl%nbsize, sl%mycol, 0, sl%npcol )

    call descinit( desca,n,n,sl%mbsize,sl%nbsize,0,0,icontxt,LDR,ierr )
    call descinit( descz,n,n,sl%mbsize,sl%nbsize,0,0,icontxt,LDR,ierr )

    flag_init = .true.

!    call write_border( 1, " init_lusl(end)" )

  end subroutine init_lusl


  subroutine d_gram_schmidt_lusl( u )
    implicit none
    real(8),intent(inout)  :: u(:,:)
    real(8),allocatable :: S(:,:),Ssub(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    integer :: m,n,ierr,i,j,it
    integer,allocatable :: ipiv(:)
    include 'mpif.h'

!    call write_border( 1, " d_gram_schmidt_lusl(start)" )
    
!    call timer( ttmp,barrier='on' ); time(:,:)=0.0d0; it=0

    m  = size( u, 1 )
    n  = size( u, 2 )

    call init_lusl( n )

    allocate( S(n,n)        ); S=zero
    allocate( Ssub(LDR,LDC) ); Ssub=zero

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="init"; it=it+1

    call dsyrk( UPLO, 'C', n, m, dV, u, m, zero, S, n )

    call MPI_Allreduce &
         ( MPI_IN_PLACE, S, size(S), MPI_REAL8, MPI_SUM, comm_grid, ierr )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dsyrnk"; it=it+1

    call distribute_matrix( sl, S, Ssub )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dist"; it=it+1

    call PDPOTRF( UPLO, n, Ssub, 1, 1, desca, ierr ) ! Cholesky factorization

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="pdpotrf"; it=it+1

    call PDTRTRI( UPLO, 'N', n, Ssub, 1, 1, desca, ierr )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="pdtrtri"; it=it+1

    call gather_matrix( sl, Ssub, S, icontxt )

    deallocate( Ssub )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="gather"; it=it+1

!    allocate( utmp(m,n) ); utmp=u
!    call DGEMM( 'N', 'N', m, n, n, one, utmp, m, S, n, zero, u, m )
!    u = utmp
!    deallocate( utmp )
!    call timer( ttmp,time(:,it) ); tlabel(it)="dgemm"; it=it+1

    if ( UPLO == 'U' .or. UPLO == 'u' ) then
       call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )
    else
       call dtrmm( 'R', 'L', 'T', 'N', m, n, one, S, n, u, m )
    end if

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dtrmm"; it=it+1

    deallocate( S )

!    call timer( ttmp,time(:,it) ); tlabel(it)="finalize"; it=it+1

!    if ( myrank == 0 ) then
!    call write_timer( time(:,0:it-1), it-1, tlabel(0:it-1) )
!    end if

!    call write_border( 1, " d_gram_schmidt_lusl(end)" )

  end subroutine d_gram_schmidt_lusl

#ifdef TSET
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
#endif

end module gram_schmidt_lusl_module
