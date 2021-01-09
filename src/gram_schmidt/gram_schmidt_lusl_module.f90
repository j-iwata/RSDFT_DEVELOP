module gram_schmidt_lusl_module

  use io_tools_module, only: IOTools_readIntegerKeyword, IOTools_readStringKeyword
  use sl_tools_module, only: slinfo, distribute_matrix, gather_matrix, restore_param_sl
!  use watch_module, only: timer => watchb, write_timer => write_watchb
  use rgrid_variables, only: dV  
  use parallel_module, only: comm_grid
  use dsyrk_module, only: calc_dsyrk3, ialgo_dsyrk, nblk_dsyrk, calc_zherk3
  use trmm_module, only: calc_ztrmm3
  use rsdft_mpi_module, only: rsdft_allreduce

  implicit none

  private
  public :: gram_schmidt_lusl

  integer :: myrank, nprocs
  integer :: icontxt
  integer :: icontxt_sys
  integer :: desca(9), descz(9)
  integer :: LDR,LDC
  character(1) :: UPLO='U'

  integer,public :: iparam_gs_lusl(9)=1
  integer :: ialgo_dtrmm = 1
  integer :: nblk_dtrmm = 1

  logical :: flag_init = .false.

  interface gram_schmidt_lusl
     module procedure d_gram_schmidt_lusl, z_gram_schmidt_lusl
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

    ialgo_dsyrk = iparam_gs_lusl(1)
    nblk_dsyrk  = iparam_gs_lusl(3)
    ialgo_dtrmm = iparam_gs_lusl(2)
    nblk_dtrmm  = iparam_gs_lusl(4)

    flag_init = .true.

!    call write_border( 1, " init_lusl(end)" )

  end subroutine init_lusl


  subroutine calc_dsyrk( u, S )
    implicit none
    real(8),intent(in) :: u(:,:)
    real(8),intent(inout) :: S(:,:)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: utmp(:,:)
    integer :: m,n
    m=size(u,1)
    n=size(u,2)
    call dsyrk( 'U', 'C', n, m, dV, u, m, zero, S, n )
!    allocate( utmp(n,m) ); utmp=zero !transpose(u)
!    call tr
!    call DGEMM( 'N', 'N', n, n, m, dV, utmp, n, u, m, zero, S, n )
!    call DGEMM( 'T', 'N', n, n, m, dV, u, m, u, m, zero, S, n )
!    deallocate( utmp )
  contains
    subroutine tr
       utmp=transpose(u)
    end subroutine tr
  end subroutine calc_dsyrk

  subroutine calc_dsyrk2( u, S )
    implicit none
    real(8),intent(in) :: u(:,:)
    real(8),intent(inout) :: S(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: utmp(:,:)
    integer :: m,n,nblk,iblk,jblk
    integer :: i0,i1,j0,j1,ni,nj

    m=size(u,1)
    n=size(u,2)
    nblk = max( n/nblk_dsyrk, 1 )

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk, n )
       nj = j1 - j0 + 1

       do iblk = 1, jblk, nblk

          i0 = iblk
          i1 = min( iblk+nblk, n )
          ni = i1 - i0 + 1

          if ( iblk == jblk ) then
             if ( ialgo_dsyrk == 1 ) then
                call DGEMM( 'T','N',ni,nj,m,dV,u(1,i0),m,u(1,j0),m,zero,S(i0,j0),n )
             else
                call DSYRK( 'U', 'C', ni, m, dV, u(1,i0), m, zero, S(i0,i0), n )
             end if
          else
             call DGEMM( 'T','N',ni,nj,m,dV,u(1,i0),m,u(1,j0),m,zero,S(i0,j0),n )
          end if

       end do !iblk

    end do !jblk

!    allocate( utmp(n,m) ) !; utmp=transpose(u)
!    call tr
!    call DGEMM( 'N', 'N', n, n, m, dV, utmp, n, u, m, zero, S, n )
!    call DGEMM( 'T', 'N', n, n, m, dV, u, m, u, m, zero, S, n )
!    deallocate( utmp )
!  contains
!    subroutine tr
!       utmp=transpose(u)
!    end subroutine tr
  end subroutine calc_dsyrk2

  subroutine d_gram_schmidt_lusl( u )
    implicit none
    real(8),intent(inout)  :: u(:,:)
    real(8),allocatable :: S(:,:),Ssub(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    integer :: m,n,ierr,i,j,it
    integer,allocatable :: ipiv(:)
    include 'mpif.h'

    call write_border( 1, " d_gram_schmidt_lusl(start)" )
    
!    call timer( ttmp,barrier='on' ); time(:,:)=0.0d0; it=0

    m  = size( u, 1 )
    n  = size( u, 2 )

    call init_lusl( n )

    allocate( S(n,n)        ); S=zero
    allocate( Ssub(LDR,LDC) ); Ssub=zero

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="init"; it=it+1

!    call DSYRK( 'U', 'C', n, m, dV, u, m, zero, S, n )
!    call DGEMM('T','N',n,n,m,dV,u,m,u,m,zero,S,n)
!    call calc_dsyrk( u, S )
!    call calc_dsyrk2( u, S )
    call calc_dsyrk3( u, S )

    do j=1,n
    do i=j+1,n
       S(i,j)=zero
    end do
    end do

!    write(*,*) sum((Ssub-S)**2),count(S/=zero),count(Ssub/=zero)

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
!
!    if ( UPLO == 'U' .or. UPLO == 'u' ) then
!       call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )
!    else
!       call dtrmm( 'R', 'L', 'T', 'N', m, n, one, S, n, u, m )
!    end if

!    call calc_dtrmm( u, S )
    call calc_dtrmm3( u, S )
!    call calc_dtrmm4( u, S )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dtrmm"; it=it+1

    deallocate( S )

!    call timer( ttmp,time(:,it) ); tlabel(it)="finalize"; it=it+1

!    if ( myrank == 0 ) then
!    call write_timer( time(:,0:it-1), it-1, tlabel(0:it-1) )
!    end if

    call write_border( 1, " d_gram_schmidt_lusl(end)" )

  end subroutine d_gram_schmidt_lusl


  subroutine calc_dtrmm( u, S )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: S(:,:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n
    real(8),allocatable :: utmp(:,:)
    m=size(u,1)
    n=size(u,2)
    call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )
!    call alloc_utmp
!    utmp=u
!    call DGEMM( 'N', 'N', m, n, n, one, utmp, m, S, n, zero, u, m )
!    call dealloc_utmp
  contains
    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      deallocate( utmp )
    end subroutine dealloc_utmp
  end subroutine calc_dtrmm


  subroutine calc_dtrmm3( u, S )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: S(:,:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n,j0,j1,i0,i1,iblk,jblk
    integer :: nblk,nj,ni
    real(8),allocatable :: utmp(:,:)

    m=size(u,1)
    n=size(u,2)

    nblk = max( n/nblk_dtrmm, 1 )

    call alloc_utmp

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk-1, n )
       nj = j1 - j0 + 1

       do iblk = 1, jblk, nblk

          i0 = iblk
          i1 = min( iblk+nblk-1, n )
          ni = i1 - i0 + 1

          if ( iblk == jblk ) then
             if ( ialgo_dtrmm == 1 ) then
                call DGEMM( 'N','N',m,ni,ni,one,u(1,i0),m,S(i0,i0),n,one,utmp(1,i0),m )
             else
                call DTRMM( 'R', 'U', 'N', 'N', m, ni, one, S(i0,i0), n, u(1,i0), m )
             end if
          else
             call DGEMM( 'N','N',m,nj,ni,one,u(1,i0),m,S(i0,j0),n,one,utmp(1,j0),m )
          end if

       end do !iblk

    end do !jblk

    !u=utmp
    call dealloc_utmp

  contains
    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      u=utmp
      deallocate( utmp )
    end subroutine dealloc_utmp
  end subroutine calc_dtrmm3


  subroutine calc_dtrmm4( u, S )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: S(:,:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n,j0,j1,i0,i1,iblk,jblk
    integer :: nblk,nj,ni
    real(8),allocatable :: utmp(:,:)
    real(8),allocatable :: stmp(:,:)

    m=size(u,1)
    n=size(u,2)

    nblk = max( n/2, 1 )

    call alloc_utmp

    do iblk = 1, n, nblk
       i0 = iblk
       i1 = min( iblk+nblk-1, n )
       ni = i1 - i0 + 1
       utmp(:,i0:i1) = u(:,i0:i1)
       call DTRMM( 'R', 'U', 'N', 'N', m, ni, one, S(i0,i0), n, utmp(1,i0), m )
    end do !iblk

 ! ---

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk-1, n )
       nj = j1 - j0 + 1

       do iblk = 1, jblk-1, nblk

          i0 = iblk
          i1 = min( iblk+nblk-1, n )
          ni = i1 - i0 + 1

          call DGEMM( 'N','N',m,nj,ni,one,u(1,i0),m,S(i0,j0),n,one,utmp(1,j0),m )

       end do !iblk

    end do !jblk

    !u=utmp
    call dealloc_utmp

  contains

    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      u=utmp
      deallocate( utmp )
    end subroutine dealloc_utmp

!    subroutine alloc_stmp
!      implicit none
!      allocate( stmp(i0:i1,i0:i1) ); stmp=S(i0:i1,i0:i1)
!    end subroutine alloc_stmp
!    subroutine dealloc_stmp
!      deallocate( stmp )
!    end subroutine dealloc_stmp

  end subroutine calc_dtrmm4



  subroutine z_gram_schmidt_lusl( u )
    implicit none
    complex(8),intent(inout)  :: u(:,:)
    complex(8),allocatable :: S(:,:),Ssub(:,:),utmp(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    integer :: m,n,ierr,i,j,it
    integer,allocatable :: ipiv(:)

    call write_border( 1, " z_gram_schmidt_lusl(start)" )
    
!    call timer( ttmp,barrier='on' ); time(:,:)=0.0d0; it=0

    m = size( u, 1 )
    n = size( u, 2 )

    if ( .not.flag_init ) call init_lusl( n )

    allocate( S(n,n)        ); S=zero
    allocate( Ssub(LDR,LDC) ); Ssub=zero

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="init"; it=it+1

!    call DSYRK( 'U', 'C', n, m, dV, u, m, zero, S, n )
!    call DGEMM('T','N',n,n,m,dV,u,m,u,m,zero,S,n)
!    call calc_dsyrk( u, S )
!    call calc_dsyrk2( u, S )
    call calc_zherk3( u, S )

    do j=1,n
    do i=j+1,n
       S(i,j)=zero
    end do
    end do

!    write(*,*) sum((Ssub-S)**2),count(S/=zero),count(Ssub/=zero)

    !call MPI_Allreduce &
    !     ( MPI_IN_PLACE, S, size(S), MPI_REAL8, MPI_SUM, comm_grid, ierr )
    call rsdft_allreduce( S, comm_grid )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dsyrnk"; it=it+1

    call distribute_matrix( sl, S, Ssub )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dist"; it=it+1

    call PZPOTRF( UPLO, n, Ssub, 1, 1, desca, ierr ) ! Cholesky factorization

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="pzpotrf"; it=it+1

    call PZTRTRI( UPLO, 'N', n, Ssub, 1, 1, desca, ierr )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="pdtrtri"; it=it+1

    call gather_matrix( sl, Ssub, S, icontxt )

    deallocate( Ssub )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="gather"; it=it+1

!    allocate( utmp(m,n) ); utmp=u
!    call DGEMM( 'N', 'N', m, n, n, one, utmp, m, S, n, zero, u, m )
!    u = utmp
!    deallocate( utmp )
!    call timer( ttmp,time(:,it) ); tlabel(it)="dgemm"; it=it+1
!    if ( UPLO == 'U' .or. UPLO == 'u' ) then
!       call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )
!    else
!       call dtrmm( 'R', 'L', 'T', 'N', m, n, one, S, n, u, m )
!    end if

!    call calc_ztrmm( u, S )
    call calc_ztrmm3( u, S )
!    call calc_ztrmm4( u, S )

!    call timer( ttmp,time(:,it),'on' ); tlabel(it)="dtrmm"; it=it+1

    deallocate( S )

!    call timer( ttmp,time(:,it) ); tlabel(it)="finalize"; it=it+1

!    if ( myrank == 0 ) then
!    call write_timer( time(:,0:it-1), it-1, tlabel(0:it-1) )
!    end if

    call write_border( 1, " z_gram_schmidt_lusl(end)" )

  end subroutine z_gram_schmidt_lusl


end module gram_schmidt_lusl_module
