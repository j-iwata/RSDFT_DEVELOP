module subspace_sdsl_module

  use io_tools_module
  use sl_tools_module
  use watch_module, only: timer => watchb, write_timer => write_watchb
  use rgrid_variables, only: dV
  use parallel_module, only: comm_grid, comm_band, MB_d, id_band, ir_band, id_grid
  use hamiltonian_module, only: hamiltonian
  use calc_overlap_sd_module, only: calc_overlap_sd, calc_overlap_sd_blk, nblk_ovlp, calc_overlap_sd_dsyr2k

  implicit none

  private
  public :: subspace_sdsl

  integer :: nband
  integer :: myrank, nprocs
  integer :: icontxt_sys
  integer :: LDR,LDC
  character(1) :: UPLO='L'
  integer :: myrnk_g, myrnk_b
  real(8),allocatable :: Hsub(:,:)
  real(8),allocatable :: Vsub(:,:)

  logical :: has_init_done  = .false.
  logical :: has_init1_done = .false.
  logical :: has_init2_done = .false.

  real(8) :: time(2,0:9),ttmp(2)
  character(10) :: tlabel(0:9)

  type(slinfo) :: sl

  integer,external :: NUMROC

  interface subspace_sdsl
     module procedure d_subspace_sdsl
  end interface

  integer :: iparam_sdsl(2)=0

contains

  subroutine init1_sdsl( n )

    implicit none
    integer,intent(in) :: n
    integer :: i,j,k,l,ierr
    integer,allocatable :: usermap(:,:)
    include 'mpif.h'
    integer,parameter :: unit=90

    if ( has_init1_done ) return

!    call write_border( 1, " init1_sdsl(start)" )

    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    sl%nprow=0
    sl%npcol=0
    sl%mbsize=0
    sl%nbsize=0

!    open( unit, file='sdsl.in', status='old' )
!    call IOTools_readIntegerKeyword( "NPROW" , sl%nprow, unit )
!    call IOTools_readIntegerKeyword( "NPCOL" , sl%npcol, unit )
!    call IOTools_readIntegerKeyword( "MBSIZE", sl%mbsize, unit )
!    call IOTools_readIntegerKeyword( "NBSIZE", sl%nbsize, unit )
!    call IOTools_readStringKeyword( "UPLO", UPLO, unit )
!    close( unit )

    call IOTools_readIntegerKeyword( "SDPARAM", iparam_sdsl )
    if ( iparam_sdsl(2) /= 0 ) nblk_ovlp=iparam_sdsl(2)

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

!    call write_border( 1, " init1_sdsl(end)" )

  end subroutine init1_sdsl


  subroutine init2_sdsl( n )
    implicit none
    integer,intent(in) :: n
    integer,allocatable :: usermap(:,:)
    integer :: i,j,k,mgs,mbs,ierr,npr,npc,myr,myc
    integer :: np_grid, np_band
    include 'mpif.h'

    if ( has_init2_done ) return

!    call write_border( 1, " init2_sdsl(start)" )

    call MPI_Comm_size( comm_grid, np_grid, ierr )
    call MPI_Comm_rank( comm_grid, myrnk_g, ierr )
    call MPI_Comm_size( comm_band, np_band, ierr )
    call MPI_Comm_rank( comm_band, myrnk_b, ierr )

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

    if ( mgs*npr < n ) mgs=mgs+1
    if ( mbs*npc < n ) mbs=mbs+1

    call descinit( sl%descb,n,n,mgs,mbs,0,0,sl%icontxt_b,mgs,ierr )

    deallocate( usermap )

    has_init2_done = .true.

!    call write_border( 1, " init2_sdsl(end)" )

  end subroutine init2_sdsl


  subroutine get_nband( nb )
    implicit none
    integer,intent(out) :: nb
    nb = sum(ir_band)
  end subroutine get_nband


  subroutine d_subspace_sdsl( k, s, u, e )
    implicit none
    integer,intent(in) :: k,s
    real(8),intent(inout)  :: u(:,:)
    real(8),intent(inout)  :: e(:)
    real(8),allocatable :: H(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    integer :: m,n,ierr,i,j,it,n0,n1
    integer,allocatable :: ipiv(:)
    integer :: ib0,ib1,ib2,ig1,ig2,jb1,jb2

!    call write_border( 1, " d_subspace_sdsl(start)" )

!    call timer( ttmp ); time(:,:)=0.0d0; it=0

    m = size( u, 1 )
    n = size( u, 2 )

    if ( .not.has_init_done ) then
       call get_nband( nband )
       call init1_sdsl( nband )
       call init2_sdsl( nband )
       !call init_overlap_bp( m, n, dV, comm_grid, comm_band )
       !call init_dtrmm_bp( m, n, comm_grid, comm_band )
       has_init_done = .true.
    end if

! --- Hamiltonian

    allocate( utmp(m,n) ); utmp=zero

    ig1 = id_grid(myrnk_g) + 1
    ig2 = ig1 + m - 1
    ib0 = id_band(myrnk_b)

    do ib1 = 1, n, MB_d
       ib2 = min(ib1+MB_d-1,n)
       jb1 = ib1+ib0
       jb2 = ib2+ib0
       call hamiltonian( k,s,u(:,ib1:ib2),utmp(:,ib1:ib2),jb1,jb2,ig1,ig2 )
    end do

!    call timer( ttmp,time(:,it) ); tlabel(it)="ini3"; it=it+1

    allocate( H(nband,nband) ); H=zero

    select case( iparam_sdsl(1) )
    case( 0 )
       call calc_overlap_sd( u, utmp, dV, H )
    case( 1 )
       call calc_overlap_sd_blk( u, utmp, dV, H )
    case( 2 )
       call calc_overlap_sd_dsyr2k( u, utmp, dV, H )
!    case default
!       write(*,*) "iparam_sdsl(1)=",iparam_sdsl(1)," is not available"
    end select

!    call timer( ttmp,time(:,it) ); tlabel(it)="ovlp"; it=it+1

    deallocate( utmp )

    allocate( Hsub(LDR,LDC) ); Hsub=zero

    call distribute_matrix( sl, H, Hsub )

!    call timer( ttmp,time(:,it) ); tlabel(it)="dist"; it=it+1

    deallocate( H )

    allocate( Vsub(LDR,LDC) ); Vsub=zero

    call solv_sub( Hsub, Vsub, e )

!    call timer( ttmp,time(:,it) ); tlabel(it)="dist"; it=it+1

    deallocate( Hsub )

    allocate( H(nband,nband) ); H=zero

    call gather_matrix( sl, Vsub, H, sl%icontxt_a )

!    call timer( ttmp,time(:,it) ); tlabel(it)="pdpotrf"; it=it+1

    deallocate( Vsub ) 

    call rotv_sub( u, H )

!    call timer( ttmp,time(:,it) ); tlabel(it)="pdtrtri"; it=it+1

! ---

    !deallocate( sl%map_1to2 )
    !deallocate( sl%map_2to1 )

    deallocate( H )

!    call blacs_gridexit( sl%icontxt_a )

!    call timer( ttmp, time(:,it) ); tlabel(it)="finalize"; it=it+1

!    call write_timer( time(:,0:it-1), it-1, tlabel(0:it-1) )

!    call write_border( 1, " d_subspace_sdsl(end)" )

  end subroutine d_subspace_sdsl


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


  subroutine solv_sub( Hsub, Vsub, e )
    implicit none
    real(8),intent(inout) :: Hsub(:,:)
    real(8),intent(inout) :: Vsub(:,:)
    real(8),intent(inout) :: e(:)
    integer,save :: LRWORK=0
    integer,save :: LIWORK=0
    integer :: itmp(1),ierr
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
    real(8) :: rtmp(1)
    if ( LRWORK == 0 ) then
       LIWORK = 7*nband + 8*sl%npcol + 2
       call PDSYEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                   ,sl%descz,rtmp,-1,itmp,LIWORK,ierr)
       LRWORK = nint( rtmp(1) )*2
       LIWORK = max( itmp(1), LIWORK )
    end if
    allocate( rwork(LRWORK) ); rwork=0.0d0
    allocate( iwork(LIWORK) ); iwork=0
    call PDSYEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                ,sl%descz,rwork,LRWORK,iwork,LIWORK,ierr)
    deallocate( iwork )
    deallocate( rwork )
  end subroutine solv_sub


  subroutine rotv_sub( u, R )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: R(:,:)
    real(8),allocatable :: utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    integer :: m,n,nb
    m =size(u,1)
    n =size(u,2)
    nb=size(R,1)
    allocate( utmp(m,n) ); utmp=u
    call DGEMM('N','N',m,n,nb,one,utmp,m,R,nb,zero,u,m)
    deallocate( utmp )
  end subroutine rotv_sub


end module subspace_sdsl_module
