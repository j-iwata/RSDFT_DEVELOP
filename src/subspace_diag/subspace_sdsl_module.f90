module subspace_sdsl_module

  use io_tools_module
  use sl_tools_module
  use watch_module, only: timer => watchb, write_timer => write_watchb
  use rgrid_variables, only: dV
  use parallel_module, only: comm_grid, comm_band, MB_d, id_band, ir_band, id_grid
  use hamiltonian_module, only: hamiltonian
  use calc_overlap_sd_module, only: calc_overlap_sd
  use rsdft_mpi_module, only: rsdft_allreduce
  use memory_module, only: check_memory

  implicit none

  private
  public :: subspace_sdsl

  integer :: nband
  integer :: myrank, nprocs
  integer :: icontxt_sys
  integer :: LDR,LDC
  integer :: NP,NQ,TRILWMIN
  character(1) :: UPLO='L'
  integer :: myrnk_g, myrnk_b, np_grid, np_band
#ifdef _DRSDFT_
  real(8),allocatable :: Hsub(:,:)
  real(8),allocatable :: Vsub(:,:)
#else
  complex(8),allocatable :: Hsub(:,:)
  complex(8),allocatable :: Vsub(:,:)
#endif
  logical :: has_init_done  = .false.
  logical :: has_init1_done = .false.
  logical :: has_init2_done = .false.

  real(8) :: time(2,0:9),ttmp(2)
  character(10) :: tlabel(0:9)

  type(slinfo) :: sl

  integer,external :: NUMROC

  integer :: iparam_sdsl(2)=0
  integer :: nband_backup=0

contains

  subroutine init1_sdsl( n )

    implicit none
    integer,intent(in) :: n
    integer :: i,j,k,p,c,r,ierr
    integer :: i_gridmap_group, num_gridmap_groups
    integer,allocatable :: usermap(:,:,:)
    include 'mpif.h'
    integer,parameter :: unit=90

    if ( nband_backup == n ) then
      return
    else
      if ( nband_backup /= 0 ) call blacs_gridexit(sl%icontxt_a)
      nband_backup=n
    end if

    call write_border( 0, " init1_sdsl(start)" )

    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size( comm_grid, np_grid, ierr )
    call MPI_Comm_rank( comm_grid, myrnk_g, ierr )
    call MPI_Comm_size( comm_band, np_band, ierr )
    call MPI_Comm_rank( comm_band, myrnk_b, ierr )

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
!    call IOTools_readIntegerKeyword( "SDPARAM", iparam_sdsl )
!    if ( iparam_sdsl(2) /= 0 ) nblk_ovlp=iparam_sdsl(2)

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
    sl%map_1to2(1,myrank)=myrnk_g+1
    sl%map_1to2(2,myrank)=myrnk_b+1
    call rsdft_allreduce( sl%map_1to2 )

    num_gridmap_groups = nprocs/(np_grid*np_band)

    allocate( usermap(0:sl%nprow-1,0:sl%npcol-1,num_gridmap_groups) ); usermap=0

    do p = 0, nprocs-1
      i = sl%map_1to2(1,p)-1
      j = sl%map_1to2(2,p)-1
      k = i + j*np_grid
      r = k/sl%npcol
      c = k - r*sl%npcol
      if ( k <= sl%nprow*sl%npcol-1 ) then
        i_gridmap_group = p/(np_grid*np_band) + 1
        usermap(r,c,i_gridmap_group)=p
        !if(myrank==0)write(*,'(i4,2i6,2x,i6,2x,2i6,2x,i6)') p, i, j, k,r,c, c+r*sl%npcol
      end if
    end do

    do i = 1, num_gridmap_groups
      i_gridmap_group = myrank/(np_grid*np_band) + 1
      if ( i == i_gridmap_group ) then
        call blacs_gridmap( sl%icontxt_a, usermap(0,0,i), sl%nprow, sl%nprow, sl%npcol )
      end if
    end do

    deallocate( usermap )
    deallocate( sl%map_1to2 )

    call blacs_gridinfo( sl%icontxt_a, sl%nprow, sl%npcol, sl%myrow, sl%mycol )

    LDR = NUMROC( n, sl%mbsize, sl%myrow, 0, sl%nprow )
    LDC = NUMROC( n, sl%nbsize, sl%mycol, 0, sl%npcol )

    call descinit( sl%desca,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )
    call descinit( sl%descz,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )

    NP = NUMROC( n, sl%mbsize, sl%myrow, 0, sl%nprow )
    NQ = NUMROC( n, sl%mbsize, sl%mycol, 0, sl%npcol )
    TRILWMIN = 3*n + max( sl%mbsize*(NP+1), 3*sl%mbsize )

    has_init1_done = .true.

    call write_border( 0, " init1_sdsl(end)" )

  end subroutine init1_sdsl


  subroutine init2_sdsl( n )
    implicit none
    integer,intent(in) :: n
    integer,allocatable :: usermap(:,:,:)
    integer :: i,j,k,p,mgs,mbs,ierr,npr,npc,myr,myc
    integer :: np_grid, np_band
    integer :: num_gridmap_groups,i_gridmap_group
    include 'mpif.h'

    if ( has_init2_done ) return

    !call write_border( 1, " init2_sdsl(start)" )

    num_gridmap_groups = nprocs/(np_grid*np_band)

    allocate( usermap(0:np_grid-1,0:np_band-1,num_gridmap_groups) ); usermap=0

    do p = 0, nprocs-1
      i = sl%map_1to2(1,p)-1
      j = sl%map_1to2(2,p)-1
      i_gridmap_group = p/(np_grid*np_band) + 1
      usermap(i,j,i_gridmap_group)=p
    end do

    sl%icontxt_b = icontxt_sys
    do i = 1, num_gridmap_groups
      i_gridmap_group = myrank/(np_grid*np_band)
      call blacs_gridmap( sl%icontxt_b, usermap(0,0,i_gridmap_group), np_grid, np_grid, np_band )
    end do

    deallocate( usermap )

    call blacs_gridinfo( sl%icontxt_b, npr,npc,myr,myc )

    mgs = n/npr
    mbs = n/npc

    if ( mgs*npr < n ) mgs=mgs+1
    if ( mbs*npc < n ) mbs=mbs+1

    !call descinit( sl%descb,n,n,mgs,mbs,0,0,sl%icontxt_b,mgs,ierr )

    has_init2_done = .true.

    !call write_border( 1, " init2_sdsl(end)" )

  end subroutine init2_sdsl


  subroutine get_nband( nb )
    implicit none
    integer,intent(out) :: nb
    nb = sum(ir_band)
  end subroutine get_nband


  subroutine subspace_sdsl( k, s, u, e )
    implicit none
    integer,intent(in) :: k,s
#ifdef _DRSDFT_
    real(8),intent(inout)  :: u(:,:)
    real(8),allocatable :: H(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),intent(inout)  :: u(:,:)
    complex(8),allocatable :: H(:,:),utmp(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    real(8),intent(inout)  :: e(:)
    integer :: m,n,ierr,i,j,it,n0,n1
    integer,allocatable :: ipiv(:)
    integer :: ib0,ib1,ib2,ig1,ig2,jb1,jb2

    call write_border( 1, " subspace_sdsl(start)" )

!    call timer( ttmp ); time(:,:)=0.0d0; it=0

    m = size( u, 1 )
    n = size( u, 2 )

    call get_nband( nband )
    call init1_sdsl( nband )
    !call init2_sdsl( nband )

! --- Hamiltonian

    !call check_memory( 'dz', m, n )
    allocate( utmp(m,n) ); utmp=zero

    ig1 = id_grid(myrnk_g) + 1
    ig2 = ig1 + m - 1
    ib0 = id_band(myrnk_b)

    do ib1 = 1, n, MB_d
       ib2 = min(ib1+MB_d-1,n)
       jb1 = ib1+ib0
       jb2 = ib2+ib0
       call hamiltonian( u(:,ib1:ib2), utmp(:,ib1:ib2), ib1,k,s )
    end do

!    call timer( ttmp,time(:,it) ); tlabel(it)="ini3"; it=it+1

    !call check_memory( 'dz', nband, nband )
    allocate( H(nband,nband) ); H=zero

    select case( iparam_sdsl(1) )
    case( 0 )
       call calc_overlap_sd( u, utmp, dV, H )
!    case( 1 )
!       call calc_overlap_sd_blk( u, utmp, dV, H )
!    case( 2 )
!       call calc_overlap_sd_dsyr2k( u, utmp, dV, H )
!    case default
!       write(*,*) "iparam_sdsl(1)=",iparam_sdsl(1)," is not available"
    end select

!    call timer( ttmp,time(:,it) ); tlabel(it)="ovlp"; it=it+1

    deallocate( utmp )

    !call check_memory( 'dz', LDR, LDC )
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

    call write_border( 1, " subspace_sdsl(end)" )

  end subroutine subspace_sdsl


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
#ifdef _DRSDFT_
    real(8),intent(inout) :: Hsub(:,:)
    real(8),intent(inout) :: Vsub(:,:)
#else
    complex(8),intent(inout) :: Hsub(:,:)
    complex(8),intent(inout) :: Vsub(:,:)
#endif
    real(8),intent(inout) :: e(:)
    integer,save :: LIWORK=0
    integer,save :: LRWORK=0
    integer,save :: LZWORK=0
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
    complex(8),allocatable :: zwork(:)
    integer :: itmp(1),ierr
    real(8) :: rtmp(1)
    complex(8) :: ztmp(1)
    logical :: disp_sw
    if ( LRWORK == 0 ) then
      LRWORK = max( 1+6*nband+2*NP*NQ, TRILWMIN ) + 2*nband
      LIWORK = 7*nband + 8*sl%npcol + 2
      call check_disp_switch( disp_sw, 0 )
#ifdef _DRSDFT_
      call PDSYEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                   ,sl%descz,rtmp,-1,itmp,LIWORK,ierr)
      if ( disp_sw ) then
        write(*,'("LRWORK,rtmp(1),LIWORK,itmp(1)",4i10)') LRWORK,nint(rtmp(1)),LIWORK,itmp(1)
      end if
#else
      call PZHEEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                   ,sl%descz,ztmp,-1,rtmp,-1,itmp,-1,ierr)
      LZWORK = nint( real(ztmp(1)) )
      if ( disp_sw ) then
        write(*,'("LRWORK,rtmp(1),LIWORK,itmp(1)LZWORK,ztmp(1)",6i10)') &
          LRWORK,nint(rtmp(1)),LIWORK,itmp(1),LZWORK,nint(real(ztmp(1)))
      end if
#endif
      if ( nint(rtmp(1)) > LRWORK ) then
        LRWORK=rtmp(1)
      else
        ! array size query seems not to work well in the mkl routine. 
        ! The following margin is for avoiding the bug
        LRWORK = max( nint(rtmp(1))*10, LRWORK )
        if ( disp_sw ) write(*,'("LRWORK is replaced: LRWORK=",i10)') LRWORK 
      end if
      LIWORK = max( itmp(1), LIWORK )
    end if
    allocate( rwork(LRWORK) ); rwork=0.0d0
    allocate( iwork(LIWORK) ); iwork=0
#ifdef _DRSDFT_
    call PDSYEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                ,sl%descz,rwork,LRWORK,iwork,LIWORK,ierr)
#else
    allocate( zwork(LZWORK) ); zwork=(0.0d0,0.0d0)
    call PZHEEVD('V',UPLO,nband,Hsub,1,1,sl%desca,e,Vsub,1,1 &
                 ,sl%descz,zwork,LZWORK,rwork,LRWORK,iwork,LIWORK,ierr)
    deallocate( zwork )
#endif
    deallocate( iwork )
    deallocate( rwork )
  end subroutine solv_sub


  subroutine rotv_sub( u, R )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: R(:,:)
    real(8),allocatable :: utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
#else
    complex(8),intent(inout) :: u(:,:)
    complex(8),intent(in) :: R(:,:)
    complex(8),allocatable :: utmp(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
#endif
    integer :: m,n,nb
    m =size(u,1)
    n =size(u,2)
    nb=size(R,1)
    allocate( utmp(m,n) ); utmp=u
#ifdef _DRSDFT_
    call DGEMM('N','N',m,n,nb,one,utmp,m,R,nb,zero,u,m)
#else
    call ZGEMM('N','N',m,n,nb,one,utmp,m,R,nb,zero,u,m)
#endif
    deallocate( utmp )
  end subroutine rotv_sub


end module subspace_sdsl_module
