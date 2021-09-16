module ps_nloc2_variables

  use parallel_module, only: MPI_REAL8,RSDFT_MPI_COMPLEX16,nprocs_g
  use memory_module, only: check_memory
  use io_tools_module, only: IOTools_readReal8Keyword,IOTools_readIntegerKeyword

  integer :: Mlma,nzlma
  integer,allocatable :: JJ_MAP(:,:,:),MJJ_MAP(:),MJJ(:)
  integer :: MMJJ,MAXMJJ,nrlma_xyz(6)
  real(8),allocatable :: uV(:,:)
  integer,allocatable :: amap(:),lmap(:),mmap(:),iorbmap(:),iamap(:)
  integer,allocatable :: lma_nsend(:),iuV(:),nl_rank_map(:)
  integer,allocatable :: num_2_rank(:,:)
  integer,allocatable :: sendmap(:,:),recvmap(:,:)
  integer,allocatable :: JJP(:,:)
  integer :: nl_max_send
#ifdef _DRSDFT_
  real(8),allocatable :: uVk(:,:,:),sbufnl(:,:),rbufnl(:,:)
  real(8),allocatable :: sbufnl1(:),rbufnl1(:)
  real(8),allocatable :: sbufnl3(:,:,:),rbufnl3(:,:,:)
  real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  real(8),parameter :: zero=0.d0
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),allocatable :: uVk(:,:,:),sbufnl(:,:),rbufnl(:,:)
  complex(8),allocatable :: sbufnl1(:),rbufnl1(:)
  complex(8),allocatable :: sbufnl3(:,:,:),rbufnl3(:,:,:)
  complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
  integer,parameter :: TYPE_MAIN=RSDFT_MPI_COMPLEX16
#endif

  real(8),allocatable :: d_sbufnl(:,:), d_rbufnl(:,:)
  real(8),allocatable :: d_sbufnl1(:), d_rbufnl1(:)
  real(8),allocatable :: d_sbufnl3(:,:,:), d_rbufnl3(:,:,:)
  real(8),allocatable :: d_uVk(:,:,:), d_uVunk(:,:), d_uVunk0(:,:)

  complex(8),allocatable :: z_sbufnl(:,:), z_rbufnl(:,:)
  complex(8),allocatable :: z_sbufnl1(:), z_rbufnl1(:)
  complex(8),allocatable :: z_sbufnl3(:,:,:), z_rbufnl3(:,:,:)
  complex(8),allocatable :: z_uVk(:,:,:), z_uVunk(:,:), z_uVunk0(:,:)

  complex(8),allocatable :: xVk(:,:,:), yVk(:,:,:), zVk(:,:,:)

! spin-orbit
  real(8),allocatable :: uVso(:,:)
  complex(8),allocatable :: uVk_so00(:,:,:)
  complex(8),allocatable :: uVk_so11(:,:,:)
  complex(8),allocatable :: uVk_so12(:,:,:)
  complex(8),allocatable :: uVk_so21(:,:,:)
  complex(8),allocatable :: uVk_so22(:,:,:)

! muti-reference & spin-orbit
  integer :: N_nzqr
  integer,allocatable :: nzqr_pair(:,:)
  real(8),allocatable :: Dij00(:)
  real(8),allocatable :: Kij00(:)

  public :: backup_uVunk_ps_nloc2
  public :: restore_uVunk_ps_nloc2
  public :: prep_backup_uVunk_ps_nloc2
  logical,public :: flag_backup_uVunk_ps_nloc2=.false.

  real(8),allocatable,private :: d_uVunk_backup(:,:,:,:)
  complex(8),allocatable,private :: z_uVunk_backup(:,:,:,:)

  interface backup_uVunk_ps_nloc2
    module procedure d_backup_uVunk_ps_nloc2, z_backup_uVunk_ps_nloc2
  end interface

  interface restore_uVunk_ps_nloc2
    module procedure d_restore_uVunk_ps_nloc2, z_restore_uVunk_ps_nloc2
  end interface

  logical,public :: FLAG_KEEP_JJ_MAP=.false.
  logical,public :: FLAG_KEEP_uV=.false.

contains

  subroutine read_fmax_conv_ps_nloc2
    implicit none
    real(8) :: fmax_conv
    integer :: swopt, swband
    fmax_conv=0.0d0
    swopt=0
    swband=0
    call IOTools_readReal8Keyword( "FMAXCONV", fmax_conv )
    call IOTools_readIntegerKeyword( "SWOPT" , swopt )
    call IOTools_readIntegerKeyword( "SWBAND", swband )
    if ( fmax_conv > 0.0d0 .or. swopt /= 0 ) FLAG_KEEP_JJ_MAP=.true.
    if ( swband /= 0 ) then
      FLAG_KEEP_JJ_MAP=.true.
      FLAG_KEEP_uV=.true.
    end if
  end subroutine read_fmax_conv_ps_nloc2

  subroutine allocate_ps_nloc2( MB_d, itype_nl_sendrecv )
    implicit none
    integer,intent(in) :: MB_d
    integer,optional,intent(in) :: itype_nl_sendrecv
    integer :: n, itype, n2,n3
    logical :: disp_sw
    call check_disp_switch( disp_sw, 0 )
    n=maxval( lma_nsend )*4*MB_d
    if ( allocated(rbufnl)  ) deallocate(rbufnl)
    if ( allocated(sbufnl)  ) deallocate(sbufnl)
    if ( allocated(rbufnl1) ) deallocate(rbufnl1)
    if ( allocated(sbufnl1) ) deallocate(sbufnl1)
    if ( allocated(rbufnl3) ) deallocate(rbufnl3)
    if ( allocated(sbufnl3) ) deallocate(sbufnl3)
    if ( allocated(uVunk)   ) deallocate(uVunk)
    if ( allocated(uVunk0)  ) deallocate(uVunk0)
    itype=0
    if ( present(itype_nl_sendrecv) ) itype=itype_nl_sendrecv
    select case( itype )
    case default
      allocate( sbufnl(n,0:nprocs_g-1) ) ; sbufnl=zero
      allocate( rbufnl(n,0:nprocs_g-1) ) ; rbufnl=zero
    case( 1 )
      allocate( rbufnl1(n) ) ; rbufnl1=zero
      allocate( sbufnl1(n) ) ; sbufnl1=zero
      !if ( disp_sw ) write(*,*) "(rbufnl1,sbufnl1)"
      !call check_memory( 8.0d0, n, 2 )
    case( 2 )
      n2=maxval(nrlma_xyz)
      n3=6
      allocate( rbufnl3(n,n2,n3) ) ; rbufnl3=zero
      allocate( sbufnl3(n,n2,n3) ) ; sbufnl3=zero
      !if ( disp_sw ) write(*,*) "(rbufnl3,sbufnl3)"
      !call check_memory( 8.0d0, n*n2*n3, 2 )
    end select
    allocate( uVunk(nzlma,MB_d)  ) ; uVunk=zero
    allocate( uVunk0(nzlma,MB_d) ) ; uVunk0=zero
    !if ( disp_sw ) write(*,*) "(uVunk,uVunk0)"
    !call check_memory( 8.0d0, nzlma, MB_d, 2 )
  end subroutine allocate_ps_nloc2

  subroutine d_allocate_ps_nloc2( k0, k1 )
    use parallel_module, only: MB_d_nl
    implicit none
    integer,intent(in) :: k0, k1
    integer :: n0,n1,n2,n3
    logical :: disp_on
    real(8),parameter :: z0=0.0d0
    if ( allocated(JJP)       ) deallocate(JJP)
    if ( allocated(d_uVk)     ) deallocate(d_uVk)
    if ( allocated(d_rbufnl1) ) deallocate(d_rbufnl1)
    if ( allocated(d_sbufnl1) ) deallocate(d_sbufnl1)
    if ( allocated(d_rbufnl3) ) deallocate(d_rbufnl3)
    if ( allocated(d_sbufnl3) ) deallocate(d_sbufnl3)
    if ( allocated(d_uVunk)   ) deallocate(d_uVunk)
    if ( allocated(d_uVunk0)  ) deallocate(d_uVunk0)
    n0=0; if ( nzlma > 0 ) n0=maxval( MJJ(1:nzlma) )
    n1=maxval( lma_nsend )*4*MB_d_nl
    n2=maxval(nrlma_xyz)
    n3=6
    ! call check_disp_switch( disp_on, 0 )
    !if ( disp_on ) write(*,*) "(z_rbufnl1,z_sbufnl1)"
    !call check_memory( 8.0d0, n, 2 )
    allocate( JJP(n0,nzlma)           ); JJP=0
    allocate( d_uVk(n0,nzlma,k0:k1)   ); d_uVk=z0
    allocate( d_rbufnl1(n1)           ); d_rbufnl1=z0
    allocate( d_sbufnl1(n1)           ); d_sbufnl1=z0
    allocate( d_rbufnl3(n1,n2,n3)     ); d_rbufnl3=z0
    allocate( d_sbufnl3(n1,n2,n3)     ); d_sbufnl3=z0
    allocate( d_uVunk(nzlma,MB_d_nl)  ); d_uVunk=z0
    allocate( d_uVunk0(nzlma,MB_d_nl) ); d_uVunk0=z0
  end subroutine d_allocate_ps_nloc2

  subroutine z_allocate_ps_nloc2( k0, k1 )
    use parallel_module, only: MB_d_nl
    implicit none
    integer,intent(in) :: k0, k1
    integer :: n0,n1,n2,n3
    logical :: disp_on
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    if ( allocated(JJP)       ) deallocate(JJP)
    if ( allocated(z_uVk)     ) deallocate(z_uVk)
    if ( allocated(z_rbufnl1) ) deallocate(z_rbufnl1)
    if ( allocated(z_sbufnl1) ) deallocate(z_sbufnl1)
    if ( allocated(z_rbufnl3) ) deallocate(z_rbufnl3)
    if ( allocated(z_sbufnl3) ) deallocate(z_sbufnl3)
    if ( allocated(z_uVunk)   ) deallocate(z_uVunk)
    if ( allocated(z_uVunk0)  ) deallocate(z_uVunk0)
    n0=0; if ( nzlma > 0 ) n0=maxval( MJJ(1:nzlma) )
    n1=maxval( lma_nsend )*4*MB_d_nl
    n2=maxval(nrlma_xyz)
    n3=6
    ! call check_disp_switch( disp_on, 0 )
    !if ( disp_sw ) write(*,*) "(z_rbufnl1,z_sbufnl1)"
    !call check_memory( 16.0d0, n, 2 )
    allocate( JJP(n0,nzlma)           ); JJP=0
    allocate( z_uVk(n0,nzlma,k0:k1)   ); z_uVk=z0
    allocate( z_rbufnl1(n1)           ); z_rbufnl1=z0
    allocate( z_sbufnl1(n1)           ); z_sbufnl1=z0
    allocate( z_rbufnl3(n1,n2,n3)     ); z_rbufnl3=z0
    allocate( z_sbufnl3(n1,n2,n3)     ); z_sbufnl3=z0
    allocate( z_uVunk(nzlma,MB_d_nl)  ); z_uVunk=z0
    allocate( z_uVunk0(nzlma,MB_d_nl) ); z_uVunk0=z0
  end subroutine z_allocate_ps_nloc2

  subroutine checkMapsBeforeForce(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: i
    write(5000+myrank,'(2A7)') 'nzlma','amap'
    write(5100+myrank,'(2A7)') 'nzlma','lmap'
    write(5200+myrank,'(2A7)') 'nzlma','mmap'
    write(5300+myrank,'(2A7)') 'nzlma','iorbmap'
    do i=1,nzlma
      write(5000+myrank,'(2I7)') nzlma,amap(i)
      write(5100+myrank,'(2I7)') nzlma,lmap(i)
      write(5200+myrank,'(2I7)') nzlma,mmap(i)
      write(5300+myrank,'(2I7)') nzlma,iorbmap(i)
    enddo
  END subroutine checkMapsBeforeForce


  subroutine prep_backup_uVunk_ps_nloc2( mb0,mb1,mk0,mk1,ms0,ms1 )
    implicit none
    integer,intent(in) :: mb0,mb1,mk0,mk1,ms0,ms1
    flag_backup_uVunk_ps_nloc2=.true.
    if ( allocated(d_uVunk_backup) ) then
       if ( size(d_uVunk_backup,1) == nzlma ) return
       deallocate(d_uVunk_backup)
    end if
    allocate( d_uVunk_backup(nzlma,mb0:mb1,mk0:mk1,ms0:ms1) )
    d_uVunk_backup =0.0d0
  end subroutine prep_backup_uVunk_ps_nloc2


  subroutine d_backup_uVunk_ps_nloc2( wtmp )
    implicit none
    real(8),intent(in) :: wtmp(:,:,:,:)
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
!$omp workshare
    d_uVunk_backup(:,:,:,:)=wtmp(:,:,:,:)
!$omp end workshare
  end subroutine d_backup_uVunk_ps_nloc2

  subroutine z_backup_uVunk_ps_nloc2( wtmp )
    implicit none
    complex(8),intent(in) :: wtmp(:,:,:,:)
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
!$omp workshare
    z_uVunk_backup(:,:,:,:)=wtmp(:,:,:,:)
!$omp end workshare
  end subroutine z_backup_uVunk_ps_nloc2


  subroutine d_restore_uVunk_ps_nloc2( uVunk_out,n,k,s )
    implicit none
    real(8),intent(inout) :: uVunk_out(:,:)
    integer,intent(in) :: n,k,s
    integer :: n0, n1
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
    n0 = n
    n1 = n0 + size( uVunk_out, 2 ) - 1
!$omp workshare
    uVunk_out(:,:) = d_uVunk_backup(:,n0:n1,k,s)
!$omp end workshare
  end subroutine d_restore_uVunk_ps_nloc2

  subroutine z_restore_uVunk_ps_nloc2( uVunk_out,n,k,s )
    implicit none
    complex(8),intent(inout) :: uVunk_out(:,:)
    integer,intent(in) :: n,k,s
    integer :: n0,n1
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
    n0 = n
    n1 = n0 + size( uVunk_out, 2 ) - 1
!$omp workshare
    uVunk_out(:,:) = z_uVunk_backup(:,n0:n1,k,s)
!$omp end workshare
  end subroutine z_restore_uVunk_ps_nloc2


end module ps_nloc2_variables
