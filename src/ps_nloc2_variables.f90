module ps_nloc2_variables

  use parallel_module, only: MPI_REAL8,RSDFT_MPI_COMPLEX16,nprocs_g
  use memory_module, only: check_memory
  use io_tools_module, only: IOTools_readReal8Keyword

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
  real(8),allocatable :: xVk(:,:,:),yVk(:,:,:),zVk(:,:,:)
  real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  real(8),parameter :: zero=0.d0
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),allocatable :: uVk(:,:,:),sbufnl(:,:),rbufnl(:,:)
  complex(8),allocatable :: sbufnl1(:),rbufnl1(:)
  complex(8),allocatable :: sbufnl3(:,:,:),rbufnl3(:,:,:)
  complex(8),allocatable :: xVk(:,:,:),yVk(:,:,:),zVk(:,:,:)
  complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
  integer,parameter :: TYPE_MAIN=RSDFT_MPI_COMPLEX16
#endif

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
#ifdef _DRSDFT_
  real(8),allocatable,private :: uVunk_backup(:,:,:,:)
#else
  complex(8),allocatable,private :: uVunk_backup(:,:,:,:)
#endif

  logical,public :: FLAG_KEEP_JJ_MAP=.false.

contains

  subroutine read_fmax_conv_ps_nloc2
    implicit none
    real(8) :: fmax_conv
    fmax_conv=0.0d0
    call IOTools_readReal8Keyword( "FMAXCONV", fmax_conv )
    if ( fmax_conv > 0.0d0 ) FLAG_KEEP_JJ_MAP=.true.
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
      if ( disp_sw ) write(*,*) "(rbufnl1,sbufnl1)"
      call check_memory( 8.0d0, n, 2 )
    case( 2 )
      n2=maxval(nrlma_xyz)
      n3=6
      allocate( rbufnl3(n,n2,n3) ) ; rbufnl3=zero
      allocate( sbufnl3(n,n2,n3) ) ; sbufnl3=zero
      if ( disp_sw ) write(*,*) "(rbufnl3,sbufnl3)"
      call check_memory( 8.0d0, n*n2*n3, 2 )
    end select
    allocate( uVunk(nzlma,MB_d)  ) ; uVunk=zero
    allocate( uVunk0(nzlma,MB_d) ) ; uVunk0=zero
    if ( disp_sw ) write(*,*) "(uVunk,uVunk0)"
    call check_memory( 8.0d0, nzlma, MB_d, 2 )
  end subroutine allocate_ps_nloc2

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
    if ( allocated(uVunk_backup) ) then
       if ( size(uVunk_backup,1) == nzlma ) return
       deallocate(uVunk_backup)
    end if
    allocate( uVunk_backup(nzlma,mb0:mb1,mk0:mk1,ms0:ms1) )
    uVunk_backup =zero
  end subroutine prep_backup_uVunk_ps_nloc2


  subroutine backup_uVunk_ps_nloc2( wtmp5 )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(in) :: wtmp5(0:,:,:,:,:)
#else
    complex(8),intent(in) :: wtmp5(0:,:,:,:,:)
#endif
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
!$OMP workshare
    uVunk_backup(:,:,:,:)=wtmp5(0,:,:,:,:)
!$OMP end workshare
  end subroutine backup_uVunk_ps_nloc2


  subroutine restore_uVunk_ps_nloc2( uVunk_out, mb0,mb1,k,s )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(out) :: uVunk_out(:,:)
#else
    complex(8),intent(out) :: uVunk_out(:,:)
#endif
    integer,intent(in) :: mb0,mb1,k,s
    if ( .not.flag_backup_uVunk_ps_nloc2 ) return
!$OMP workshare
    uVunk_out(:,:) = uVunk_backup(:,mb0:mb1,k,s)
!$OMP end workshare
  end subroutine restore_uVunk_ps_nloc2


end module ps_nloc2_variables
