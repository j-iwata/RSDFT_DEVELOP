module sl_tools_module

  use rsdft_mpi_module, only: rsdft_allreduce

  implicit none

  private
  public :: fill_matrix
  public :: gather_matrix
  public :: gatherA_matrix
  public :: distribute_matrix
  public :: slinfo
  public :: sl_block_map
  public :: backup_param_sl
  public :: restore_param_sl

  type slinfo
    integer :: myrow,  mycol
    integer :: nprow,  npcol
    integer :: mbsize, nbsize
    integer,allocatable :: map_1to2(:,:)
    integer,allocatable :: map_2to1(:,:)
    integer :: icontxt_a, icontxt_b
    integer :: desca(9), descb(9), descz(9)
  end type slinfo

  type(slinfo) :: sl_backup

  type blkinfo
    integer :: iband(2), jband(2)
    integer :: istor(2), jstor(2)
    integer :: iproc, jproc
  end type blkinfo

  type(blkinfo),allocatable :: blk_map(:,:)

  INTERFACE fill_matrix
     MODULE PROCEDURE d_fill_matrix, z_fill_matrix
  END INTERFACE

  interface gather_matrix
    module procedure d_gather_matrix, z_gather_matrix
  end interface

  INTERFACE gatherA_matrix
     MODULE PROCEDURE d_gatherA_matrix
  END INTERFACE

  interface distribute_matrix
    module procedure d_distribute_matrix, z_distribute_matrix
  end interface

contains


  SUBROUTINE d_fill_matrix( h )
    implicit none
    real(8),intent(INOUT) :: h(:,:)
    integer :: n,i,j
    n=size(h,1)
    do j=1,n
       do i=j+1,n
          h(i,j) = h(j,i)
       end do
    end do
  END SUBROUTINE d_fill_matrix

  SUBROUTINE z_fill_matrix( h )
    implicit none
    complex(8),intent(INOUT) :: h(:,:)
    integer :: n,i,j
    n=size(h,1)
    do j=1,n
       do i=j+1,n
          h(i,j) = conjg( h(j,i) )
       end do
    end do
  END SUBROUTINE z_fill_matrix


  subroutine d_distribute_matrix( sl, H, Hsub )
    implicit none
    type(slinfo),intent(in) :: sl
    real(8),intent(in)      :: H(:,:)
    real(8),intent(out)     :: Hsub(:,:)
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    m = size( H, 1 )
    n = size( H, 2 )
    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
      j1=min(j0+sl%nbsize-1,n)
      jj=j1-j0+1
      jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
      ip=-1
      is= 0
      do i0=1,m,sl%mbsize
        i1=min(i0+sl%mbsize-1,m)
        ii=i1-i0+1
        ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
        do j=1,jj
        do i=1,ii
          Hsub(is+i,js+j) = H(i0+i-1,j0+j-1)
        end do
        end do
        is = is + ii
      end do
      js = js + jj
    end do
  end subroutine d_distribute_matrix

  subroutine z_distribute_matrix( sl, H, Hsub )
    implicit none
    type(slinfo),intent(in) :: sl
    complex(8),intent(in)   :: H(:,:)
    complex(8),intent(out)  :: Hsub(:,:)
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    m = size( H, 1 )
    n = size( H, 2 )
    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
      j1=min(j0+sl%nbsize-1,n)
      jj=j1-j0+1
      jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
      ip=-1
      is= 0
      do i0=1,m,sl%mbsize
        i1=min(i0+sl%mbsize-1,m)
        ii=i1-i0+1
        ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
        do j=1,jj
        do i=1,ii
          Hsub(is+i,js+j) = H(i0+i-1,j0+j-1)
        end do
        end do
        is = is + ii
      end do
      js = js + jj
    end do
  end subroutine z_distribute_matrix


  subroutine d_gather_matrix( sl, Hsub, H, icontxt, mpicomm )
    implicit none
    type(slinfo) :: sl
    real(8),intent(in)  :: Hsub(:,:)
    real(8),intent(out) :: H(:,:)
    integer,optional,intent(in) :: icontxt
    integer,optional,intent(in) :: mpicomm
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj

    H(:,:) = 0.0d0

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
      j1=min(j0+sl%nbsize-1,n)
      jj=j1-j0+1
      jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
      ip=-1
      is= 0
      do i0=1,m,sl%mbsize
        i1=min(i0+sl%mbsize-1,m)
        ii=i1-i0+1
        ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
        do j=1,jj
        do i=1,ii
          H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
        end do
        end do
        is = is + ii
      end do
      js = js + jj
    end do

    if ( present(icontxt) ) call dgsum2d( icontxt, 'All', ' ', m, n, H, m, -1, -1 )
    if ( present(mpicomm) ) call rsdft_allreduce( H, mpicomm )

  end subroutine d_gather_matrix

  subroutine z_gather_matrix( sl, Hsub, H, icontxt, mpicomm )
    implicit none
    type(slinfo) :: sl
    complex(8),intent(in) :: Hsub(:,:)
    complex(8),intent(out) :: H(:,:)
    integer,optional,intent(in) :: icontxt
    integer,optional,intent(in) :: mpicomm
    integer :: m,n,m0,n0,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj

    H(:,:) = (0.0d0,0.0d0)

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
      j1=min(j0+sl%nbsize-1,n)
      jj=j1-j0+1
      jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
      ip=-1
      is= 0
      do i0=1,m,sl%mbsize
        i1=min(i0+sl%mbsize-1,m)
        ii=i1-i0+1
        ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
        do j=1,jj
        do i=1,ii
          H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
        end do
        end do
        is = is + ii
      end do
      js = js + jj
    end do

    if ( present(icontxt) ) call zgsum2d( icontxt, 'All', ' ', m, n, H, m, -1, -1 )
    if ( present(mpicomm) ) call rsdft_allreduce( H, mpicomm )

  end subroutine z_gather_matrix


  SUBROUTINE d_gatherA_matrix( sl, Hsub, H, icontxt )
    implicit none
    type(slinfo) :: sl
    real(8),intent(IN) :: Hsub(:,:)
    real(8),intent(OUT) :: H(:,:)
    integer,optional,intent(IN) :: icontxt
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    include 'mpif.h'

    H(:,:) = 0.0d0

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
!            if ( i0+i-1 < j0+j-1 ) cycle
             H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do

    if ( present(icontxt) ) then
       call dgsum2d( icontxt, 'All', ' ', m, n, H, m, -1, -1 )
    else
       call rsdft_allreduce( H )
    end if

  END SUBROUTINE d_gatherA_matrix


  SUBROUTINE sl_block_map( m,n,sl )
    implicit none
    integer,intent(IN) :: m,n
    type(slinfo),intent(IN) :: sl
    integer :: i0,j0,i1,j1,ip,jp,ii,jj
    integer :: ib,jb,mblk,nblk
    integer,allocatable :: is(:,:),js(:,:)

    mblk = (m+sl%mbsize-1)/sl%mbsize
    nblk = (n+sl%nbsize-1)/sl%nbsize

    if ( allocated(blk_map) ) deallocate(blk_map)
    allocate( blk_map(mblk,nblk) )

    allocate( is(0:sl%nprow-1,0:sl%npcol-1) ) ; is=0
    allocate( js(0:sl%nprow-1,0:sl%npcol-1) ) ; js=0

    jp=-1
    jb= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jb=jb+1
       jp=mod(jp+1,sl%npcol)
       ip=-1
       is(:,jp)=0
       ib= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ib=ib+1
          ip=mod(ip+1,sl%nprow)
          blk_map(ib,jb)%iband(1)=i0
          blk_map(ib,jb)%iband(2)=i1
          blk_map(ib,jb)%jband(1)=j0
          blk_map(ib,jb)%jband(2)=j1
          blk_map(ib,jb)%istor(1)=is(ip,jp)+1
          blk_map(ib,jb)%istor(2)=is(ip,jp)+ii
          blk_map(ib,jb)%jstor(1)=js(ip,jp)+1
          blk_map(ib,jb)%jstor(2)=js(ip,jp)+jj
          blk_map(ib,jb)%iproc=ip
          blk_map(ib,jb)%jproc=jp
          is(ip,jp)=is(ip,jp)+ii
       end do
       js(:,jp)=js(:,jp)+jj
    end do
    deallocate( js, is )
  END SUBROUTINE sl_block_map


  subroutine backup_param_sl( nprow, npcol, mbsize, nbsize )
    implicit none
    integer,intent(in) :: nprow,npcol,mbsize,nbsize
    sl_backup%nprow =nprow
    sl_backup%npcol =npcol
    sl_backup%mbsize=mbsize
    sl_backup%nbsize=nbsize
  end subroutine backup_param_sl

  subroutine restore_param_sl( sl )
    implicit none
    type(slinfo),intent(inout) :: sl
    sl%nprow =sl_backup%nprow
    sl%npcol =sl_backup%npcol
    sl%mbsize=sl_backup%mbsize
    sl%nbsize=sl_backup%nbsize
  end subroutine restore_param_sl


  subroutine prep_param_sl( MB, np_grid, np_band )
    implicit none
    integer,intent(inout) :: MB
    integer,intent(in) :: np_grid, np_band
    integer :: i,j,n,NPCOL0,myrank,ierr
    integer :: NPROW,NPCOL,MBSIZE,NBSIZE
    logical :: disp_sw
    include 'mpif.h'

    call write_border( 0, " prep_param_sl(start)" )

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    call check_disp_switch( disp_sw, 0 )

! --- NPROW,NPCOL,MBSIZE,NBSIZE ---

    NPROW=0; NPCOL=0; MBSIZE=0; NBSIZE=0

    if ( NPROW<1 .or. NPCOL<1 ) then
      NPCOL0 = np_band*np_grid
      NPCOL  = np_band*np_grid
      NPROW  = 1
      do i=2,np_grid
        j=NPCOL0/i
        if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
          NPCOL=j
          NPROW=i
        end if
      end do
    else
      if ( NPROW*NPCOL > np_band*np_grid ) then
        write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL,np_band,np_grid,myrank
        stop
      end if
    end if

    if ( MBSIZE<1 .or. NBSIZE<1 ) then
      i=(MB+NPROW-1)/NPROW
      j=(MB+NPCOL-1)/NPCOL
      MBSIZE=min(i,j)
      NBSIZE=MBSIZE
    else
      if ( MBSIZE/=NBSIZE ) then
        write(*,*) "MBSIZE,NBSIZE=",MBSIZE,NBSIZE,myrank
        stop
      end if
    end if

    if ( disp_sw ) then
      write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
      write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
      n=max( NBSIZE*NPCOL, MB )
      n=min( n, (NBSIZE+1)*NPCOL )
      if ( disp_sw ) then
        write(*,*) "NBSIZE*NPCOL/=MB!"
        write(*,*) "recommended value for MB =",n
        write(*,*) "MB is replaced"
      end  if
      MB=n
    end if

    call backup_param_sl( NPCOL, NPROW, MBSIZE, NBSIZE )

    call write_border( 0, " prep_param_sl(end)" )

  end subroutine prep_param_sl


end module sl_tools_module
