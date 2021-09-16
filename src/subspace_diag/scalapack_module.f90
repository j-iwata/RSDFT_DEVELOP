module scalapack_module

  use parallel_module, only: node_partition, myrank, disp_switch_parallel, &
                             np_grid, np_band, myrank_g, myrank_b, &
                             myrank_k, myrank_s, id_class
  use io_tools_module, only: IOTools_readIntegerKeyword

  implicit none

  private
  public :: prep_scalapack
  public :: init_scalapack
  public :: NPROW,NPCOL,MBSIZE,NBSIZE,LLD_R,LLD_C,DESCA,DESCB,DESCZ, &
            NP0,NQ0,NPX,NQX,UPLO,usermap
  public :: read_scalapack
  public :: allocated_workarray_scalapack

  integer :: NPROW=0
  integer :: NPCOL=0
  integer :: MBSIZE=0
  integer :: NBSIZE=0
  integer :: LLD_R,LLD_C
  integer :: NP0,NQ0,NPX,NQX
  integer :: DESCA(9),DESCB(9),DESCZ(9)
  integer,allocatable :: usermap(:,:,:)
  character(1) :: UPLO

  logical :: iblacs = .false.
  logical :: flag_read = .true.
  logical :: is_workarray_allocated = .false.

  integer :: ialgo_sl2 = 0

contains


  subroutine init_scalapack( MB )
    implicit none
    integer,intent(inout) :: MB
    integer :: NPCOL0,i,j,n,loop

    if ( ialgo_sl2 /= 0 ) then
      call init_scalapack_2( MB, node_partition(4), product(node_partition(1:3)) )
      return
    end if

    call write_border( 0, " init_scalapack(start)" )

    MBSIZE = 0
    NBSIZE = 0
    iblacs = .false.

    if ( NPROW<1 .or. NPCOL<1 ) then
       NPCOL0 = node_partition(1)*node_partition(2) &
               *node_partition(3)*node_partition(4)
       NPCOL  = NPCOL0
       NPROW  = 1
       do i=2,node_partition(1)*node_partition(2)*node_partition(3)
          j=NPCOL0/i
          if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
             NPCOL=j
             NPROW=i
          end if
       end do
    else
       n=node_partition(1)*node_partition(2)*node_partition(3)*node_partition(4)
       if ( NPROW*NPCOL > n ) then
          write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL &
               ,node_partition(4) &
               ,node_partition(1)*node_partition(2)*node_partition(3),myrank
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
          write(*,*) "MBSIZE /= NBSIZE may not work well"
          stop "stop@init_scalapack"
       end if
    end if

    if ( disp_switch_parallel ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       if ( disp_switch_parallel ) then
         write(*,*) "NBSIZE*NPCOL/=MB!"
         write(*,*) "recommended value for MB =",n
         write(*,*) "replace MB"
       end if
       MB=n
       MBSIZE=0
       NBSIZE=0
    end if

    is_workarray_allocated = .false.

    call write_border( 0, " init_scalapack(end)" )

  end subroutine init_scalapack


  subroutine init_scalapack_2( nband_in, nprocs_band, nprocs_grid )
    use gcd_lcm_pf_module, only: prime_factorization, lcm
    use sl_variables, only: sl0, sl1, sl2
    use pinfo_module, only: get_world_rank_pinfo
    implicit none
    integer,intent(in) :: nband_in
    integer,intent(in) :: nprocs_band, nprocs_grid
    integer,allocatable :: ifact(:), lst(:)
    integer :: i,j,n,m,grp(2),grp_r(2),sqrt_nprocs_band,itmp(4)
    integer :: num_blocks,block_size,nband,idiv_switch,LDR,LDC
    integer,external :: NUMROC
    logical :: disp_on

    call flush_barrier(6)
    call write_border( 0, ' init_scalapack_2(start)' )
    call check_disp_switch( disp_on, 0 )

    nband = nband_in

    if ( disp_on ) then
      write(*,*) "nband=",nband
      write(*,*) "nprocs_band=",nprocs_band
      write(*,*) "nprocs_grid=",nprocs_grid
    end if

    do n = 1, nprocs_band
      if ( n*n <= nprocs_band ) then
        m = n
      else
        exit
      end if
    end do
    sqrt_nprocs_band = m
    if ( disp_on ) then
      write(*,'("The largest n*n process grid (within nprocs_band): " &
                ,i3," x ",i3," = ",i4)') m,m,m*m
    end if

    allocate( ifact(nprocs_band) ); ifact=0
    call prime_factorization( nprocs_band, ifact )
    j=1
    do i = 1, nprocs_band
      if ( ifact(i) > 0 ) then
        j=j*i**ifact(i)
        ! if ( disp_on ) write(*,*) i, ifact(i), j
      end if
    end do
    !
    n=sum(ifact)
    allocate( lst(n) ); lst=0
    m=0
    do i = 1, size(ifact)
      do j = 1, ifact(i)
        m = m + 1
        lst(m) = i
      end do
    end do
    grp = (/ 1, 1 /)
    grp_r = (/ 1, nprocs_band /)
    call div_2products( lst, grp, grp_r )
    num_blocks = lcm( grp_r(1), grp_r(2) )
    block_size = nband / num_blocks
    if ( disp_on ) then
      write(*,*) "grp_r(1:2),product(grp_r)",grp_r,product(grp_r)
      write(*,*) "Number of blocks:",num_blocks
      write(*,*) "Block size:",block_size,mod(nband,num_blocks)
    end if
    !
    if ( allocated(sl0%usermap) ) deallocate(sl0%usermap)
    sl0%distribution_method = ialgo_sl2
    !
    select case( sl0%distribution_method )
    case( 1 )
      sl0%nprow  = sqrt_nprocs_band
      sl0%npcol  = sl0%nprow
      sl0%mbsize = nband / sqrt_nprocs_band
      sl0%nbsize = sl0%mbsize
      allocate( sl0%usermap(0:sl0%nprow-1,0:sl0%npcol-1) ); sl0%usermap=-1
      n=-1
      do j = 0, sl0%npcol-1
        do i = 0, sl0%nprow-1
          n=n+1
          sl0%usermap(i,j) = get_world_rank_pinfo( 0, n, myrank_k, myrank_s )
          ! write(*,*) n,i,j,sl0%usermap(i,j)
        end do
      end do
    case( 2 )
      sl0%nprow  = grp_r(1)
      sl0%npcol  = grp_r(2)
      sl0%mbsize = block_size
      sl0%nbsize = sl0%mbsize
      allocate( sl0%usermap(0:sl0%nprow-1,0:sl0%npcol-1) ); sl0%usermap=-1
      n=-1
      do j = 0, sl0%npcol-1
        do i = 0, sl0%nprow-1
          n=n+1
          sl0%usermap(i,j) = get_world_rank_pinfo( 0, n, myrank_k, myrank_s )
          ! write(*,*) n,i,j,sl0%usermap(i,j)
        end do
      end do
    case( 3 )
      sl0%nprow  = grp_r(1)
      sl0%npcol  = grp_r(2)
      sl0%mbsize = nband / maxval( grp_r )
      sl0%nbsize = sl0%mbsize
      allocate( sl0%usermap(0:sl0%nprow-1,0:sl0%npcol-1) ); sl0%usermap=-1
      n=-1
      do j = 0, sl0%npcol-1
        do i = 0, sl0%nprow-1
          n=n+1
          sl0%usermap(i,j) = get_world_rank_pinfo( 0, n, myrank_k, myrank_s )
          ! write(*,*) n,i,j,sl0%usermap(i,j)
        end do
      end do
    case( 4 )
      sl0%nprow  = min( nprocs_grid, nband )
      sl0%npcol  = min( nprocs_band, nband )
      num_blocks = lcm( sl0%nprow, sl0%npcol )
      ! sl0%mbsize = nband / max( sl0%nprow, sl0%npcol )
      sl0%mbsize = nband / num_blocks
      sl0%nbsize = sl0%mbsize
      allocate( sl0%usermap(0:sl0%nprow-1,0:sl0%npcol-1) ); sl0%usermap=-1
      n=-1
      do j = 0, sl0%npcol-1
        do i = 0, sl0%nprow-1
          n=n+1
          sl0%usermap(i,j) = get_world_rank_pinfo( i, j, myrank_k, myrank_s )
          ! write(*,*) n,i,j,sl0%usermap(i,j)
        end do
      end do
    end select

    sl0%nband = nband

    if ( disp_on ) then
      write(*,'("distribution method id =",i6)') sl0%distribution_method
      write(*,'("nprow , npcol ",2i6)') sl0%nprow , sl0%npcol
      write(*,'("mbsize, nbsize",2i6)') sl0%mbsize, sl0%nbsize
    end if

    deallocate( lst )
    deallocate( ifact )

    call blacs_get( 0, 0, sl0%icontxt_sys )
    sl0%icontxt_a = sl0%icontxt_sys
    call blacs_gridmap(sl0%icontxt_a,sl0%usermap,sl0%nprow,sl0%nprow,sl0%npcol)
    call blacs_gridinfo(sl0%icontxt_a,m,n,sl0%myrow,sl0%mycol)
    write(*,*) "myrow, mycol",sl0%myrow, sl0%mycol, myrank
    !
    ! ----------
    !
    sl1%nband = sl0%nband
    itmp(1:4) = (/ sl0%nprow, sl0%npcol, sl0%mbsize, sl0%nbsize /)
    call IOTools_readIntegerKeyword( 'SCL3', itmp )
    sl1%nprow = itmp(1)
    sl1%npcol = itmp(2)
    sl1%mbsize = itmp(3)
    sl1%nbsize = sl1%mbsize
    if ( disp_on ) then
      write(*,'("ScaLapack parameters")')
      write(*,'("nprow1 , npcol1 ",2i6)') sl1%nprow , sl1%npcol
      write(*,'("mbsize1, nbsize1",2i6)') sl1%mbsize, sl1%nbsize
    end if
    if ( allocated(sl1%usermap) ) deallocate( sl1%usermap )
    allocate( sl1%usermap(0:sl1%nprow-1,0:sl1%npcol-1) ); sl1%usermap=-1
    n=-1
    do j = 0, sl1%npcol-1
      do i = 0, sl1%nprow-1
        n=n+1
        sl1%usermap(i,j) = n
        ! sl1%usermap(i,j) = get_world_rank_pinfo( i, j, myrank_k, myrank_s )
        ! write(*,*) n,i,j,sl1%usermap(i,j)
      end do
    end do
    call blacs_get( 0, 0, sl1%icontxt_sys )
    sl1%icontxt_a = sl1%icontxt_sys
    call blacs_gridmap(sl1%icontxt_a,sl1%usermap,sl1%nprow,sl1%nprow,sl1%npcol)
    sl1%lwork=0
    sl1%lrwork=0
    sl1%liwork=0
    sl1%idiag = ''
    sl1%uplo = ''
    call flush_barrier(6)
    if ( any(sl1%usermap == myrank) ) then
      call blacs_gridinfo(sl1%icontxt_a,m,n,sl1%myrow,sl1%mycol)
      write(*,*) "myrow1, mycol1",sl1%myrow, sl1%mycol, m,n,myrank
      LDR = NUMROC( nband, sl1%mbsize, sl1%myrow, 0, sl1%nprow )
      LDC = NUMROC( nband, sl1%nbsize, sl1%mycol, 0, sl1%npcol )
      call descinit( sl1%desca,nband,nband,sl1%mbsize,sl1%nbsize,0,0,sl1%icontxt_a,LDR,i )
      call descinit( sl1%descz,nband,nband,sl1%mbsize,sl1%nbsize,0,0,sl1%icontxt_a,LDR,i )
      sl1%ldr = LDR
      sl1%ldc = LDC
      sl1%idiag = 'PZHEEVD'
      sl1%uplo = 'L'
    end if
    !
    ! ----------
    !
    sl2%nband  = sl0%nband
    sl2%npcol  = nprocs_band
    sl2%nprow  = min( nprocs_grid, sl2%nband )
    sl2%nbsize = ( sl2%nband + sl2%npcol - 1 ) / sl2%npcol
    sl2%mbsize = ( sl2%nband + sl2%nprow - 1 ) / sl2%nprow
    if ( allocated(sl2%usermap) ) deallocate(sl2%usermap)
    allocate( sl2%usermap(0:sl2%nprow-1,0:sl2%npcol-1) ); sl2%usermap=-1
    do j = 0, sl2%npcol-1
      do i = 0, sl2%nprow-1
        sl2%usermap(i,j) = get_world_rank_pinfo( i, j, myrank_k, myrank_s )
      end do
    end do
    call blacs_get( 0, 0, sl2%icontxt_sys )
    sl2%icontxt_a = sl2%icontxt_sys
    call blacs_gridmap(sl2%icontxt_a,sl2%usermap,sl2%nprow,sl2%nprow,sl2%npcol)
    sl2%lwork=0
    sl2%lrwork=0
    sl2%liwork=0
    sl2%idiag = ''
    sl2%uplo = ''
    call flush_barrier(6)
    if ( any(sl2%usermap == myrank) ) then
      call blacs_gridinfo(sl2%icontxt_a,m,n,sl2%myrow,sl2%mycol)
      write(*,*) "myrow2, mycol2",sl2%myrow, sl2%mycol, m,n,myrank
      LDR = NUMROC( nband, sl2%mbsize, sl2%myrow, 0, sl2%nprow )
      LDC = NUMROC( nband, sl2%nbsize, sl2%mycol, 0, sl2%npcol )
      sl2%ldr = LDR
      sl2%ldc = LDC
    else
      sl2%myrow = -1
      sl2%mycol = -1
    end if

    call flush_barrier(6)
    call write_border( 0, ' init_scalapack_2(end)' )

  contains

    recursive subroutine div_2products( lst, grp, grp_r )
      implicit none
      integer,intent(in) :: lst(:)
      integer,intent(inout) :: grp(:), grp_r(:)
      integer :: n,i,idiff_r,idiff_t
      n = size(lst)
      do i = 1, 2
        grp(i) = grp(i)*lst(n)
        if ( n > 2 ) then
          call div_2products( lst(1:n-1), grp, grp_r )
        else if ( n == 2 ) then
          idiff_t = abs( grp(1) - grp(2) )
          idiff_r = abs( grp_r(1) - grp_r(2) )
          if ( idiff_t < idiff_r ) grp_r=grp
        end if
        grp(i) = grp(i)/lst(n)
      end do
    end subroutine div_2products

  end subroutine init_scalapack_2


  subroutine read_scalapack
    implicit none
    integer :: itmp(2)
    itmp(:)=-1
    call IOTools_readIntegerKeyword( "SCL", itmp )
    if ( itmp(1) > -1 ) NPROW=itmp(1)
    if ( itmp(2) > -1 ) NPCOL=itmp(2)
    itmp(1)=-1
    call IOTools_readIntegerKeyword( "SCL2", itmp(1) )
    if ( itmp(1) > -1 ) ialgo_sl2 = itmp(1)
  end subroutine read_scalapack


  subroutine prep_scalapack( MB )
    implicit none
    integer,intent(inout) :: MB
    integer :: ierr,NPCOL0,i,j,n,is,ik,m,ib,l,i1,i2,i3,i7
    integer :: MXLLD,MYROW,MYCOL,mm,mchk
    integer,save :: icount_visit=0, ICTXT=0, ICTXT0=0
    integer :: NUMROC
    include 'mpif.h'

#ifndef _LAPACK_

    if ( iblacs ) return

    call write_border( 0, " prep_scalapack(start)" )

    if ( icount_visit>0 ) then
       call blacs_gridexit(ICTXT)
    end if
    icount_visit=icount_visit+1

    iblacs = .true.

! --- NPROW,NPCOL,MBSIZE,NBSIZE ---

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

    if ( disp_switch_parallel ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       if ( disp_switch_parallel ) then
          write(*,*) "NBSIZE*NPCOL/=MB!"
          write(*,*) "recommended value for MB =",n
          write(*,*) "MB is replaced"
       end  if
       MB=n
    end if

! --- preparation for ScaLAPACK ---

    if ( .not.allocated(usermap) ) then

       allocate( usermap(0:NPROW-1,0:NPCOL-1,2) )
       usermap(:,:,:)=MPI_PROC_NULL

       n=-1
       do i7=0,node_partition(7)-1
       do is=0,node_partition(6)-1
       do ik=0,node_partition(5)-1
          m=-1 ; mchk=-NPROW*NPCOL
          do ib=0,node_partition(4)-1
             l=-1
             do i3=0,node_partition(3)-1
             do i2=0,node_partition(2)-1
             do i1=0,node_partition(1)-1
                n=n+1
                m=m+1 ; if ( mod(m,NPROW*NPCOL) == 0 ) mchk=mchk+NPROW*NPCOL
                l=l+1
                i=mod(m+NPROW,NPROW)
                j=m/NPROW
                mm=myrank_g+np_grid*myrank_b+1
                if ( id_class(myrank,5)==ik .and. &
                     id_class(myrank,6)==is .and. &
                     id_class(myrank,7)==i7 .and. mm > mchk ) then
                   usermap(i,j,1)=n
                   usermap(i,j,2)=l
                end if
             end do ! i1
             end do ! i2
             end do ! i3
          end do ! ib
       end do ! ik
       end do ! is
       end do ! i7

    end if

    call blacs_get(0,0,ICTXT0)
    ICTXT = ICTXT0

    call blacs_gridmap(ICTXT,usermap(0,0,1),NPROW,NPROW,NPCOL)
    call blacs_gridinfo(ICTXT,NPROW,NPCOL,MYROW,MYCOL)

    LLD_R = NUMROC(MB,MBSIZE,MYROW,0,NPROW)
    LLD_C = NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)
    MXLLD = LLD_R

    call descinit(DESCA,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
    call descinit(DESCB,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
    call descinit(DESCZ,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)

    NP0=NUMROC(MB,MBSIZE,0,0,NPROW)
    NQ0=NUMROC(MB,NBSIZE,0,0,NPCOL)

    NPX=NUMROC(MB,MBSIZE,MYROW,0,NPROW)
    NQX=NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)

    if ( disp_switch_parallel ) then
       write(*,*) "ICTXTX,ICTXT0    =",ICTXT,ICTXT0
       write(*,*) "NPROW,NPCOL      =",NPROW,NPCOL
       write(*,*) "MBSIZE,NBSIZE    =",MBSIZE,NBSIZE
       write(*,*) "MYROW,MYCOL      =",MYROW,MYCOL
       write(*,*) "LLD_R,LLD_C,MXLLD=",LLD_R,LLD_C,MXLLD
       write(*,*) "NP0,NQ0          =",NP0,NQ0
       write(*,*) "NPX,NQX          =",NPX,NQX
       write(*,*) "iblacs           =",iblacs
    end if

    call write_border( 0, " prep_scalapack(end)" )

#endif

  end subroutine prep_scalapack


  logical function allocated_workarray_scalapack( flag )
    implicit none
    logical,optional,intent(in) :: flag
    if ( present(flag) ) is_workarray_allocated = flag
    allocated_workarray_scalapack = is_workarray_allocated
  end function allocated_workarray_scalapack

end module scalapack_module
