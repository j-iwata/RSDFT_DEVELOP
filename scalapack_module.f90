MODULE scalapack_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: prep_scalapack,prep_0_scalapack, UPLO &
           ,NPROW,NPCOL,MBSIZE,NBSIZE,LLD_R,LLD_C,DESCA,DESCB,DESCZ &
           ,NP0,NQ0,NPX,NQX,usermap

  integer :: NPROW=0,NPCOL=0
  integer :: MBSIZE=0,NBSIZE=0
  integer :: LLD_R,LLD_C
  integer :: NP0,NQ0,NPX,NQX
  integer :: DESCA(9),DESCB(9),DESCZ(9)
  integer,allocatable :: usermap(:,:,:)
  character(1) :: UPLO

  logical :: iblacs = .false.
  logical :: flag_read = .true.

CONTAINS

  SUBROUTINE read_scalapack(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: ierr,i
    character(3) :: cbuf,ckey
    NPROW=0
    NPCOL=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:3) == "SCL" ) then
             backspace(unit)
             read(unit,*) cbuf,NPROW,NPCOL
          end if
       end do
999    continue
       write(*,*) "NPROW,NPCOL=",NPROW,NPCOL
    end if
    call mpi_bcast(NPROW,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NPCOL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE read_scalapack

  SUBROUTINE prep_scalapack(MB,disp_switch)
    integer,intent(INOUT) :: MB
    logical,intent(IN) :: disp_switch
    integer :: ierr,NPCOL0,i,j,n,is,ik,m,ib,l,i1,i2,i3
    integer :: MXLLD,MYROW,MYCOL,mm,mchk
    integer,save :: icount_visit=0, ICTXT=0, ICTXT0=0
    integer :: NUMROC

    if ( iblacs ) return

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

    if ( DISP_SWITCH ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       if ( disp_switch ) then
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
!             j=mod(m+NPCOL,NPCOL)
!!             i=m/NPCOL
!             i=mod(m/NPCOL,NPROW)
             i=mod(m+NPROW,NPROW)
             j=m/NPROW
             mm=myrank_g+nprocs_g*myrank_b+1
             if ( id_class(myrank,5)==ik .and. id_class(myrank,6)==is .and. &
                  mm > mchk ) then
                usermap(i,j,1)=n
                usermap(i,j,2)=l
             end if
          end do
          end do
          end do
       end do
    end do
    end do

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

    if ( DISP_SWITCH ) then
       write(*,*) "ICTXTX,ICTXT0    =",ICTXT,ICTXT0
       write(*,*) "NPROW,NPCOL      =",NPROW,NPCOL
       write(*,*) "MBSIZE,NBSIZE    =",MBSIZE,NBSIZE
       write(*,*) "MYROW,MYCOL      =",MYROW,MYCOL
       write(*,*) "LLD_R,LLD_C,MXLLD=",LLD_R,LLD_C,MXLLD
       write(*,*) "NP0,NQ0          =",NP0,NQ0
       write(*,*) "NPX,NQX          =",NPX,NQX
       write(*,*) "iblacs           =",iblacs
    end if

  END SUBROUTINE prep_scalapack


  SUBROUTINE prep_0_scalapack(MB,disp_switch)
    integer,intent(INOUT) :: MB
    logical,intent(IN) :: disp_switch
    integer :: NPCOL0,i,j,n

    MBSIZE = 0
    NBSIZE = 0
    iblacs = .false.

    if ( flag_read ) then
       NPROW = 0
       NPCOL = 0
       call read_scalapack(myrank,1)
       flag_read = .false.
    end if

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
          stop
       end if
    end if

    if ( DISP_SWITCH ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

    if ( NBSIZE*NPCOL/=MB ) then
       write(*,*) "NBSIZE*NPCOL/=MB!"
       n=max( NBSIZE*NPCOL, MB )
       n=min( n, (NBSIZE+1)*NPCOL )
       write(*,*) "recommended value for MB =",n
       write(*,*) "replace MB"
       MB=n
    end if

  END SUBROUTINE prep_0_scalapack

END MODULE scalapack_module
