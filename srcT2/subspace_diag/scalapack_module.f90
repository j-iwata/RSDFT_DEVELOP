MODULE scalapack_module

  use io_tools_module
  use vector_tools_module, only: vinfo
  use sl_tools_module, only: slinfo, gather_matrix, distribute_matrix, gatherA_matrix
#ifdef _EIGEN_
  use eigen_libs_mod
  use eigen_blacs_mod
#endif

  implicit none

  PRIVATE

  PUBLIC :: read_scalapack
  PUBLIC :: prep_scalapack
  PUBLIC :: init_scalapack
  PUBLIC :: d_gather_matrix_scalapack
  PUBLIC :: d_gatherA_matrix_scalapack
  PUBLIC :: d_scatter_matrix_scalapack
  PUBLIC :: trans_blkcy2cy
  PUBLIC :: trans_cy2blkcy
  PUBLIC :: fill_Matrix_scalapack
  PUBLIC :: end_eigen_free

  PUBLIC :: UPLO,NPROW,NPCOL,MBSIZE,NBSIZE,LLD_R,LLD_C &
           ,DESCA,DESCB,DESCZ,NP0,NQ0,NPX,NQX,usermap,sl

  PUBLIC :: DESCA_e, DESCB_e, DESCZ_e
  PUBLIC :: lda_r_e,lda_c_e
  PUBLIC :: LLD_R_e, LLD_C_e, MXLLD_e
  PUBLIC :: iflag_e, idiag, imate, irotv


  integer :: NPROW=0
  integer :: NPCOL=0
  integer :: MBSIZE=0
  integer :: NBSIZE=0
  integer :: LLD_R,LLD_C
  integer :: NP0,NQ0,NPX,NQX
  integer :: DESCA(9),DESCB(9),DESCZ(9)
  integer,allocatable :: usermap(:,:,:)
  character(1) :: UPLO

  integer :: DESCA_e(9),DESCB_e(9),DESCZ_e(9)
  integer :: lda_r_e,lda_c_e
  integer :: LLD_R_e, LLD_C_e, MXLLD_e
  integer,allocatable :: usermap2(:,:,:)
  integer :: myrow2, mycol2, nprow2, npcol2, ICTXT_E

  logical :: iblacs = .false.
  logical :: flag_read = .true.
  logical :: disp_sw

  integer :: node_partition(7)

  type(slinfo) :: sl
  integer :: icontxt

  logical :: iflag_e

  character(8) :: idiag="pdsyevd"
  integer :: imate=1
  integer :: irotv=2
  integer :: SW_USERMAP=0

CONTAINS


  SUBROUTINE init_scalapack( MB, np )
    implicit none
    integer,intent(INOUT) :: MB
    integer,intent(IN) :: np(4)
    integer :: NPCOL0,i,j,n,loop
!   logical :: disp_sw

    call write_border( 0, " init_scalapack(start)" )

!   MBSIZE = 0
!   NBSIZE = 0
    iblacs = .false.

    node_partition(:)=1
    node_partition(1:4)=np(1:4)

    if ( NPROW < 1 .or. NPCOL < 1 ) then

       NPCOL0 = np(1)*np(2)*np(3)*np(4)
       NPCOL  = NPCOL0
       NPROW  = 1
       do i=2,np(1)*np(2)*np(3)
          j=NPCOL0/i
          if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
             NPCOL=j
             NPROW=i
          end if
       end do

    else

       n=np(1)*np(2)*np(3)*np(4)
       if ( NPROW*NPCOL > n ) then
          write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL,np(4),np(1)*np(2)*np(3)
          call stop_program( "init_scalapack" )
       end if

    end if

    if ( MBSIZE < 1 .or. NBSIZE < 1 ) then

       i=(MB+NPROW-1)/NPROW
       j=(MB+NPCOL-1)/NPCOL
       MBSIZE=min(i,j)
       NBSIZE=MBSIZE

    else

       if ( MBSIZE /= NBSIZE ) then
          write(*,*) "MBSIZE,NBSIZE=",MBSIZE,NBSIZE
          write(*,*) "MBSIZE /= NBSIZE may not work well"
          call stop_program( "stop@init_scalapack" )
       end if

    end if

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

!   if ( NBSIZE*NPCOL /= MB ) then
!      write(*,*) "NBSIZE*NPCOL/=MB!"
!      n=max( NBSIZE*NPCOL, MB )
!      n=min( n, (NBSIZE+1)*NPCOL )
!      write(*,*) "recommended value for MB =",n
!      write(*,*) "replace MB"
!      MB=n
!      MBSIZE=0
!      NBSIZE=0
!   end if

    call write_border( 0, " init_scalapack(end)" )

  END SUBROUTINE init_scalapack


  SUBROUTINE read_scalapack
    implicit none
!   integer :: itmp(2)
    integer :: itmp(6)
    itmp(:)=-1
    call IOTools_readIntegerKeyword( "SCL", itmp )
    if ( itmp(1) > -1 ) NPROW=itmp(1)
    if ( itmp(2) > -1 ) NPCOL=itmp(2)
    if ( itmp(3) > -1 ) MBSIZE=itmp(3)
    if ( itmp(4) > -1 ) NBSIZE=itmp(4)
    if ( itmp(5) > -1 ) NPROW2=itmp(5)
    if ( itmp(6) > -1 ) NPCOL2=itmp(6)
    call IOTools_readStringKeyword(  "IDIAG", idiag )
    call IOTools_readIntegerKeyword( "IMATE", imate )
    call IOTools_readIntegerKeyword( "IROTV", irotv )
    call IOTools_readIntegerKeyword( "USERMAP", SW_USERMAP )
  END SUBROUTINE read_scalapack


  SUBROUTINE prep_scalapack( MB, v )
    implicit none
    integer,intent(INOUT) :: MB
    type(vinfo),intent(IN) :: v(2)
    integer :: ierr,NPCOL0,i,j,n,is,ik,m,ib,l,i1,i2,i3,i7
    integer :: MXLLD,MYROW,MYCOL,mm,mchk
    integer,save :: icount_visit=0, ICTXT=0, ICTXT0=0
    integer :: NUMROC,np_grid,np_band,myrnk_g,myrnk_b
!   logical :: disp_sw
    include 'mpif.h'

    integer :: ip_grid, ip_band

    integer :: eigen_comm
!   character(8) :: idiag
    integer :: MXLLD_e

    integer :: nnod, x_nnod, y_nnod
    integer :: inod, x_inod, y_inod

    integer,allocatable ::ranks1(:), ranks2(:)
    integer :: orig_group, group1, group2
    logical :: idisp=.false.
    integer :: ip, myrank
    integer :: comm1, comm2, comm_e
    integer :: nprocs, nprocs_e, nprocs_n

    if ( iblacs ) return

    call write_border( 0, " prep_scalapack(start)" )

    if ( icount_visit > 0 ) call blacs_gridexit(ICTXT)
    icount_visit=icount_visit+1

    iblacs = .true.

! set idiag & imate
!    call set_parameter

    np_grid = v(1)%pinfo%np
    np_band = v(2)%pinfo%np
    myrnk_g = v(1)%pinfo%me
    myrnk_b = v(2)%pinfo%me

! --- NPROW,NPCOL,MBSIZE,NBSIZE ---

    if ( NPROW < 1 .or. NPCOL < 1 ) then

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
          write(*,*) "NPROW,NPCOL,np_band,np_grid=",NPROW,NPCOL,np_band,np_grid
          call stop_program( "prep_scalapack(1)" )
       end if

    end if

    if ( MBSIZE < 1 .or. NBSIZE < 1 ) then

       i=(MB+NPROW-1)/NPROW
       j=(MB+NPCOL-1)/NPCOL
       MBSIZE=min(i,j)
       NBSIZE=MBSIZE

    else

       if ( MBSIZE/=NBSIZE ) then
          write(*,*) "MBSIZE,NBSIZE=",MBSIZE,NBSIZE
          call stop_program( "prep_scalapack(2)" )
       end if

    end if

!   call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
       write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
    end if

!   if ( NBSIZE*NPCOL /= MB ) then
!      n=max( NBSIZE*NPCOL, MB )
!      n=min( n, (NBSIZE+1)*NPCOL )
!      if ( disp_sw ) then
!         write(*,*) "NBSIZE*NPCOL/=MB!"
!         write(*,*) "recommended value for MB =",n
!         write(*,*) "MB is replaced"
!      end  if
!      MB=n
!   end if

! --- preparation for ScaLAPACK ---

    if ( .not.allocated(usermap) ) then

       allocate( usermap(0:NPROW-1,0:NPCOL-1,2) )
       usermap(:,:,:)=MPI_PROC_NULL

       select case( SW_USERMAP )
       case( 1 )

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
                mm=myrnk_g+np_grid*myrnk_b+1
!                if ( id_class(myrnk,5)==ik .and. &
!                     id_class(myrnk,6)==is .and. &
!                     id_class(myrnk,7)==i7 .and. mm > mchk ) then
!                   usermap(i,j,1)=n
!                   usermap(i,j,2)=l
!                end if
                if ( mm > mchk ) then
                   usermap(i,j,1)=n ! sequential # in MPI_COMM_WROLD
                   usermap(i,j,2)=l ! sequential # in comm_grid
                end if

             end do ! i1
             end do ! i2
             end do ! i3
          end do ! ib
       end do ! ik
       end do ! is
       end do ! i7

!----- Debug -----!
!   print original usermap 
!
       if ( disp_sw ) then
          write(*,*) "usermap(1)="
          do i1=0,NPROW-1
             write(*,*) usermap(i1,0:NPCOL-1,1)
          enddo
          write(*,*) "usermap(2)="
          do i1=0,NPROW-1
             write(*,*) usermap(i1,0:NPCOL-1,2)
          enddo
       end if
!----- Debug -----!

       case default ! SW_USERMAP

!----- re-define usermap for block-cyclic partitioning
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

                ip_grid = mod(m, np_grid)
                ip_band = m/np_grid

                i=mod(m+NPROW,NPROW)
                j=ip_band+np_band*int(mod(m,np_grid)/NPROW)

                usermap(i,j,1)=n
                usermap(i,j,2)=l

             end do ! i1
             end do ! i2
             end do ! i3
          end do ! ib
       end do ! ik
       end do ! is
       end do ! i7

!----- Debug -----!
!   print usermap for block-cyclic partitioning
!
       if ( disp_sw ) then
          write(*,*) "usermap(1)="
          do i1=0,NPROW-1
             write(*,*) usermap(i1,0:NPCOL-1,1)
          enddo
          write(*,*) "usermap(2)="
          do i1=0,NPROW-1
             write(*,*) usermap(i1,0:NPCOL-1,2)
          enddo
       end if
!----- Debug -----!

       end select

    end if

!----------------------------------------------------------------------
! create rank information for grouping 
!----------------------------------------------------------------------
      call mpi_comm_rank (mpi_comm_world, myrank, ierr )
      call mpi_comm_size( mpi_comm_world, nprocs, ierr )
      if (myrank==0) idisp=.true.

      if(nprow2==-1) nprow2=nprow   ! default  (if < 0)
      if(npcol2==-1) npcol2=npcol   ! default  (if < 0)

      nprocs_e=nprow2*npcol2        ! exec eigen_s
      nprocs_n=nprocs-nprocs_e      ! not exec eigen_s
      if (idisp) write(*,*) 'norocs_e,n=', nprocs_e, nprocs_n

      allocate(ranks1(nprocs_e), ranks2(nprocs_n))
      call create_ranks(ranks1,nprocs_e,ranks2,nprocs_n,nprocs, nprow2*npcol2)
      call mpi_barrier(mpi_comm_world,ierr)
   
      if (idisp) call show_dim1(ranks1, nprocs_e)
      if (idisp) call show_dim1(ranks2, nprocs_n)
      call mpi_barrier(mpi_comm_world,ierr)

!----------------------------------------------------------------------
! create new communicators
!----------------------------------------------------------------------
      call mpi_comm_group(mpi_comm_world, orig_group, ierr);
      call mpi_group_incl(orig_group, nprocs_e, ranks1, group1, ierr);
      call mpi_group_incl(orig_group, nprocs_n, ranks2, group2, ierr);
      call mpi_comm_create(mpi_comm_world, group1, comm1, ierr);
      call mpi_comm_create(mpi_comm_world, group2, comm2, ierr);

      ! set communicator for Eigen
      if ( myrank<nprocs_e ) then
         comm_e = comm1
         iflag_e=.true.            ! exec eigen
      else
         comm_e = comm2
         iflag_e=.false.           ! no exec eigen
      endif
#ifdef _EIGEN_
      ! eigen initialization
      call eigen_init( comm_e )
#endif
!----------------------------------------------------------------------
!    create descriptor for Block-cyclic distribution
!----------------------------------------------------------------------

!     ICTXT = MPI_COMM_WORLD
      call blacs_get(0, 0, ICTXT0)
      ICTXT = ICTXT0

      call blacs_gridmap(ICTXT,usermap(0,0,1),NPROW,NPROW,NPCOL)
      call blacs_gridinfo(ICTXT,NPROW,NPCOL,MYROW,MYCOL)

      LLD_R = NUMROC(MB,MBSIZE,MYROW,0,NPROW)
      LLD_C = NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)
!     MXLLD = LLD_R
      MXLLD = max(LLD_R,LLD_C)
      LLD_R=MXLLD
      LLD_C=MXLLD

      call descinit(DESCA,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
      call descinit(DESCB,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)
      call descinit(DESCZ,MB,MB,MBSIZE,NBSIZE,0,0,ICTXT,MXLLD,ierr)

      NP0=NUMROC(MB,MBSIZE,0,0,NPROW)
      NQ0=NUMROC(MB,NBSIZE,0,0,NPCOL)

      NPX=NUMROC(MB,MBSIZE,MYROW,0,NPROW)
      NQX=NUMROC(MB,NBSIZE,MYCOL,0,NPCOL)


!----------------------------------------------------------------------
!      create descriptor for cyclic-cyclic distribution (Eigen)
!----------------------------------------------------------------------
#ifdef _EIGEN_
      call eigen_get_procs(nnod, x_nnod, y_nnod)
      call eigen_get_id(inod, x_inod, y_inod)

      if (idisp) write(*,'(a9,i3)')         '#   nnod=', nnod
      if (idisp) write(*,'(a9,i3,x,a7,i3)') '# x_nnod=', x_nnod, 'y_nnod=',y_nnod

      if (x_nnod/=nprow2) then
         nprow2=x_nnod
         npcol2=y_nnod
         if (idisp) write(*,*) "# NPROW and NPCOL were updated by x_nnod and y_nnod)"
      endif

      ! create process mapping information (usermap for Eigen)
      call create_usermap2

       ! set context for eigen and map process

       !ICTXT_E = MPI_COMM_WORLD !eigen_get_blacs_context()
       ICTXT_E = eigen_get_blacs_context()
       call blacs_gridmap(ICTXT_E,usermap2(0,0,1),NPROW2,NPROW2,NPCOL2)

       LLD_R_e = NUMROC(MB,1,x_inod,0,nprow2)
       LLD_C_e = NUMROC(MB,1,y_inod,0,npcol2)

       ! get size of work array

       call eigen_get_matdims(MB, lda_r_e, lda_c_e)

       LLD_R_e=max(lda_r_e,LLD_R_e)
       LLD_C_e=max(lda_c_e,LLD_C_e)
       LLD_R_e=max(LLD_R_e,LLD_C_e)
       MXLLD_e = LLD_R_e

       if (iflag_e) then
          call descinit(DESCA_e,MB,MB,1,1,0,0,ICTXT_E,MXLLD_e,ierr)
          call descinit(DESCB_e,MB,MB,1,1,0,0,ICTXT_E,MXLLD_e,ierr)
          call descinit(DESCZ_e,MB,MB,1,1,0,0,ICTXT_E,MXLLD_e,ierr)
       else
          desca_e=-1
          descb_e=-1
          descz_e=-1
       endif
       if ( disp_sw ) then
          write(*,*) "LLD_R_e,LLD_C_e,MXLLD_e=",LLD_R_e,LLD_C_e,MXLLD_e
          write(*,*) "lda_r_e=",lda_r_e,lda_c_e
          write(*,*) "nnod   =",nnod
       end if
#endif
!-------------------------------------------------------------------------


    if ( disp_sw ) then
       write(*,*) "ICTXTX,ICTXT0    =",ICTXT,ICTXT0
       write(*,*) "NPROW,NPCOL      =",NPROW,NPCOL
       write(*,*) "MBSIZE,NBSIZE    =",MBSIZE,NBSIZE
       write(*,*) "MYROW,MYCOL      =",MYROW,MYCOL
       write(*,*) "LLD_R,LLD_C,MXLLD=",LLD_R,LLD_C,MXLLD
       write(*,*) "NP0,NQ0          =",NP0,NQ0
       write(*,*) "NPX,NQX          =",NPX,NQX
       write(*,*) "iblacs           =",iblacs
       write(*,*) "NPROW2, NPCOL2   =",NPROW2, NPCOL2
    end if

    icontxt   = ICTXT
    sl%myrow  = MYROW
    sl%mycol  = MYCOL
    sl%nprow  = NPROW
    sl%npcol  = NPCOL
    sl%mbsize = MBSIZE
    sl%nbsize = NBSIZE
    sl%icontxt_a = icontxt

    call write_border( 0, " prep_scalapack(end)" )

  END SUBROUTINE prep_scalapack


  SUBROUTINE d_gather_matrix_scalapack( Hsub, Htot )
    implicit none
    real(8),intent(IN)  :: Hsub(:,:)
    real(8),intent(OUT) :: Htot(:,:)
    call gather_matrix( sl, Hsub, Htot, icontxt )
  END SUBROUTINE d_gather_matrix_scalapack


  SUBROUTINE d_gatherA_matrix_scalapack( Hsub, Htot )
    implicit none
    real(8),intent(IN)  :: Hsub(:,:)
    real(8),intent(OUT) :: Htot(:,:)
    call gatherA_matrix( sl, Hsub, Htot, icontxt )
  END SUBROUTINE d_gatherA_matrix_scalapack


  SUBROUTINE d_scatter_matrix_scalapack( Vtot, Vsub )
    implicit none
    real(8),intent(IN)  :: Vtot(:,:)
    real(8),intent(OUT) :: Vsub(:,:)
    call distribute_matrix( sl, Vtot, Vsub )
  END SUBROUTINE d_scatter_matrix_scalapack

  SUBROUTINE trans_blkcy2cy(Hsub, Hsub_e, MB)
    implicit none
    real(8),intent(IN)     :: Hsub(:,:)
    real(8),intent(INOUT)  :: Hsub_e(:,:)
    integer,intent(IN) :: MB
    integer :: ICTXT

    ICTXT = icontxt
    ! form block-cyclic to cyclic
    call pdgemr2d(MB,MB,Hsub,1,1,DESCA,Hsub_e,1,1,DESCA_e,ICTXT)

  END SUBROUTINE trans_blkcy2cy

  SUBROUTINE trans_cy2blkcy(Vsub_e, Vsub, MB)
    implicit none
    real(8),intent(IN)     :: Vsub_e(:,:)
    real(8),intent(INOUT)  :: Vsub(:,:)
    integer,intent(IN) :: MB
    integer :: ICTXT

    ICTXT = icontxt
    ! form cyclic to block-cyclic
    call pdgemr2d(MB,MB,Vsub_e,1,1,DESCZ_e,Vsub,1,1,DESCZ,ICTXT)

  END SUBROUTINE trans_cy2blkcy

  SUBROUTINE fill_Matrix_scalapack(Hsub, MB)
    implicit none
    real(8),intent(INOUT)     :: Hsub(:,:)
    integer,intent(IN) :: MB
    include 'mpif.h'
    integer :: i, j, ii, jj, n, ms, me, ns, ne, iroot1, itag, ierr, nreq
    integer :: MBLK, nn, mm, i0, j0
    integer :: IPROW, IPCOL, myrank
    integer :: nsend_me, nrecv_me, nsize, icycle1
    integer :: i1, i2, j1, j2, m
    integer,allocatable :: istatus(:,:)
    integer,allocatable :: ireq(:)
    integer,allocatable :: irecv_me(:,:), isend_me(:,:)
    real(8),allocatable :: wtmp(:,:,:), Hsub_w(:,:)
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    integer :: nwork

! alloc work array
    MBLK=max(MBSIZE, NBSIZE)
!   allocate (wtmp(MBLK,MBLK,1000));wtmp=0.d0
    nwork=(MB-1)/MBLK+1
    nwork=((nwork-1)/nprow+1)*((nwork-1)/npcol+1)
    allocate (wtmp(MBLK,MBLK,nwork));wtmp=0.d0
    allocate (Hsub_w(MBLK,MBLK));Hsub_w=0.d0
!   allocate (ireq(1000));ireq=0
!   allocate (istatus(MPI_STATUS_SIZE,1000));istatus=0
!   allocate( irecv_me(99,0:8),isend_me(99,0:8) )
    allocate (ireq(nwork));ireq=0
    allocate (istatus(MPI_STATUS_SIZE,nwork));istatus=0
    allocate( irecv_me(nwork,0:8),isend_me(nwork,0:8) )

    nrecv_me      = 0
    nsend_me      = 0
    irecv_me(:,:) =-1
    isend_me(:,:) =-1


    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    j0=0
    do ns=1,MB,MBLK
       ne=min(ns+MBLK-1,MB)
       nn=ne-ns+1
       IPCOL  = mod( (ns-1)/MBLK, NPCOL )

       i0=0 
       do ms=1,MB,MBLK
          me = min( ms+MBLK-1, MB)
          mm = me-ms+1
          IPROW = mod( (ms-1)/MBLK, NPROW )

          iroot1 = usermap(IPROW,IPCOL,1)
          if (iroot1==myrank) then

             ! get source/destination id
             i=mod( (ns-1)/MBSIZE, NPROW )
             j=mod( (ms-1)/NBSIZE, NPCOL )
             n=usermap(i,j,1)

             if (ms>ns) then
             ! send
                nsend_me=nsend_me+1
                isend_me(nsend_me,0)=n
                isend_me(nsend_me,1)=ms
                isend_me(nsend_me,2)=me
                isend_me(nsend_me,3)=ns
                isend_me(nsend_me,4)=ne
                isend_me(nsend_me,5)=i0+1
                isend_me(nsend_me,6)=i0+mm
                isend_me(nsend_me,7)=j0+1
                isend_me(nsend_me,8)=j0+nn
                i0=i0+mm

             else if (ms<ns) then
             ! recv
                nrecv_me=nrecv_me+1
                irecv_me(nrecv_me,0)=n
                irecv_me(nrecv_me,1)=ms
                irecv_me(nrecv_me,2)=me
                irecv_me(nrecv_me,3)=ns
                irecv_me(nrecv_me,4)=ne
                irecv_me(nrecv_me,5)=i0+1
                irecv_me(nrecv_me,6)=i0+mm
                irecv_me(nrecv_me,7)=j0+1
                irecv_me(nrecv_me,8)=j0+nn
                i0=i0+mm

             else if (ms==ns) then
             ! copy
!               Hsub_w(1:mm,1:nn) = Hsub(i0:i0+mm-1,j0:j0+nn-1)
                Hsub_w(1:mm,1:nn) = Hsub(i0+1:i0+mm,j0+1:j0+nn)
                do jj=1,nn
                  do ii=jj+1,mm
!                    Hsub (i0+jj-1,j0+ii-1) = Hsub_w(ii,jj)
                     Hsub (i0+jj,j0+ii) = Hsub_w(ii,jj)
                  enddo
                enddo
                i0=i0+mm
             else
             endif

          endif  !(iroot1==myrank)

       enddo
       if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn
 
    enddo

!  write(*,*) "debug> myrank,nrecv_me=",myrank,nrecv_me
!  write(*,*) "debug> myrank,nsend_me=",myrank,nsend_me

    n=0
    do i=1,nsend_me
       m=(isend_me(i,2)-isend_me(i,1)+1)*(isend_me(i,4)-isend_me(i,3)+1)
       n=max(m,n)
    end do
    allocate( vtmp2(n,nsend_me) )

    n=0
    do i=1,nrecv_me
       m=(irecv_me(i,2)-irecv_me(i,1)+1)*(irecv_me(i,4)-irecv_me(i,3)+1)
       n=max(m,n)
    end do
    allocate( wtmp2(n,nrecv_me) )

    if ( nsend_me>0 .or. nrecv_me>0 ) then
       nreq=0
       do i=1,nsend_me
          n =isend_me(i,0)
          ms=isend_me(i,1)
          me=isend_me(i,2)
          ns=isend_me(i,3)
          ne=isend_me(i,4)
          i1=isend_me(i,5)
          i2=isend_me(i,6)
          j1=isend_me(i,7)
          j2=isend_me(i,8)
          j=0
          do i0=i1,i2
          do j0=j1,j2
             j=j+1
             vtmp2(j,i)=Hsub(i0,j0)
          end do
          end do
          itag=ns+10*ne+100*ms+1000*me
          nreq=nreq+1
          call mpi_isend(vtmp2(1,i),j,MPI_REAL8,n,itag &
               ,mpi_comm_world,ireq(nreq),ierr)
       end do


       do i=1,nrecv_me
          n =irecv_me(i,0)
          ms=irecv_me(i,1)
          me=irecv_me(i,2)
          ns=irecv_me(i,3)
          ne=irecv_me(i,4)
          i1=irecv_me(i,5)
          i2=irecv_me(i,6)
          j1=irecv_me(i,7)
          j2=irecv_me(i,8)
          j=(me-ms+1)*(ne-ns+1)
          itag=ms+10*me+100*ns+1000*ne
          nreq=nreq+1
          call mpi_irecv(wtmp2(1,i),j,MPI_REAL8,n,itag &
               ,mpi_comm_world,ireq(nreq),ierr)
       end do
       call mpi_waitall(nreq,ireq,istatus,ierr)


       do i=1,nrecv_me
          i1=irecv_me(i,5)
          i2=irecv_me(i,6)
          j1=irecv_me(i,7)
          j2=irecv_me(i,8)
          j=0
          do j0=j1,j2
             do i0=i1,i2
                j=j+1
                Hsub(i0,j0)=wtmp2(j,i)
             end do
          end do
       end do

    end if

   write(*,*) "debug> myrank,nrecv_me=",myrank,nrecv_me,nreq
   write(*,*) "debug> myrank,nsend_me=",myrank,nsend_me,nreq

    call mpi_barrier(mpi_comm_world,ierr)
    deallocate (ireq, istatus, wtmp, Hsub_w)
    deallocate( irecv_me, isend_me )


  END SUBROUTINE fill_Matrix_scalapack

  subroutine create_ranks(a1,n1,a2,n2,np,ncon)
      implicit none
      integer,intent(inout) :: a1(*), a2(*)
      integer,intent(in) :: n1, n2 , np, ncon
      integer :: i, i0, j0

      i0=0; j0=0
      do i=1,np
         if (i<=ncon) then
            i0=i0+1
            a1(i0)=i-1
         else
            j0=j0+1
            a2(j0)=i-1
         endif
      enddo

      return
  end subroutine create_ranks

  subroutine show_dim1(a, n)
      implicit none
      integer,intent(in) :: a(n), n
      character(13) :: nfmt      ! '(ni3)' n:1-99999

      if (n<=0) return
      write(nfmt,'("("i0"i3)")') n
      write(*,nfmt) a
      return
  end subroutine show_dim1

  subroutine create_usermap2
      implicit none
      include 'mpif.h'
      integer :: i, j, ip

      ! create process mapping information (usermap for Eigen)
      if ( .not.allocated(usermap2) ) then
         allocate( usermap2(0:NPROW2-1,0:NPCOL2-1,2) )
         usermap2(:,:,:)=MPI_PROC_NULL
         ip=-1
         do j=1,npcol2
            do i=1,nprow2
               ip=ip+1
               usermap2(i-1,j-1,1)=ip
            enddo
         enddo
      endif

      ! print usermap for Eigen partitioning
      if ( disp_sw ) then
         write(*,*) "usermap2(1)="
         do i=0,NPROW2-1
            write(*,*) usermap2(i,0:NPCOL2-1,1)
         enddo
      end if

      return
  end subroutine create_usermap2

  subroutine end_eigen_free
#ifdef _EIGEN_
     call eigen_free()
#endif
  end subroutine end_eigen_free

  subroutine set_parameter
!   idiag = "eigen_s "
    idiag = "pdsyevd"
!   idiag = "pzheevd"
!
!   imate = 0  ! tri-matirix (Original)
    imate = 1  ! tri-matirix (Original)
!   imate = 2  ! full-matrix
!   imate = 3  ! tri-matrix + symetric exchange

!   irotv = 0
!   irotv = 1
    irotv = 2
  end subroutine set_parameter

END MODULE scalapack_module
