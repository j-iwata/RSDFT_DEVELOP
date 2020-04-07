MODULE subspace_mate_sl_module

  use rgrid_module, only: dV,zdV
  use parallel_module
  use hamiltonian_module
  use wf_module, only: unk,esp,hunk,iflag_hunk,gather_b_wf
  use scalapack_module
  use subspace_diag_variables, only: MB_diag,NBLK1,zero,TYPE_MAIN,Hsub,mat_block,nLB
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use watch_module
  use vector_tools_module, only: vinfo
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none

  PRIVATE
  PUBLIC :: subspace_mate_sl
  PUBLIC :: subspace_mate_sl_bp1
  PUBLIC :: subspace_mate_sl_bp2
  PUBLIC :: subspace_mate_sl_bp3

#ifdef _DRSDFT_
  real(8),allocatable :: hpsi(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='T', TRANSB='N'
#else
  complex(8),allocatable :: hpsi(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='C', TRANSB='N'
#endif

CONTAINS


  SUBROUTINE subspace_mate_sl(k,s)
    integer,intent(IN) :: k,s

    integer :: i,i0,i1,i2,ib1,ib2,j,j0,j1,j2,m,me,mm,ms,MB
    integer :: MBLK,MBLKH,ML0,n1,n2,mms,mme,nnn,nns,n,ne,nn,ns
!   integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer,allocatable :: istatus(:,:),ireq(:)
    integer :: nreq,itag,ierr,nwork
    integer :: IPROW,IPCOL,iroot1,iroot2,mrnk,nrecv_me,nsend_me
    integer,allocatable :: ir(:),id(:),irecv_me(:,:),isend_me(:,:)
    complex(8) :: ztmp
    type(time) :: t

    call write_border( 1, " subspace_mate_sl(start)" )
    call start_timer( t )

    UPLO = 'L'

    n1    = ML_0
    n2    = ML_1
    ML0   = ML_1-ML_0+1
    mrnk  = id_class(myrank,4)
    MB    = MB_diag

    NBLK1 = 4

    MBLK  = min( MBSIZE,NBSIZE )

    call gather_b_wf( k, s )

!   allocate( irecv_me(99,0:8),isend_me(99,0:8) )
    nwork=(MB-1)/MBLK+1
    nwork=((nwork-1)/nprow+1)*((nwork-1)/npcol+1)
    allocate( istatus(MPI_STATUS_SIZE,nwork), ireq(nwork) )
    allocate( irecv_me(nwork,0:8),isend_me(nwork,0:8) )

    nrecv_me      = 0
    nsend_me      = 0
    irecv_me(:,:) =-1
    isend_me(:,:) =-1

    allocate( hpsi(n1:n2,MBLK) )
    hpsi=zero

    j0=0

    do ns=MB_0,MB_1,MBLK

       ne=min(ns+MBLK-1,MB_1)
       nn=ne-ns+1

       IPCOL = mod( (ns-1)/NBSIZE, NPCOL )

       if ( iflag_hunk >= 1 ) then
          hpsi(:,1:nn)=hunk(:,ns:ne,k,s)
       else
          do ib1=ns,ne,MB_d
             ib2=min(ib1+MB_d-1,ne)
             call hamiltonian &
                  (k,s,unk(:,ib1:ib2,k,s),hpsi(:,ib1-ns+1:ib2-ns+1),n1,n2,ib1,ib2)
          end do
       end if

       if ( j0<LLD_C ) i0=0

       if ( j0==0 ) then
          do ms=1,ns-1,MBSIZE
             me=min(ms+MBSIZE-1,ns-1)
             mm=me-ms+1
             IPROW  = mod( (ms-1)/MBSIZE, NPROW )
             iroot1 = usermap(IPROW,IPCOL,1)
             if ( iroot1 == myrank ) i0=i0+mm
          end do
       end if

       do ms=ns,ne,MBLK
          me=min(ms+MBLK-1,ne)
          mm=me-ms+1
          IPROW  = mod( (ms-1)/MBSIZE, NPROW )
          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)
          allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
          vtmp2=zero ; wtmp2=zero
          MBLKH=max(MBLK/2,NBLK1)
          call mate_sub(ms,me,ns,ne,MBLKH,ns,ms,me,k,s)
          call mpi_reduce &
               (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)
          if ( iroot1 == myrank ) then
             Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
             i0=i0+mm
          end if
          deallocate( wtmp2,vtmp2 )
       end do

       if ( ne+1<=MB_1 ) then
          do ms=ne+1,MB_1,MBLK
             me=min(ms+MBLK-1,MB_1)
             mm=me-ms+1
             IPROW  = mod( (ms-1)/MBSIZE, NPROW )
             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)
             allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
             vtmp2=zero ; wtmp2=zero
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,1),ML0,zero,vtmp2(ms,ns),mm)
#else
             call zgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,1),ML0,zero,vtmp2(ms,ns),mm)
#endif
             call mpi_reduce &
                  (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                i0=i0+mm
             end if
             deallocate( wtmp2,vtmp2 )
          end do
       end if

       nns = max(ns,MB_1-mat_block(mrnk,2)+1)
       nnn = ne-nns+1
       mme = min(MB,MB_1+mat_block(mrnk,1))

       do ms=MB_1+1,mme,MBLK
          me=min(ms+MBLK-1,mme)
          mm=me-ms+1
          IPROW  = mod( (ms-1)/MBSIZE, NPROW )
          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)
          j1=j0
          if ( ns<nns ) then
             if ( iroot1 == myrank ) then
                i=mod( (ns-1)/MBSIZE, NPROW )
                j=mod( (ms-1)/NBSIZE, NPCOL )
                n=usermap(i,j,1)
                if ( n == myrank ) stop 'subspace_mate(1)'
                nrecv_me=nrecv_me+1
                irecv_me(nrecv_me,0)=n
                irecv_me(nrecv_me,1)=ms
                irecv_me(nrecv_me,2)=me
                irecv_me(nrecv_me,3)=ns
                irecv_me(nrecv_me,4)=min(ne,nns-1)
                irecv_me(nrecv_me,5)=i0+1
                irecv_me(nrecv_me,6)=i0+mm
                irecv_me(nrecv_me,7)=j1+1
                irecv_me(nrecv_me,8)=j1+min(ne,nns-1)-ns+1
                j1=j1+min(ne,nns-1)-ns+1
             end if
          end if
          if ( nns<=ne ) then
             allocate( vtmp2(ms:me,nns:ne),wtmp2(ms:me,nns:ne) )
             vtmp2=zero ; wtmp2=zero
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,mm,nnn,ML0, dV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,nns-ns+1),ML0,zero,vtmp2,mm)
#else
             call zgemm(TRANSA,TRANSB,mm,nnn,ML0,zdV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,nns-ns+1),ML0,zero,vtmp2,mm)
#endif
             call mpi_reduce &
                  (vtmp2,wtmp2,mm*nnn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                Hsub(i0+1:i0+mm,j1+1:j1+nnn)=wtmp2(ms:me,nns:ne)
             end if
             deallocate( wtmp2,vtmp2 )
          end if
          if ( iroot1 == myrank ) then
             i0=i0+mm
          end if
       end do

       if ( mme+1<=MB ) then
          do ms=mme+1,MB,MBLK
             me=min(ms+MBLK-1,MB)
             mm=me-ms+1
             IPROW  = mod( (ms-1)/MBSIZE, NPROW )
             iroot1 = usermap(IPROW,IPCOL,1)
             if ( iroot1 == myrank ) then
                i=mod( (ns-1)/MBSIZE, NPROW )
                j=mod( (ms-1)/NBSIZE, NPCOL )
                n=usermap(i,j,1)
                if ( n == myrank ) stop 'subspace_mate(2)'
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
             end if
          end do
       end if

       if ( MB+1<=MB_1+mat_block(mrnk,1) ) then
          i0=0
          do mms=MB+1,MB_1+mat_block(mrnk,1),MBLK
             mme=min(mms+MBLK-1,MB_1+mat_block(mrnk,1))
             ms=mod(mms+MB-1,MB)+1
             me=mod(mme+MB-1,MB)+1
             mm=me-ms+1
             IPROW  = mod( (ms-1)/MBSIZE,NPROW )
             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)
             allocate( vtmp2(ms:me,nns:ne),wtmp2(ms:me,nns:ne) )
             vtmp2=zero ; wtmp2=zero
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,mm,nnn,ML0, dV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,nns-ns+1),ML0,zero,vtmp2,mm)
#else
             call zgemm(TRANSA,TRANSB,mm,nnn,ML0,zdV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,nns-ns+1),ML0,zero,vtmp2,mm)
#endif
             call mpi_reduce(vtmp2,wtmp2,mm*nnn,TYPE_MAIN,mpi_sum &
                  ,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                i=mod( (nns-1)/MBSIZE, NPROW )
                j=mod( (ms-1)/NBSIZE, NPCOL )
                n=usermap(i,j,1)
                if ( n == myrank ) stop 'subspace_mate(3)'
                nsend_me=nsend_me+1
                isend_me(nsend_me,0)=n
                isend_me(nsend_me,1)=ms
                isend_me(nsend_me,2)=me
                isend_me(nsend_me,3)=nns
                isend_me(nsend_me,4)=ne
                isend_me(nsend_me,5)=i0+1
                isend_me(nsend_me,6)=i0+mm
                isend_me(nsend_me,7)=j0+1
                isend_me(nsend_me,8)=j0+nnn
                Hsub(i0+1:i0+mm,j0+1:j0+nnn)=wtmp2(ms:me,nns:ne)
                i0=i0+mm
             end if
             deallocate( wtmp2,vtmp2 )
          end do
       end if

       if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

    end do ! ns

    deallocate( hpsi )

! --- ME-send ---

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
          call mpi_isend(vtmp2(1,i),j,TYPE_MAIN,n,itag &
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
          call mpi_irecv(wtmp2(1,i),j,TYPE_MAIN,n,itag &
               ,mpi_comm_world,ireq(nreq),ierr)
       end do
       call mpi_waitall(nreq,ireq,istatus,ierr)

       do i=1,nrecv_me
          i1=irecv_me(i,5)
          i2=irecv_me(i,6)
          j1=irecv_me(i,7)
          j2=irecv_me(i,8)
          if ( TYPE_MAIN==mpi_complex16 ) then
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   ztmp=Hsub(i0,j0)
                   ztmp=conjg(ztmp)
                   Hsub(i0,j0)=ztmp
                end do
             end do
          else
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   Hsub(i0,j0)=wtmp2(j,i)
                end do
             end do
          end if
       end do

    end if

    deallocate( wtmp2 )
    deallocate( vtmp2 )
    deallocate( irecv_me, isend_me )

    call result_timer( "mate_sl", t )
    call write_border( 1, " subspace_mate_sl(end)" )

  END SUBROUTINE subspace_mate_sl


  RECURSIVE SUBROUTINE mate_sub(mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,k,s)
    integer,intent(IN)    :: mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,k,s
    integer :: n1,n2,ML0,n,ns,ne,nn,m,ms,me,mm,mms,MBLKH,i,ld

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    ld  = ld1-ld0+1

    do ns=nn1,nn2,MBLK
       ne=min(ns+MBLK-1,nn2)
       nn=ne-ns+1
       if ( nn<=0 ) cycle
       do mms=mm1,mm2,MBLK
          ms=max(ns,mms)
          me=min(mms+MBLK-1,mm2)
          mm=me-ms+1
          if ( mm<=0 ) cycle
          if ( ms>=ne ) then
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)
#else
             call zgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)
#endif
          else if ( mm<=NBLK1 ) then
             do n=ns,ne
#ifdef _DRSDFT_
                call dgemv(TRANSA,ML0,ne-n+1, dV,unk(n1,n,k,s) &
                     ,ML0,hpsi(n1,n-ns0+1),1,zero,vtmp2(n,n),1)
#else
                call zgemv(TRANSA,ML0,ne-n+1,zdV,unk(n1,n,k,s) &
                     ,ML0,hpsi(n1,n-ns0+1),1,zero,vtmp2(n,n),1)
#endif
             end do
          else
             MBLKH=max(MBLK/2,NBLK1)
             call mate_sub(ms,me,ns,ne,MBLKH,ns0,ld0,ld1,k,s)
          end if

       end do ! mms
    end do ! ns

    return
  END SUBROUTINE mate_sub

!-------------------------------------------------

  SUBROUTINE subspace_mate_sl_bp1( k, s, u, v )
    implicit none
    integer,intent(IN) :: k,s
    real(8),intent(IN) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i0,i1,i2,ib1,ib2,j,j0,j1,j2,m,me,mm,ms,MB
    integer :: MB_0,MB_1,myrank,n1,n2
    integer :: MBLK,MBLKH,ML0,mms,mme,nnn,nns,n,ne,nn,ns
    include 'mpif.h'
!   integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer,allocatable :: istatus(:,:),ireq(:)
    integer :: nreq,itag,ierr,nwork
    integer :: IPROW,IPCOL,iroot1,iroot2,mrnk,nrecv_me,nsend_me,comm_grid
    integer,allocatable :: ir(:),id(:),irecv_me(:,:),isend_me(:,:)
    complex(8) :: ztmp
    integer :: ncycle, icycle, nbss, ib, ncycle2
    integer :: comm_band, np_band, myrank_b, irank_b
    integer :: ix0, msV, meV, nsV, neV
    real(8),allocatable :: hu2(:,:), uB(:,:)

    call write_border( 1, " subspace_mate_sl_bp1(start)" )

    UPLO  = 'L'
    MB    = MB_diag
    MB_0  = v(2)%pinfo%id(v(2)%pinfo%me)+1
    MB_1  = MB_0+v(2)%pinfo%ir(v(2)%pinfo%me)-1
    NBLK1 = 4
    MBLK  = min( MBSIZE,NBSIZE )
    dV    = v(1)%factor
    ML0   = size( u, 1 )
    n1    = v(1)%pinfo%id(v(1)%pinfo%me)+1
    n2    = n1+v(1)%pinfo%ir(v(1)%pinfo%me)-1
    mrnk  = v(2)%pinfo%me
    comm_grid = v(1)%pinfo%comm

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

!   allocate( irecv_me(99,0:8),isend_me(99,0:8) )
    nwork=(MB-1)/MBLK+1
    nwork=((nwork-1)/nprow+1)*((nwork-1)/npcol+1)
    allocate( istatus(MPI_STATUS_SIZE,nwork), ireq(nwork) )
    allocate( irecv_me(nwork,0:8),isend_me(nwork,0:8) )

    irecv_me(:,:) =-1
    isend_me(:,:) =-1

!--------------------------------------------------------

    np_band   = v(2)%pinfo%np
    myrank_b  = v(2)%pinfo%me
    comm_band = v(2)%pinfo%comm
    ncycle    = (MB-1)/MBLK+1
    ncycle2   = (ncycle-1)/np_band+1

    allocate( uB(ML0,MBLK) ); uB=zero
    allocate( hu2(ML0,ncycle2*MBLK) ); hu2=zero

    i0=0
    ix0=0
    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band ) + 1
       meV = min( msV+MBLK-1, MB )

       if ( irank_b == myrank_b ) then
          do ib1=msV,meV
             ib2=ib1
             call hamiltonian(k,s,u(:,ib1:ib2),hu2(:,ix0+ib1-msV+1:ix0+ib2-msV+1),n1,n2,ib1,ib2)
          end do
          ix0=ix0+meV-msV+1
       end if

    enddo

    nrecv_me=0
    nsend_me=0

    ix0=0
    do icycle=1,ncycle
       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band )+ 1
       meV = min( msV+MBLK-1, MB )

       j0=0

!      if ( irank_b == myrank_b ) uB(:,:)=u(:,msV:meV)
       if ( irank_b == myrank_b ) uB(:,1:mm)=u(:,msV:meV)

       call rsdft_bcast( uB, irank_b, comm_band)

       ix0=0
       do ib=1,ncycle2

          nbss = (ib-1)*np_band + myrank_b + 1

          if ( nbss <= ncycle .and. nbss < icycle ) then

             if (nLB(icycle)+nbss>0) then

                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1

                nsV = MBLK*(ib-1)+1
                neV = min(nsV+MBLK-1,MB)

                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                iroot2 = usermap(IPROW,IPCOL,2)
                allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
                vtmp2=zero ; wtmp2=zero
                call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                          ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
                ix0=ix0+nn
                call mpi_reduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN &
                               ,mpi_sum,iroot2,comm_grid,ierr)
                if ( iroot1 == myrank ) then
                   Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                   j0=j0+nn
                end if
                deallocate( wtmp2,vtmp2 )
             else   ! nLB(icycle)<0
                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1
                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                if ( iroot1 == myrank ) then     ! <------- (N)
                   i=mod( (ns-1)/MBLK, NPROW )
                   j=mod( (ms-1)/MBLK, NPCOL )
                   n=usermap(i,j,1)
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
                   j0=j0+nn    
                end if
                ix0=ix0+nn
             endif

          else if ( nbss <= ncycle .and. nbss > icycle ) then

             if ( nbss-nLB(icycle)>=0 .and. nLB(icycle)>0 ) then

                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1

                nsV = MBLK*(ib-1)+1
                neV = min(nsV+MBLK-1,MB)

                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                iroot2 = usermap(IPROW,IPCOL,2)
                allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
                vtmp2=zero ; wtmp2=zero
                ix0=MBLK*(ib-1)
                call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                          ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
                ix0=ix0+nn
                call mpi_reduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN &
                               ,mpi_sum,iroot2,comm_grid,ierr)

                if ( iroot1 == myrank ) then
                   i=mod( (ns-1)/MBSIZE, NPROW )
                   j=mod( (ms-1)/NBSIZE, NPCOL )
                   n=usermap(i,j,1)
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
                   Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                   j0=j0+nn
                end if
                deallocate( wtmp2,vtmp2 )

             endif
         
          else if ( nbss <= ncycle .and. nbss == icycle ) then

             ns = MBLK*(nbss-1)+1
             ne = min(ns+MBLK-1,MB)
             nn=ne-ns+1

             nsV = MBLK*(ib-1)+1
             neV = min(nsV+MBLK-1,MB)

             IPCOL  = mod( (ns-1)/MBLK, NPCOL )
             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)
             allocate( vtmp2(ms:me,ns:ne) ) ; vtmp2=zero
             allocate( wtmp2(ms:me,ns:ne) ) ; wtmp2=zero
             MBLKH=max(MBLK/2,NBLK1)
             call mate_sub_bp(ms,me,ns,ne,MBLKH,ns,ms,me &
                             ,uB,hu2(:,ix0+1:ix0+nn))
             ix0=ix0+nn              
             call mpi_reduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN &
                            ,mpi_sum,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                j0=j0+nn
             end if
             deallocate( wtmp2,vtmp2 )

          end if

       end do ! ib

       if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

    end do ! icycle  (column loop)

!   exchange symetrical data (B_n,m <->B_m,n)
! --- ME-send ---

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
!-debug
             Hsub(i0,j0)=0.d0
!-debug
          end do
          end do
          itag=ns+10*ne+100*ms+1000*me
          nreq=nreq+1
          call mpi_isend(vtmp2(1,i),j,TYPE_MAIN,n,itag &
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
          call mpi_irecv(wtmp2(1,i),j,TYPE_MAIN,n,itag &
                        ,mpi_comm_world,ireq(nreq),ierr)
       end do
       call mpi_waitall(nreq,ireq,istatus,ierr)

       do i=1,nrecv_me
          i1=irecv_me(i,5)
          i2=irecv_me(i,6)
          j1=irecv_me(i,7)
          j2=irecv_me(i,8)
          if ( TYPE_MAIN==mpi_complex16 ) then
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   ztmp=Hsub(i0,j0)
                   ztmp=conjg(ztmp)
                   Hsub(i0,j0)=ztmp
                end do
             end do
          else
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   Hsub(i0,j0)=wtmp2(j,i)
                end do
             end do
          end if
       end do

    end if

    deallocate( wtmp2, vtmp2 )

    deallocate (uB)
    deallocate (hu2)
    deallocate( irecv_me, isend_me )

    call write_border( 1, " subspace_mate_sl_bp1(end)" )

  END SUBROUTINE subspace_mate_sl_bp1

  SUBROUTINE subspace_mate_sl_bp2( k, s, u, v )
    implicit none
    integer,intent(IN) :: k,s
    real(8),intent(IN) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i0,i1,i2,ib1,ib2,j,j0,j1,j2,m,me,mm,ms,MB
    integer :: MB_0,MB_1,myrank,n1,n2
    integer :: MBLK,MBLKH,ML0,mms,mme,nnn,nns,n,ne,nn,ns
    include 'mpif.h'
!   integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer,allocatable :: istatus(:,:),ireq(:)
    integer :: nreq,itag,ierr,nwork
    integer :: IPROW,IPCOL,iroot1,iroot2,mrnk,nrecv_me,nsend_me,comm_grid
    integer,allocatable :: ir(:),id(:),irecv_me(:,:),isend_me(:,:)
    complex(8) :: ztmp
    integer :: ncycle, icycle, nbss, ib, ncycle2
    integer :: comm_band, np_band, myrank_b, irank_b
    integer :: ix0, msV, meV, nsV, neV
    real(8),allocatable :: hu2(:,:), uB(:,:)

    call write_border( 1, " subspace_mate_sl_bp2(start)" )

    UPLO  = 'L'
    MB    = MB_diag
    MB_0  = v(2)%pinfo%id(v(2)%pinfo%me)+1
    MB_1  = MB_0+v(2)%pinfo%ir(v(2)%pinfo%me)-1
    NBLK1 = 4
    MBLK  = min( MBSIZE,NBSIZE )
    dV    = v(1)%factor
    ML0   = size( u, 1 )
    n1    = v(1)%pinfo%id(v(1)%pinfo%me)+1
    n2    = n1+v(1)%pinfo%ir(v(1)%pinfo%me)-1
    mrnk  = v(2)%pinfo%me
    comm_grid = v(1)%pinfo%comm

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

!   allocate( irecv_me(99,0:8),isend_me(99,0:8) )
    nwork=(MB-1)/MBLK+1
    nwork=((nwork-1)/nprow+1)*((nwork-1)/npcol+1)
    allocate( istatus(MPI_STATUS_SIZE,nwork), ireq(nwork) )
    allocate( irecv_me(nwork,0:8),isend_me(nwork,0:8) )

    irecv_me(:,:) =-1
    isend_me(:,:) =-1

!--------------------------------------------------------

    np_band   = v(2)%pinfo%np
    myrank_b  = v(2)%pinfo%me
    comm_band = v(2)%pinfo%comm
    ncycle    = (MB-1)/MBLK+1
    ncycle2   = (ncycle-1)/np_band+1

    allocate( uB(ML0,MBLK) ); uB=zero
    allocate( hu2(ML0,ncycle2*MBLK) ); hu2=zero
 
    i0=0
    ix0=0
    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band )+ 1
       meV = min( msV+MBLK-1, MB )
 
       if ( irank_b == myrank_b) then
          do ib1=msV,meV
             ib2=ib1
             call hamiltonian(k,s,u(:,ib1:ib2),hu2(:,ix0+ib1-msV+1:ix0+ib2-msV+1),n1,n2,ib1,ib2)
          end do
          ix0=ix0+meV-msV+1
       end if

    enddo

    nrecv_me=0
    nsend_me=0

    ix0=0
    do icycle=1,ncycle
       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band )+ 1
       meV = min( msV+MBLK-1, MB )

       j0=0

!      if ( irank_b == myrank_b ) uB(:,:)=u(:,msV:meV)
       if ( irank_b == myrank_b ) uB(:,1:mm)=u(:,msV:meV)

       call rsdft_bcast( uB, irank_b, comm_band)

       ix0=0
       do ib=1,ncycle2

          nbss = (ib-1)*np_band + myrank_b + 1

          if ( nbss <= ncycle ) then

             ns = MBLK*(nbss-1)+1
             ne = min(ns+MBLK-1,MB)
             nn=ne-ns+1

             nsV = MBLK*(ib-1)+1
             neV = min(nsV+MBLK-1,MB)

             IPCOL  = mod( (ns-1)/MBLK, NPCOL )
             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)
             allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
             vtmp2=zero ; wtmp2=zero
             call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                       ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
             ix0=ix0+nn
             call mpi_reduce &
                  (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                j0=j0+nn
             end if
             deallocate( wtmp2,vtmp2 )

          end if

       end do

       if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

    end do ! icycle  (column loop)

!   exchange symetrical data (B_n,m <->B_m,n)
! --- ME-send ---

    deallocate (uB)
    deallocate (hu2)
    deallocate( irecv_me, isend_me )

    call write_border( 1, " subspace_mate_sl_bp2(end)" )

  END SUBROUTINE subspace_mate_sl_bp2

  SUBROUTINE subspace_mate_sl_bp3( k, s, u, v )
    implicit none
    integer,intent(IN) :: k,s
    real(8),intent(IN) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i0,i1,i2,ib1,ib2,j,j0,j1,j2,m,me,mm,ms,MB
    integer :: MB_0,MB_1,myrank,n1,n2
    integer :: MBLK,MBLKH,ML0,mms,mme,nnn,nns,n,ne,nn,ns
    include 'mpif.h'
!   integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer,allocatable :: istatus(:,:),ireq(:)
    integer :: nreq,itag,ierr,nwork
    integer :: IPROW,IPCOL,iroot1,iroot2,mrnk,nrecv_me,nsend_me,comm_grid
    integer,allocatable :: ir(:),id(:),irecv_me(:,:),isend_me(:,:)
    complex(8) :: ztmp
    integer :: ncycle, icycle, nbss, ib, ncycle2
    integer :: comm_band, np_band, myrank_b, irank_b
    integer :: ix0, msV, meV, nsV, neV
    real(8),allocatable :: hu2(:,:), uB(:,:)

    call write_border( 1, " subspace_mate_sl_bp3(start)" )

    UPLO  = 'L'
    MB    = MB_diag
    MB_0  = v(2)%pinfo%id(v(2)%pinfo%me)+1
    MB_1  = MB_0+v(2)%pinfo%ir(v(2)%pinfo%me)-1
    NBLK1 = 4
    MBLK  = min( MBSIZE,NBSIZE )
    dV    = v(1)%factor
    ML0   = size( u, 1 )
    n1    = v(1)%pinfo%id(v(1)%pinfo%me)+1
    n2    = n1+v(1)%pinfo%ir(v(1)%pinfo%me)-1
    mrnk  = v(2)%pinfo%me
    comm_grid = v(1)%pinfo%comm

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

!   allocate( irecv_me(99,0:8),isend_me(99,0:8) )
    nwork=(MB-1)/MBLK+1
    nwork=((nwork-1)/nprow+1)*((nwork-1)/npcol+1)
    allocate( istatus(MPI_STATUS_SIZE,nwork), ireq(nwork) )
    allocate( irecv_me(nwork,0:8),isend_me(nwork,0:8) )

    irecv_me(:,:) =-1
    isend_me(:,:) =-1

!--------------------------------------------------------

    np_band   = v(2)%pinfo%np
    myrank_b  = v(2)%pinfo%me
    comm_band = v(2)%pinfo%comm
    ncycle    = (MB-1)/MBLK+1
    ncycle2   = (ncycle-1)/np_band+1

    allocate( uB(ML0,MBLK) ); uB=zero
    allocate( hu2(ML0,ncycle2*MBLK) ); hu2=zero
 

    i0=0
    ix0=0
    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band )+ 1
       meV = min( msV+MBLK-1, MB )
 
       if ( irank_b == myrank_b) then
          do ib1=msV,meV
             ib2=ib1
             call hamiltonian(k,s,u(:,ib1:ib2),hu2(:,ix0+ib1-msV+1:ix0+ib2-msV+1),n1,n2,ib1,ib2)
          end do
          ix0=ix0+meV-msV+1
       end if

    enddo

    nrecv_me=0
    nsend_me=0

    ix0=0
    do icycle=1,ncycle
       irank_b = mod( icycle-1, np_band )

       ms = MBLK*(icycle-1)+1
       me = min( ms+MBLK-1, MB)
       mm = me-ms+1
       IPROW = mod( (ms-1)/MBLK, NPROW )

       msV = MBLK*int((icycle-1)/np_band )+ 1
       meV = min( msV+MBLK-1, MB )

       j0=0

!      if ( irank_b == myrank_b ) uB(:,:)=u(:,msV:meV)
       if ( irank_b == myrank_b ) uB(:,1:mm)=u(:,msV:meV)

       call rsdft_bcast( uB, irank_b, comm_band)

       ix0=0
       do ib=1,ncycle2

          nbss = (ib-1)*np_band + myrank_b + 1

          if ( nbss <= ncycle .and. nbss < icycle ) then
             if (nLB(icycle)+nbss>0) then

                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1

                nsV = MBLK*(ib-1)+1
                neV = min(nsV+MBLK-1,MB)

                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                iroot2 = usermap(IPROW,IPCOL,2)
                allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
                vtmp2=zero ; wtmp2=zero
                call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                     ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
                ix0=ix0+nn
                call mpi_reduce &
                     (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

                if ( iroot1 == myrank ) then
                   i=mod( (ns-1)/MBSIZE, NPROW )
                   j=mod( (ms-1)/NBSIZE, NPCOL )
                   n=usermap(i,j,1)
                   if ( n == myrank ) then
                      write(*,*) "subspace_mat(3):",n,i,j,ns,ms,icycle,ib
                   endif
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
                   Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                   j0=j0+nn
                end if
                deallocate( wtmp2,vtmp2 )
             else   ! nLB(icycle)<0
                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1
                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                if ( iroot1 == myrank ) then     ! <------- (N)
                   i=mod( (ns-1)/MBLK, NPROW )
                   j=mod( (ms-1)/MBLK, NPCOL )
                   n=usermap(i,j,1)
                   if ( n == myrank ) then
                      write(*,*) "subspace_mat(2):",n,i,j,ns,ms,icycle,ib
                   endif

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
                   j0=j0+nn    
                end if
                ix0=ix0+nn
             endif

          else if ( nbss <= ncycle .and. nbss > icycle ) then

             if (nbss-nLB(icycle)>=0 .and. nLB(icycle)>0 ) then

                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1

                nsV = MBLK*(ib-1)+1
                neV = min(nsV+MBLK-1,MB)

                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                iroot2 = usermap(IPROW,IPCOL,2)
                allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
                vtmp2=zero ; wtmp2=zero
                ix0=MBLK*(ib-1)
                call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                     ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
                ix0=ix0+nn
                call mpi_reduce &
                     (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

                if ( iroot1 == myrank ) then
                   i=mod( (ns-1)/MBSIZE, NPROW )
                   j=mod( (ms-1)/NBSIZE, NPCOL )
                   n=usermap(i,j,1)
                   if ( n == myrank ) then
                      write(*,*) "subspace_mat(3):",n,i,j,ns,ms,icycle,ib
                   endif
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
                   Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                   j0=j0+nn
                end if
                deallocate( wtmp2,vtmp2 )

             else

                ns = MBLK*(nbss-1)+1
                ne = min(ns+MBLK-1,MB)
                nn=ne-ns+1
                IPCOL  = mod( (ns-1)/MBLK, NPCOL )
                iroot1 = usermap(IPROW,IPCOL,1)
                if ( iroot1 == myrank ) then     ! <------- (N)
                   i=mod( (ns-1)/MBLK, NPROW )
                   j=mod( (ms-1)/MBLK, NPCOL )
                   n=usermap(i,j,1)
                   if ( n == myrank ) then
                      write(*,*) "subspace_mat(2):",n,i,j,ns,ms,icycle,ib
                   endif

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
                   j0=j0+nn
                end if
                ix0=ix0+nn

             end if

          else if ( nbss <= ncycle .and. nbss == icycle ) then

             ns = MBLK*(nbss-1)+1
             ne = min(ns+MBLK-1,MB)
             nn=ne-ns+1

             nsV = MBLK*(ib-1)+1
             neV = min(nsV+MBLK-1,MB)

             IPCOL  = mod( (ns-1)/MBLK, NPCOL )
             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)

             allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) )
             vtmp2=zero ; wtmp2=zero
             ix0=MBLK*(ib-1)
             call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,uB &
                       ,ML0,hu2(:,ix0+1),ML0,zero,vtmp2(ms,ns),mm)
             ix0=ix0+nn
             call mpi_reduce &
                  (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
                Hsub(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                j0=j0+nn
             endif
             deallocate( wtmp2,vtmp2 )

          end if

       end do

       if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

    end do ! icycle  (column loop)

!debug
!    write(*,*) "myrank=",myrank," nsend=",nsend_me," nrecv=",nrecv_me
!debug
!        exchange symetrical data (B_n,m <->B_m,n)
! --- ME-send ---

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
!-debug
!            Hsub(i0,j0)=0.d0
!-debug
          end do
          end do
          itag=ns+10*ne+100*ms+1000*me
          nreq=nreq+1
          call mpi_isend(vtmp2(1,i),j,TYPE_MAIN,n,itag &
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
          call mpi_irecv(wtmp2(1,i),j,TYPE_MAIN,n,itag &
               ,mpi_comm_world,ireq(nreq),ierr)
       end do
       call mpi_waitall(nreq,ireq,istatus,ierr)

       do i=1,nrecv_me
          i1=irecv_me(i,5)
          i2=irecv_me(i,6)
          j1=irecv_me(i,7)
          j2=irecv_me(i,8)
          if ( TYPE_MAIN==mpi_complex16 ) then
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   ztmp=Hsub(i0,j0)
                   ztmp=conjg(ztmp)
                   Hsub(i0,j0)=ztmp
                end do
             end do
          else
             j=0
             do j0=j1,j2
                do i0=i1,i2
                   j=j+1
                   Hsub(i0,j0)=wtmp2(j,i)
                end do
             end do
          end if
       end do

    end if

    deallocate( wtmp2, vtmp2 )
    deallocate (uB)
    deallocate (hu2)
    deallocate( irecv_me, isend_me )

    call write_border( 1, " subspace_mate_sl_bp3(end)" )

  END SUBROUTINE subspace_mate_sl_bp3

  RECURSIVE SUBROUTINE mate_sub_bp(mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,u,hu)
    implicit none
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1
    real(8),intent(IN) :: u(:,:),hu(:,:)
    integer :: n1,n2,ML0,n,ns,ne,nn,m,ms,me,mm,mms,MBLKH,i,ld
    integer :: ms0

    ML0 = size(u,1)
    ld  = ld1-ld0+1
    ms0 = ld0

    do ns=nn1,nn2,MBLK

       ne=min(ns+MBLK-1,nn2)
       nn=ne-ns+1

       if ( nn <= 0 ) cycle

       do mms=mm1,mm2,MBLK

          ms=max(ns,mms)
          me=min(mms+MBLK-1,mm2)
          mm=me-ms+1

          if ( mm <= 0 ) cycle

          if ( ms >= ne ) then
             call dgemm(TRANSA,TRANSB,mm,nn,ML0,dV,u(1,ms-ms0+1) &
                  ,ML0,hu(1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)
          else if ( mm <= NBLK1 ) then
             do n=ns,ne
                call dgemv(TRANSA,ML0,ne-n+1,dV,u(1,n-ms0+1) &
                     ,ML0,hu(1,n-ns0+1),1,zero,vtmp2(n,n),1)
             end do
          else
             MBLKH=max(MBLK/2,NBLK1)
             call mate_sub_bp(ms,me,ns,ne,MBLKH,ns0,ld0,ld1,u,hu)
          end if

       end do ! mms

    end do ! ns

    return
  END SUBROUTINE mate_sub_bp


END MODULE subspace_mate_sl_module
