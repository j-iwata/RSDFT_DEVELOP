module subspace_gthr_sl_module

  use wf_module, only: unk,esp !,hunk,iflag_hunk,USE_WORKWF_AT_ROTV
  use scalapack_module
  use parallel_module
  use subspace_diag_variables, only: MB_diag,NBLK2,zero,one,Hsub,Vsub,TYPE_MAIN
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use bcast_module
  use watch_module
  use rsdft_mpi_module

  implicit none

  private
  public :: subspace_gthr_sl

#ifdef _DRSDFT_
  real(8),allocatable :: utmp(:,:),utmp2(:,:)
#else
  complex(8),allocatable :: utmp(:,:),utmp2(:,:)
#endif

contains


  subroutine subspace_gthr_sl(k,s)
    implicit none
    integer,intent(in) :: k,s
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm,loop
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB,n_loop
    type(time) :: t
    real(8),allocatable :: H(:,:)
    integer :: LWORK
    real(8),allocatable :: rwork(:)
    real(8) :: rtmp(1)

    call write_border( 1, " subspace_gthr_sl(start)" )
    call start_timer( t )

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    MB  = MB_diag

    allocate( H(MB,MB) ); H=0.0d0

    j0=0
    do ns=MB_0,MB_1,NBSIZE
       ne=min(ns+NBSIZE-1,MB_1)
       nn=ne-ns+1

       IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

       i0=0
       do ms=1,MB,MBSIZE
          me=min(ms+MBSIZE-1,MB)
          mm=me-ms+1

          IPROW=mod( (ms-1)/MBSIZE,NPROW )

          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)

          if ( mm<1 .or. nn<1 ) cycle

          allocate( utmp2(ms:me,ns:ne) )

          if ( iroot1 == myrank ) then
!$OMP parallel workshare
             utmp2(ms:me,ns:ne)=Hsub(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
             i0=i0+mm
          end if

#ifdef _DRSDFT_
          call d_rsdft_bcast(utmp2,mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)
#else
          call z_rsdft_bcast(utmp2,mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)
#endif

          H(ms:me,ns:ne) = utmp2(ms:me,ns:ne)

          deallocate( utmp2 )

       end do ! ms

       if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

    end do ! ns

    if ( myrank == 0 ) then
       ii=0
       do i2=1,MB
          do i1=1,i2
             if ( H(i1,i2) /= 0.0d0 ) ii=ii+1 
          end do
       end do
       mm=0
       do i2=1,MB
          do i1=i2,MB
             if ( H(i1,i2) /= 0.0d0 ) mm=mm+1
          end do
       end do
       write(*,*) "count(H/=0.0d0),U,L=",count(H/=0.0d0),ii,mm
    end if

! ---

    !if ( LWORK==0 ) LWORK=3*MB-1
    call DSYEV('V','L',MB,H,MB,esp(1,k,s),rtmp,-1,ierr)
    LWORK=nint(rtmp(1))
    allocate( rwork(LWORK) ); rwork=0.0d0

    call DSYEV('V','L',MB,H,MB,esp(1,k,s),rwork,LWORK,ierr)

    deallocate( rwork )

! ---

    j0=0
    do ns=MB_0,MB_1,NBSIZE
       ne=min(ns+NBSIZE-1,MB_1)
       nn=ne-ns+1

       IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

       i0=0
       do ms=1,MB,MBSIZE
          me=min(ms+MBSIZE-1,MB)
          mm=me-ms+1

          IPROW=mod( (ms-1)/MBSIZE,NPROW )

          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)

          if ( mm<1 .or. nn<1 ) cycle

          if ( iroot1 == myrank ) then
!$OMP parallel workshare
             Vsub(i0+1:i0+mm,j0+1:j0+nn)=H(ms:me,ns:ne)
!$OMP end parallel workshare
             i0=i0+mm
          end if

       end do ! ms

       if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

    end do ! ns

! ---

    deallocate( H )

    call result_timer( "gthr_sl", t )
    call write_border( 1, " subspace_gthr_sl(end)" )

  end subroutine subspace_gthr_sl



  subroutine subspace_rotv2_sl( u )
    implicit none
    real(8),intent(inout) :: u(:,:)
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm,loop
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB,MB0
    type(time) :: t
    real(8),parameter :: zero=0.0d0, one=1.0d0

    call write_border( 1, " subspace_rotv2(start)" )
    call start_timer( t )

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    MB  = MB_diag
    MB0 = MB_1-MB_0+1

    NBLK2 = ML0 !maxval(ircnt)

    allocate( utmp(NBLK2,MB_0:MB_1) ); utmp=zero

    allocate( utmp2(MB,MB_0:MB_1) ); utmp2=zero

    j0=0
    do ns=MB_0,MB_1,NBSIZE
       ne=min(ns+NBSIZE-1,MB_1)
       nn=ne-ns+1

       IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

       i0=0
       do ms=1,MB,MBSIZE
          me=min(ms+MBSIZE-1,MB)
          mm=me-ms+1

          IPROW=mod( (ms-1)/MBSIZE,NPROW )

          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)

          if ( mm<1 .or. nn<1 ) cycle

          if ( iroot1 == myrank ) then
!$OMP parallel workshare
             utmp2(ms:me,ns:ne) = Vsub(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
             i0=i0+mm
          end if

       end do ! ms

       if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

    end do ! ns

    call rsdft_allreduce_sum( utmp2, comm_grid )
    call rsdft_allreduce_sum( utmp2, comm_band )

    call DGEMM('N','N',ML0,MB0,MB,one,u,ML0,utmp2(1,MB_0),MB,zero,utmp,ML0)

    u(:,MB_0:MB_1) = utmp

    deallocate( utmp2 )
    deallocate( utmp  )

    call result_timer( "rotv2_sl", t )
    call write_border( 1, " subspace_rotv2(end)" )

  end subroutine subspace_rotv2_sl


END MODULE subspace_gthr_sl_module
