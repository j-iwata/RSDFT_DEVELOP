MODULE subspace_rotv_sl_module

  use rgrid_module, only: dV,zdV
  use wf_module, only: unk
!  use electron_module
  use scalapack_module
  use parallel_module
  use subspace_diag_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1

  implicit none

  PRIVATE
  PUBLIC :: subspace_rotv_sl

#ifdef _DRSDFT_
  real(8),allocatable :: utmp(:,:),utmp2(:,:)
#else
  complex(8),allocatable :: utmp(:,:),utmp2(:,:)
#endif

CONTAINS


  SUBROUTINE subspace_rotv_sl(k,s)
    integer,intent(IN) :: k,s
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    MB  = MB_diag

    NBLK2 = maxval(ircnt)

    allocate( utmp(NBLK2,MB_0:MB_1) )
    utmp=zero

    do i=1,maxval(ircnt),NBLK2
       i1=n1+i-1
       i2=min(i1+NBLK2-1,n2)
       ii=i2-i1+1

       utmp(:,:)=zero

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
                utmp2(ms:me,ns:ne)=Vsub(i0+1:i0+mm,j0+1:j0+nn)
                i0=i0+mm
             end if

             call mpi_bcast(utmp2(ms,ns),mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)

             if ( ii>0 ) then
#ifdef _DRSDFT_
                call dgemm('N','N',ii,nn,mm,one,unk(i1,ms,k,s) &
                     ,ML0,utmp2(ms,ns),mm,one,utmp(1,ns),NBLK2)
#else
                call zgemm('N','N',ii,nn,mm,one,unk(i1,ms,k,s) &
                     ,ML0,utmp2(ms,ns),mm,one,utmp(1,ns),NBLK2)
#endif
             end if

             deallocate( utmp2 )

          end do ! ms

          if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

       end do ! ns

       if ( ii>0 ) then
          unk(i1:i2,MB_0:MB_1,k,s)=utmp(1:ii,MB_0:MB_1)
       end if

    end do ! ii

    deallocate( utmp )

  END SUBROUTINE subspace_rotv_sl

END MODULE subspace_rotv_sl_module
