MODULE subspace_rotv_sl_module

  use scalapack_module
  use subspace_diag_variables
  use vector_tools_module, only: vinfo

  implicit none

  PRIVATE
  PUBLIC :: subspace_rotv_sl
  PUBLIC :: subspace_rotv_sl_bp1
  PUBLIC :: subspace_rotv_sl_bp2

  real(8),allocatable :: utmp(:,:),utmp2(:,:)

CONTAINS


  SUBROUTINE subspace_rotv_sl( u, v )

    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm,MB_0,MB_1
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB,comm_grid
    integer :: myrank
    include 'mpif.h'

    call write_border( 1, " subspace_rotv_sl(start)" )

    ML0  = size( u, 1 )
    MB   = MB_diag
    MB_0 = v(2)%pinfo%id(v(2)%pinfo%me) + 1
    MB_1 = MB_0 + v(2)%pinfo%ir(v(2)%pinfo%me) - 1
    comm_grid = v(1)%pinfo%comm

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    NBLK2 = maxval( v(1)%pinfo%ir )

    allocate( utmp(NBLK2,MB_0:MB_1) )
    utmp=zero

    do i1=1,maxval( v(1)%pinfo%ir ),NBLK2

       i2=min(i1+NBLK2-1,ML0)
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

             if ( mm < 1 .or. nn < 1 ) cycle

             allocate( utmp2(ms:me,ns:ne) )

             if ( iroot1 == myrank ) then
!$OMP parallel workshare
                utmp2(ms:me,ns:ne)=Vsub(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
                i0=i0+mm
             end if

             call mpi_bcast(utmp2,mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)

             if ( ii > 0 ) then
                call dgemm('N','N',ii,nn,mm,one,u(i1,ms) &
                     ,ML0,utmp2(ms,ns),mm,one,utmp(1,ns),NBLK2)
             end if

             deallocate( utmp2 )

          end do ! ms

          if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

       end do ! ns

       if ( ii > 0 ) then
!$OMP parallel workshare
          u(i1:i2,MB_0:MB_1)=utmp(1:ii,MB_0:MB_1)
!$OMP end parallel workshare
       end if

    end do ! ii

    deallocate( utmp )

    call write_border( 1, " subspace_rotv_sl(end)" )

  END SUBROUTINE subspace_rotv_sl

  SUBROUTINE subspace_rotv_sl_bp1( u, v )

    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm,MB_0,MB_1
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB,comm_grid
    integer :: myrank
    include 'mpif.h'
    real(8),allocatable :: utmp3(:,:),utmp4(:,:)
    integer :: comm_band, np_band, irank_b, MBLK, myrank_b
    integer :: irank_c, ncycle, icycle, ib, nbss
    integer :: MB2, jj0

    call write_border( 1, " subspace_rotv_sl_bp1(start)" )

    ML0  = size( u, 1 )
    MB   = MB_diag
    MB_0 = v(2)%pinfo%id(v(2)%pinfo%me) + 1
    MB_1 = MB_0 + v(2)%pinfo%ir(v(2)%pinfo%me) - 1
    comm_grid = v(1)%pinfo%comm

    np_band   = v(2)%pinfo%np
    myrank_b  = v(2)%pinfo%me
    comm_band = v(2)%pinfo%comm
    MBLK=min(MBSIZE,NBSIZE)

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    NBLK2 = maxval( v(1)%pinfo%ir )

    allocate( utmp3(NBLK2,MBLK) ); utmp3=zero

    MB2=(MB-1)/np_band+1
    allocate( utmp4(NBLK2,MB2) ); utmp4=zero

    ncycle=(MB-1)/MBLK+1

    do i1=1,maxval( v(1)%pinfo%ir ),NBLK2

       i2=min(i1+NBLK2-1,ML0)
       ii=i2-i1+1

       utmp4(:,:)=zero

       i0=0
       do icycle=1,ncycle
          ms=MBLK*(icycle-1)+1
          me=min(ms+MBLK-1,MB)
          mm=me-ms+1

          irank_b = mod( icycle-1, np_band )

          IPROW=mod( (ms-1)/MBLK,NPROW )

          utmp3=zero
          if ( irank_b == myrank_b ) utmp3(1:ii,1:mm)=u(i1:i2,ms:me)

          call mpi_bcast(utmp3,ii*mm,TYPE_MAIN,irank_b,comm_band,ierr)

          j0=0
          jj0=0
          do ib=1,(ncycle-1)/np_band+1
             nbss=(ib-1)*np_band+myrank_b+1
             ns=MBLK*(nbss-1)+1
             ne=min(ns+MBLK-1,MB)
             nn=ne-ns+1

             IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)

             if ( mm < 1 .or. nn < 1 ) cycle

             irank_c = mod( (ns-1)/MBLK, np_band )
             if ( irank_c == myrank_b ) then
             
                allocate( utmp2(ms:me,ns:ne) )

                if ( iroot1 == myrank ) then
!$OMP parallel workshare
                   utmp2(ms:me,ns:ne)=Vsub(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
                   j0=j0+nn
                end if

                call mpi_bcast(utmp2,mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)

                if ( ii > 0 ) then
                   call dgemm('N','N',ii,nn,mm,one,utmp3 &
                        ,ML0,utmp2(ms,ns),mm,one,utmp4(1,jj0+1),NBLK2)
                   jj0=jj0+nn
                end if

                deallocate( utmp2 )

             endif

          end do ! ms

          if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

       end do ! ns

! update u(:,:)
       if ( ii > 0 ) then
          jj0=0
          do ib=1,(ncycle-1)/np_band+1
             nbss= (ib-1)*np_band+myrank_b+1
             ns=MBLK*(nbss-1)+1
             ne=min(ns+MBLK-1,MB)
             nn=ne-ns+1
             u(i1:i2,ns:ne)=utmp4(1:ii,jj0+1:jj0+nn)
             jj0=jj0+nn
          enddo
       end if

    end do ! ii

    deallocate( utmp3 )
    deallocate( utmp4 )

    call write_border( 1, " subspace_rotv_sl_bp1(end)" )

  END SUBROUTINE subspace_rotv_sl_bp1

  SUBROUTINE subspace_rotv_sl_bp2( u, v )

    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: i,i1,i2,ii,n1,n2,i0,j0,ns,ne,nn,ms,me,mm,MB_0,MB_1
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,MB,comm_grid
    integer :: myrank
    include 'mpif.h'
    real(8),allocatable :: utmp3(:,:),utmp4(:,:)
    integer :: comm_band, np_band, irank_b, MBLK, myrank_b
    integer :: irank_c, ncycle, icycle, ib, nbss
    integer :: MB2, jj0
    integer :: msV, meV, nsV, neV

    call write_border( 1, " subspace_rotv_sl_bp2(start)" )

    ML0  = size( u, 1 )
    MB   = MB_diag
    MB_0 = v(2)%pinfo%id(v(2)%pinfo%me) + 1
    MB_1 = MB_0 + v(2)%pinfo%ir(v(2)%pinfo%me) - 1
    comm_grid = v(1)%pinfo%comm

    np_band   = v(2)%pinfo%np
    myrank_b  = v(2)%pinfo%me
    comm_band = v(2)%pinfo%comm
    MBLK=min(MBSIZE,NBSIZE)

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    NBLK2 = maxval( v(1)%pinfo%ir )

    allocate( utmp3(NBLK2,MBLK) ); utmp3=zero

    MB2=(MB-1)/np_band+1
    allocate( utmp4(NBLK2,MB2) ); utmp4=zero

    ncycle=(MB-1)/MBLK+1

    do i1=1,maxval( v(1)%pinfo%ir ),NBLK2

       i2=min(i1+NBLK2-1,ML0)
       ii=i2-i1+1

       utmp4(:,:)=zero

       i0=0
       do icycle=1,ncycle
          ms=MBLK*(icycle-1)+1
          me=min(ms+MBLK-1,MB)
          mm=me-ms+1

          msV=MBLK*int((icycle-1)/np_band )+ 1
          meV=min(msV+MBLK-1,MB)

          irank_b = mod( icycle-1, np_band )

          IPROW=mod( (ms-1)/MBLK,NPROW )

          utmp3=zero
          if ( irank_b == myrank_b ) utmp3(1:ii,1:mm)=u(i1:i2,msV:meV)

          call mpi_bcast(utmp3,ii*mm,TYPE_MAIN,irank_b,comm_band,ierr)

          j0=0
          jj0=0
          do ib=1,(ncycle-1)/np_band+1
             nbss=(ib-1)*np_band+myrank_b+1
             ns=MBLK*(nbss-1)+1
             ne=min(ns+MBLK-1,MB)
             nn=ne-ns+1

             nsV = MBLK*(ib-1)+1
             neV = min(nsV+MBLK-1,MB)

             IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)

             if ( mm < 1 .or. nn < 1 ) cycle

             irank_c = mod( (ns-1)/MBLK, np_band )
             if (irank_c==myrank_b) then
             
                allocate( utmp2(ms:me,ns:ne) )

                if ( iroot1 == myrank ) then
!$OMP parallel workshare
                   utmp2(ms:me,ns:ne)=Vsub(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
                   j0=j0+nn
                end if

                call mpi_bcast(utmp2,mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr)

                if ( ii > 0 ) then
                   call dgemm('N','N',ii,nn,mm,one,utmp3 &
                        ,ML0,utmp2(ms,ns),mm,one,utmp4(1,jj0+1),NBLK2)
                   jj0=jj0+nn
                end if

                deallocate( utmp2 )

             endif

          end do ! ms

          if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

       end do ! ns

! update u(:,:)
       if ( ii > 0 ) then
          jj0=0
          do ib=1,(ncycle-1)/np_band+1
             nbss= (ib-1)*np_band+myrank_b+1
             ns=MBLK*(nbss-1)+1
             ne=min(ns+MBLK-1,MB)
             nn=ne-ns+1

             nsV = MBLK*(ib-1)+1
             neV = min(nsV+MBLK-1,MB)

             u(i1:i2,nsV:neV)=utmp4(1:ii,jj0+1:jj0+nn)
             jj0=jj0+nn
          enddo
       end if

    end do ! ii

    deallocate( utmp3 )
    deallocate( utmp4 )

    call write_border( 1, " subspace_rotv_sl_bp2(end)" )

  END SUBROUTINE subspace_rotv_sl_bp2

END MODULE subspace_rotv_sl_module
