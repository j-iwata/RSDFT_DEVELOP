MODULE calc_overlap_mbp_module

  use parallel_module , only: comm_grid, comm_band, myrank, myrank_b, np_band
  use rsdft_mpi_module, only: rsdft_allreduce_sum
  use scalapack_module, only: usermap, NPROW, NPCOL, MBSIZE, NBSIZE
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none
  include 'mpif.h'

  PRIVATE
  PUBLIC :: calc_overlap_mbp

  integer :: MBLK0
  integer :: MBLK1
  real(8),parameter :: alpha=1.d0,beta=0.d0

  integer,allocatable :: ir(:),id(:),irecv_me(:,:),isend_me(:,:)
  integer :: nsend_me, nrecv_me

CONTAINS

  SUBROUTINE calc_overlap_mbp(m,n,nb,a,b,dv,c,iflag,ii0,jj0)
    implicit none
    integer,intent(IN)  :: m,n,nb
!   real(8),intent(IN)  :: a(m,n),b(m,n),dv
    real(8),intent(IN)  :: a(m,nb),b(m,nb),dv
!   real(8),intent(OUT) :: c(n,n)
    real(8),intent(OUT) :: c(:,:)
    integer,intent(OUT) :: ii0, jj0
    integer :: nme,i,j,k,ierr
    real(8),allocatable :: s(:),r(:)
  
!   include 'mpif.h'

!   integer :: i0,i1,j0,j1,ni,nj,MBLKH,MBLK
    integer :: MBLKH,MBLK
    real(8),allocatable :: ctmp(:,:)
  
    integer :: IPROW, IPCOL, iroot1, iroot2
    integer :: ns, ne, ms, me, nn, mm, i0, j0
    real(8),allocatable :: vtmp2(:,:),wtmp2(:,:)
    logical :: iflag
    logical,parameter :: iswitch_sm_cp=.true.
!   logical,parameter :: iswitch_sm_cp=.false.


!
! block-cycle, band-parallel
!________________________________________________________________________________
    integer :: ncycle, ncycle2, icycle, icycle2, irankb, nbss
    integer :: msV, meV, nsV, neV
    real(8), allocatable :: aa(:,:)
    real(8),parameter :: zero=0.0d0

    MBLK=max(MBSIZE,NBSIZE)

    ncycle  = (n-1)/MBLK+1
    ncycle2 = (ncycle-1)/np_band+1
    allocate ( aa(m,MBLK) ) ; aa=zero

    MBLK0 = n
    MBLK1 = 4

    ii0=0
    jj0=0

    call init_srdata

!   call calc_overlap_sub(m,n,MBLK0,a,b,c)

    i0=0
!   do ms=1,n,MBLK        ! Row-loop
!      me=min(ms+MBLK-1,n)
!      mm=me-ms+1
    do icycle=1,ncycle                   ! Row-loop
       irankb = mod( icycle-1, np_band) 
       ms = MBLK*(icycle-1)+1
       me = min(ms+MBLK-1,n)
       mm = me-ms+1
       if (mm<=0) cycle
       IPROW = mod( (ms-1)/MBLK, NPROW)
       msV = MBLK*int((icycle-1)/np_band) + 1
       meV = min( msV+MBLK-1, n )

!----  Brodcast ( irankb => other process in comm_band)
       if (irankb==myrank_b) aa(:,:) = a(:,msV:meV)
       call rsdft_bcast ( aa, irankb, comm_band )

       j0=0
       do icycle2=1, ncycle2      ! Column-loop
          nbss = (icycle2-1)*np_band + myrank_b + 1  
          ns = MBLK*(nbss-1)+1
          ne = min(ns+MBLK-1,n)
          nn=ne-ns+1
          if (nn<=0) cycle

          nsV = MBLK*(icycle2-1)+1
          neV = min(nsV+MBLK-1,n)

          IPCOL = mod( (ns-1)/MBLK, NPCOL)
          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)

          if ( ns > me ) then

!!!          cycle
!!           call dgemm ('T','N',mm,nn,m,alpha,a(1,ms),m,b(1,ns),m,beta,vtmp2(ms,ns),mm)

             if (iswitch_sm_cp) then
                if ( iroot1 == myrank ) then
                    call set_recvdata( ms,me,ns,ne,i0+1,i0+mm,j0+1,j0+nn )
                    j0=j0+nn
                end if
             else 
                allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) ); vtmp2=0.d0; wtmp2=0.d0
                call dgemm ('T','N',mm,nn,m,alpha,aa,m,b(1,nsV),m,beta,vtmp2(ms,ns),mm)
                call mpi_reduce(vtmp2,wtmp2,mm*nn,MPI_REAL8,MPI_SUM,iroot2,comm_grid,ierr)
                if ( iroot1 == myrank ) then
                    c(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)*dv
                    j0=j0+nn
                end if
                deallocate( wtmp2,vtmp2 )
              endif
          else if ( ne <= ms ) then
 
             allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) ); vtmp2=0.d0; wtmp2=0.d0
!!           call dgemm &
!!                ('T','N',mm,nn,m,alpha,a(1,ms),m,b(1,ns),m,beta,c(ms,ns),n)
!!           call dgemm ('T','N',mm,nn,m,alpha,a(1,ms),m,b(1,ns),m,beta,vtmp2(ms,ns),mm)
             call dgemm ('T','N',mm,nn,m,alpha,aa,m,b(1,nsV),m,beta,vtmp2(ms,ns),mm)
             call mpi_reduce(vtmp2,wtmp2,mm*nn,MPI_REAL8,MPI_SUM,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then

                 if (iswitch_sm_cp ) call set_senddata( ms,me,ns,ne,i0+1,i0+mm,j0+1,j0+nn )

!                c(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                 c(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)*dv
                 j0=j0+nn
             end if
             deallocate( wtmp2,vtmp2 )

          else

             allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) ); vtmp2=0.d0; wtmp2=0.d0

             if ( mm > MBLK1 ) then
!!              write(*,*) " # overlap_sub :myrank,ms,ns",myrank,ms,ns
                MBLKH=MBLK/2
!!              call calc_overlap_sub(m,mm,MBLKH,a(1,ms),b(1,ns),vtmp2)
                call calc_overlap_mbp_sub(m,mm,MBLKH,aa,b(1,nsV),vtmp2)
             else
!               write(*,*) " # mvect :myrank,ms,ns",myrank,ms,ns,msV,nsV,icycle,icycle2
                do j=ns,ne
                   do i=j ,me
!!                    vtmp2(i,j)=sum( a(:,i)*b(:,j) )
                      vtmp2(i,j)=sum( aa(:,i-ms+1)*b(:,nsV+j-ns) )
                   end do
                end do

             end if
!
! reduce comm_grid
!
             call mpi_reduce(vtmp2,wtmp2,mm*nn,MPI_REAL8,MPI_SUM,iroot2,comm_grid,ierr)
             if ( iroot1 == myrank ) then
!!               c(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
!!               c(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)*dv
!!               j0=j0+nn
                 if ( iflag ) then
                    do j=ns,ne
                       do i=j,me
                         if (i/=j) then
                            c(i0+i-ms+1,j0+j-ns+1)=wtmp2(i,j)*dv
! 2020/01/28                c(j0+j-ns+1,i0+i-ms+1)=wtmp2(i,j)*dv
                            c(i0+j-ns+1,j0+i-ms+1)=wtmp2(i,j)*dv
                         else
                            c(i0+i-ms+1,j0+j-ns+1)=wtmp2(i,j)*dv + 1.0d0
!                           c(i0+i-ms+1,j0+j-ns+1)=wtmp2(i,j)*dv
                         endif
                       enddo
                    enddo
                 else
                    do j=ns,ne
                       do i=j,me
                          c(i0+i-ms+1,j0+j-ns+1)=wtmp2(i,j)*dv
! 2020/01/28              c(j0+j-ns+1,i0+i-ms+1)=wtmp2(i,j)*dv
                          c(i0+j-ns+1,j0+i-ms+1)=wtmp2(i,j)*dv
                       enddo
                    enddo
                 endif
                 j0=j0+nn
             end if
             deallocate( wtmp2,vtmp2 )

          endif
          jj0=max(j0,jj0)
       enddo
       if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm
       ii0=max(i0,ii0)
    enddo

!--------------------------------------------------------------------------------

!    exchange symetrical data ( B_n,m <---> B_m,n )
    call exchange_srdata(c(:,:))
    call finalize_srdata

    return

  END SUBROUTINE calc_overlap_mbp

  RECURSIVE SUBROUTINE calc_overlap_mbp_sub(m,n,MBLK,a,b,c)
    implicit none
    integer,intent(IN)  :: m,n,MBLK
    real(8),intent(IN)  :: a(m,n),b(m,n)
    real(8),intent(OUT) :: c(n,n)
    integer :: i,j,i0,i1,j0,j1,ni,nj,MBLKH
    real(8),allocatable :: ctmp(:,:)

    do i0=1,n,MBLK
       i1=min(i0+MBLK-1,n)
       ni=i1-i0+1

       do j0=1,n,MBLK
          j1=min(j0+MBLK-1,n)
          nj=j1-j0+1

          if ( j0 > i1 ) then

             cycle

          else if ( j1 <= i0 ) then

             call dgemm &
                  ('T','N',ni,nj,m,alpha,a(1,i0),m,b(1,j0),m,beta,c(i0,j0),n)

          else

             if ( ni > MBLK1 ) then
                allocate( ctmp(ni,ni) ) ; ctmp=0.d0
                MBLKH=MBLK/2
                call calc_overlap_mbp_sub(m,ni,MBLKH,a(1,i0),b(1,j0),ctmp)
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
                do j=j0,j1
                do i=j ,i1
                    c(i,j)=sum( a(:,i)*b(:,j) )
                end do
                end do
             end if

          end if

       end do ! j0

    end do ! i0

  END SUBROUTINE calc_overlap_mbp_sub




  SUBROUTINE init_srdata
    implicit none

    allocate( irecv_me(99,0:8),isend_me(99,0:8) )

    nrecv_me      = 0
    nsend_me      = 0
    irecv_me(:,:) =-1
    isend_me(:,:) =-1

    return
  END SUBROUTINE init_srdata

  SUBROUTINE set_senddata( ms,me,ns,ne,i1,i2,j1,j2 )
    implicit none
    integer, intent(IN) :: ms,me,ns,ne,i1,i2,j1,j2
    integer :: i, j, n

    i=mod( (ns-1)/MBSIZE, NPROW )
    j=mod( (ms-1)/NBSIZE, NPCOL )
    n=usermap(i,j,1)

    nsend_me=nsend_me+1
    isend_me(nsend_me,0)=n
    isend_me(nsend_me,1)=ms
    isend_me(nsend_me,2)=me
    isend_me(nsend_me,3)=ns
    isend_me(nsend_me,4)=ne
    isend_me(nsend_me,5)=i1         !i0+1
    isend_me(nsend_me,6)=i2         !i0+mm
    isend_me(nsend_me,7)=j1         !j0+1
    isend_me(nsend_me,8)=j2         !j0+nn

    return
  END SUBROUTINE set_senddata

  SUBROUTINE set_recvdata( ms,me,ns,ne,i1,i2,j1,j2 )
    implicit none
    integer, intent(IN) :: ms,me,ns,ne,i1,i2,j1,j2
    integer :: i, j, n

    i=mod( (ns-1)/MBSIZE, NPROW )
    j=mod( (ms-1)/NBSIZE, NPCOL )
    n=usermap(i,j,1)

    nrecv_me=nrecv_me+1
    irecv_me(nrecv_me,0)=n
    irecv_me(nrecv_me,1)=ms
    irecv_me(nrecv_me,2)=me
    irecv_me(nrecv_me,3)=ns
    irecv_me(nrecv_me,4)=ne
    irecv_me(nrecv_me,5)=i1         !i0+1
    irecv_me(nrecv_me,6)=i2         !i0+mm
    irecv_me(nrecv_me,7)=j1         !j0+1
    irecv_me(nrecv_me,8)=j2         !j0+nn

  END SUBROUTINE set_recvdata

  SUBROUTINE exchange_srdata(Hsub)
    implicit none
    real(8), intent(INOUT) :: Hsub(:,:)
    integer :: i, j, n, m, i0, j0
    real(8),allocatable :: vtmp2(:,:), wtmp2(:,:)
    integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer ::  ms, me, ns, ne, i1, i2, j1, j2

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
    endif

  END SUBROUTINE exchange_srdata

  SUBROUTINE finalize_srdata
    deallocate( irecv_me,isend_me )
  END SUBROUTINE finalize_srdata


END MODULE calc_overlap_mbp_module
