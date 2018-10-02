!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

MODULE gram_schmidt_p_module

  use io_tools_module
  use vector_tools_module, only: normalize_vector, inner_product_vector, vinfo
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none

  PRIVATE
  PUBLIC :: read_gram_schmidt_p
  PUBLIC :: gram_schmidt_p

  integer :: NBLK=0
  integer :: NBLK1=0

  INTERFACE gram_schmidt_p
     MODULE PROCEDURE d_gram_schmidt_p
  END INTERFACE

CONTAINS


  SUBROUTINE read_gram_schmidt_p
    implicit none
    call IOTools_readIntegerKeyword( "NBLK", NBLK )
    call IOTools_readIntegerKeyword( "NBLK1", NBLK1 )
  END SUBROUTINE read_gram_schmidt_p


  SUBROUTINE d_gram_schmidt_p( u, v )
    implicit none
    include 'mpif.h'
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: irank_b,np_band,myrank_b,ns,ne,ms,me
    integer :: nbss,ib,NBAND_BLK,ncycle,icycle,MB
    integer i, comm_count, recved_flag
    integer,allocatable:: index_array(:), send_req(:), status(:,:)
    logical :: lflag, recv_req_lflag
    integer :: recv_req, ierr

    call write_border( 1, " d_gram_schmidt_p(start)" )

    MB       = size( u, 2 )
    np_band  = v(2)%pinfo%np
    myrank_b = v(2)%pinfo%me

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK+1

    comm_count = 0
    allocate(index_array(0:ncycle-1))
    allocate(send_req(0:np_band-1))
    allocate(status(MPI_STATUS_SIZE,0:np_band-1))
    send_req = MPI_REQUEST_NULL
    index_array = 0
    recv_req_lflag = .false.

    do
       ! 対角計算およびデータ送信
       do ib=myrank_b, ncycle-1, np_band
          if(index_array(ib) .eq. ib) then
             ns = NBAND_BLK*index_array(ib) + 1
             ne = min( ns+NBAND_BLK-1, MB )
#ifdef _TIMER_
             call timer_sta(2)
#endif
             call d_gram_schmidt_p_sub( u, v(1), ns,ne, ns,ne, NBLK )
#ifdef _TIMER_
             call timer_end(2)
#endif
             call MPI_Waitall(np_band, send_req, status, ierr)
             do i=0, np_band-1
                if(i.eq.myrank_b) cycle
                call MPI_Isend(u(:,ns:ne),size(u,1)*(ne-ns+1), MPI_REAL8, i, 100, v(2)%pinfo%comm, send_req(i), ierr)
             enddo
             index_array(ib) = index_array(ib)+1
             comm_count = comm_count+1
             goto 100
          endif
       enddo

       ! データ受信
       irank_b = mod(comm_count, np_band)
       if(myrank_b .ne. irank_b) then
          if(.not.recv_req_lflag) then
             ns = NBAND_BLK*(comm_count) + 1
             ne = min( ns+NBAND_BLK-1, MB )
             call MPI_IRecv(u(:,ns:ne),size(u,1)*(ne-ns+1), MPI_REAL8,irank_b, &
                  100, v(2)%pinfo%comm, recv_req, ierr)
             recv_req_lflag = .true.
             recved_flag = 0
             lflag = .false.
          endif
          if(recv_req_lflag) then
             if (recved_flag==0) then
                call MPI_Test(recv_req, recved_flag, status(1,myrank_b), ierr)
             endif
             if (recved_flag==1) then
                lflag = .true.
             else
                lflag = .false.
             endif
             call MPI_Allreduce(MPI_IN_PLACE, lflag, 1, mpi_logical, mpi_land, v(1)%pinfo%comm, ierr)
             if (lflag) then
                comm_count = comm_count+1
                recv_req_lflag = .false.
             endif
          endif
       endif

       ! 行列計算
       do ib=myrank_b, ncycle-1, np_band
          if(index_array(ib) .lt. comm_count) then
             ns = NBAND_BLK*index_array(ib) + 1
             ne = min( ns+NBAND_BLK-1, MB )

!            nbss = ib*np_band + myrank_b + 1
             nbss = ib + 1

             if ( nbss <= ncycle .and. nbss >= index_array(ib)+2 ) then

                ms = NBAND_BLK*(nbss-1)+1
                me = min(ms+NBAND_BLK-1,MB)
                if ( ms <= me ) then
#ifdef _TIMER_
                   call timer_sta(3)
#endif
                   call d_gram_schmidt_p_sub( u, v(1), ms,me, ns,ne, NBLK )
#ifdef _TIMER_
                   call timer_end(3)
#endif
                   index_array(ib) = index_array(ib)+1
                   goto 100
                end if
             endif
          endif
       enddo

100    continue
       ! 終了判定
       if(comm_count .eq. ncycle) exit
    enddo
    call MPI_Waitall(np_band, send_req, status, ierr)
    deallocate(index_array, send_req, status)

    call write_border( 1, " d_gram_schmidt_p(end)" )

    return

  END SUBROUTINE d_gram_schmidt_p

  RECURSIVE SUBROUTINE d_gram_schmidt_p_sub( u, v, mm1,mm2, nn1,nn2, MBLK )
    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK
    integer :: n,ns,ne,m,ms,me,mm,nn,MBLKH,ML0
    real(8),allocatable :: utmp1(:), utmp2(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
!#ifndef _TIMER_
!    call write_border( 1, " d_gram_schmidt_p_sub(start)" )
!#endif
    ML0 = size( u, 1 )

    do ms=mm1,mm2,MBLK

       me=min(ms+MBLK-1,mm2)
       mm=me-ms+1

       do ns=nn1,nn2,MBLK

          ne=min(ns+MBLK-1,nn2)
          ne=min(ne,me-1)
          nn=ne-ns+1

          if ( nn <= 0 ) cycle

          if ( ms >= ne+1 ) then

             allocate( utmp2(ns:ne,ms:me) ) ; utmp2=zero

             call inner_product_vector( u(:,ns:ne), u(:,ms:me), utmp2, v )

             call dgemm('N','N',ML0,mm,nn,-one,u(:,ns) &
                        ,ML0,utmp2,nn,one,u(:,ms),ML0)

             deallocate( utmp2 )

             if ( ms == ne+1 ) then

                call normalize_vector( u(:,ms), v )

             end if

          else if ( mm <= NBLK1 ) then

             allocate( utmp1(NBLK1) ) ; utmp1=zero

             do m=ms,me

                n=min(m-1,ne)

                if ( n-ns+1 > 0 ) then

                   call inner_product_vector( u(:,ns:n), u(:,m), utmp1, v )

                   call dgemv('N',ML0,n-ns+1,-one,u(:,ns) &
                              ,ML0,utmp1,1,one,u(:,m),1)

                end if

                if ( m==1 .or. ( n==m-1 .and. m/=ns ) ) then

                   call normalize_vector( u(:,m), v )

                end if

             end do ! m

             deallocate( utmp1 )

          else

             MBLKH=max(MBLK/2,NBLK1)
             call d_gram_schmidt_p_sub( u, v, ms,me, ns,ne, MBLKH )

          end if

       end do ! ns

    end do ! ms
!#ifndef _TIMER_
!    call write_border( 1, " d_gram_schmidt_p_sub(end)" )
!#endif
    return

  END SUBROUTINE d_gram_schmidt_p_sub

END MODULE gram_schmidt_p_module
