!--------------------------------------------------------------------------
! Copyright 2016 Junichi Iwata
! 
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
! 
!     http://www.apache.org/licenses/LICENSE-2.0
! 
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!--------------------------------------------------------------------------

MODULE gram_schmidt_pbp_module

  use io_tools_module
  use vector_tools_module, only: normalize_vector, inner_product_vector, vinfo
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none

  PRIVATE
  PUBLIC :: read_gram_schmidt_pbp
  PUBLIC :: gram_schmidt_pbp

  integer :: NBLK=0
  integer :: NBLK1=0

  integer :: dunit=6

  integer :: JBUFNUM=0

  INTERFACE gram_schmidt_pbp
     MODULE PROCEDURE d_gram_schmidt_pbp
  END INTERFACE

CONTAINS


  SUBROUTINE read_gram_schmidt_pbp( NBLK_OUT )
    implicit none
    integer,optional,intent(OUT) :: NBLK_OUT
    call IOTools_readIntegerKeyword( "NBLK", NBLK )
    call IOTools_readIntegerKeyword( "NBLK1", NBLK1 )
    call IOTools_readIntegerKeyword( "JBUFNUM", JBUFNUM )
    if ( present(NBLK_OUT) ) NBLK_OUT=NBLK
  END SUBROUTINE read_gram_schmidt_pbp


  SUBROUTINE d_gram_schmidt_pbp( u, v )
    implicit none
    include 'mpif.h'
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: irank_b,np_band,myrank_b,ns,ne,ms,me
    integer :: nbss,ib,NBAND_BLK,ncycle,MB
    integer i, comm_count, recved_flag
    real(8),allocatable :: utmp(:,:,:)
    integer :: nsV, neV, msV, meV, nnV, ML0
    integer,allocatable:: index_array(:), send_req(:), status(:,:)
    logical :: lflag, recv_req_lflag
    integer :: recv_req, ierr
    integer ::            jbuf_index
    integer,allocatable:: jbuf_flag(:)

    call write_border( 1, " d_gram_schmidt_pbp(start)" )

    MB       = size( u, 2 )
    np_band  = v(2)%pinfo%np
    myrank_b = v(2)%pinfo%me

    MB       = MB*np_band

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4
    if ( JBUFNUM == 0 ) JBUFNUM = 3
!   JBUFNUM = min(JBUFNUM,ncycle)

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK+1
    JBUFNUM = min(JBUFNUM,ncycle)

    comm_count = 0
    allocate(index_array(0:ncycle-1))
    allocate(send_req(0:np_band-1))
    allocate(status(MPI_STATUS_SIZE,0:np_band-1))
    send_req = MPI_REQUEST_NULL
    index_array = 0
    recv_req_lflag = .false.

    allocate (jbuf_flag(JBUFNUM))
    jbuf_flag(:) = 0
    ML0 = size( u, 1 )
    allocate (utmp(ML0,NBAND_BLK,JBUFNUM))
    do
       ! 対角計算およびデータ送信
       do ib=myrank_b, ncycle-1, np_band
          if(index_array(ib) .eq. ib) then
             jbuf_index = mod(comm_count,JBUFNUM)+1
             if(jbuf_flag(jbuf_index) .eq. 0) then
                ns  = NBAND_BLK*index_array(ib) + 1
                ne  = min( ns+NBAND_BLK-1, MB )
                nsV = int((ns-1)/(np_band*NBAND_BLK))*NBAND_BLK+1
                neV = min( nsV+NBAND_BLK-1, MB )
                nnV = neV-nsV+1
                call d_gram_schmidt_t_sub( u(:,nsV:neV), u(:,nsV:neV),v(1), &
                                      ns,ne, ns,ne, NBLK )
                call MPI_Waitall(np_band, send_req, status, ierr)
                comm_count = comm_count+1
                do i=0, np_band-1
                   if(i .eq. myrank_b) cycle
                   call MPI_Isend(u(:,nsV:neV),ML0*nnV, MPI_REAL8, i, &
                        100, v(2)%pinfo%comm, send_req(i), ierr)
                enddo
                index_array(ib) = index_array(ib)+1
                if((comm_count-1) .lt. (ncycle-np_band+myrank_b)) then
                   jbuf_flag(jbuf_index) = 1
                endif
                goto 100
             endif
          endif
       enddo

       ! データ受信
       irank_b = mod(comm_count, np_band)
       if(myrank_b .ne. irank_b) then
          if(.not.recv_req_lflag) then
             jbuf_index = mod(comm_count,JBUFNUM)+1
             if(jbuf_flag(jbuf_index) .eq. 0) then
                ns  = NBAND_BLK*(comm_count) + 1
                ne  = min( ns+NBAND_BLK-1, MB )
                nsV = int((ns-1)/(np_band*NBAND_BLK))*NBAND_BLK+1
                neV = min( nsV+NBAND_BLK-1, MB )
                nnV = neV-nsV+1
                call MPI_IRecv(utmp(:,:,jbuf_index),ML0*nnV, MPI_REAL8,irank_b, &
                     100, v(2)%pinfo%comm, recv_req, ierr)
                recv_req_lflag = .true.
                recved_flag = 0
                lflag = .false.
                if(comm_count .lt. (ncycle-np_band+myrank_b)) then
                   jbuf_flag(jbuf_index) = 1
                endif
             endif
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
             nbss = ib + 1
             if ( index_array(ib)+1 <= ib .and. ib <= ncycle-1 ) then
                ns  = NBAND_BLK*index_array(ib) + 1
                ne  = min( ns+NBAND_BLK-1, MB )
                ms  = NBAND_BLK*ib+1
                me  = min( ms+NBAND_BLK-1, MB )
                msV = int((ms-1)/(NBAND_BLK*np_band))*NBAND_BLK+1
                meV = min( msV+NBAND_BLK-1, MB )
                jbuf_index = mod(index_array(ib),JBUFNUM)+1
                if ( ms <= me ) then
                   if(mod(index_array(ib),np_band) .eq. myrank_b) then
                      nsV = int((ns-1)/(np_band*NBAND_BLK))*NBAND_BLK+1
                      neV = min( nsV+NBAND_BLK-1, MB )
                      call d_gram_schmidt_t_sub( u(:,msV:meV), u(:,nsV:neV), v(1), &
                                        ms,me, ns,ne, NBLK )
                   else
                      call d_gram_schmidt_t_sub( u(:,msV:meV), utmp(:,:,jbuf_index), v(1), &
                                        ms,me, ns,ne, NBLK )
                   endif
                   index_array(ib) = index_array(ib)+1
                   if(ib .gt. (ncycle-1-np_band)) then 
                     jbuf_flag(jbuf_index) = 0
                   endif
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

    deallocate (utmp) 
    deallocate (jbuf_flag) 

    call write_border( 1, " d_gram_schmidt_pbp(end)" )

    return

  END SUBROUTINE d_gram_schmidt_pbp

! RECURSIVE SUBROUTINE d_gram_schmidt_t_sub( u, v, mm1,mm2, nn1,nn2, MBLK )
  RECURSIVE SUBROUTINE d_gram_schmidt_t_sub( u, u2, v, mm1,mm2, nn1,nn2, MBLK )
    implicit none
    real(8),intent(INOUT) :: u(:,:)
    real(8),intent(IN)    :: u2(:,:)
!   integer, intent(IN)   :: icyc1, icyc2
    type(vinfo),intent(IN) :: v
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK
    integer :: n,ns,ne,m,ms,me,mm,nn,MBLKH,ML0
    real(8),allocatable :: utmp1(:), utmp2(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
!---
!-- MPI
    include 'mpif.h'
    integer :: MB,np_band,myrank_b,ierr
    integer :: nsR,neR, msR, meR
!---
    integer :: mm1mm2, nn1nn2

    call write_border( 1, " d_gram_schmidt_t_sub(start)" )

!----
!   call mpi_comm_size(mpi_comm_world,np_band,ierr)
!   call mpi_comm_rank(mpi_comm_world,myrank_b,ierr)
    ML0 = size( u, 1 )
!   MB = np_band*size( u, 2 )
!----

    mm1mm2=mm2-mm1+1
!   do ms=mm1,mm2,MBLK
    do ms=1,mm1mm2,MBLK

!      me=min(ms+MBLK-1,mm2)
       me=min(ms+MBLK-1,mm1mm2)
       mm=me-ms+1

!---
       msR = mm1+ms-1
       meR = mm1+me-1
!---

       nn1nn2=nn2-nn1+1
!      do ns=nn1,nn2,MBLK
       do ns=1,nn1nn2,MBLK


          ne=min(ns+MBLK-1,nn1nn2)
!         ne=min(ns+MBLK-1,nn2)
!----------------------------------------------------
!         ne=min(ne,me-1)
          if ( (nn1+ne-1) > (mm1+me-1)-1 ) then 
             ne=min(ne,me-1)
          end if
!----------------------------------------------------
          nn=ne-ns+1
!---
          nsR = nn1+ns-1
          neR = nn1+ne-1
!---

          if ( nn <= 0 ) cycle

!         if ( ms >= ne+1 ) then  
          if ( msR >= neR+1 ) then  

             allocate( utmp2(ns:ne,ms:me) ) ; utmp2=zero

!            call inner_product_vector( u(:,ns:ne), u(:,ms:me), utmp2, v )
             call inner_product_vector( u2(:,ns:ne), u(:,ms:me), utmp2, v )

!            call dgemm('N','N',ML0,mm,nn,-one,u(:,ns) &
!                       ,ML0,utmp2,nn,one,u(:,ms),ML0)
             call dgemm('N','N',ML0,mm,nn,-one,u2(:,ns) &
                        ,ML0,utmp2,nn,one,u(:,ms),ML0)

             deallocate( utmp2 )

!            if ( ms == ne+1 ) then
             if ( msR == neR+1 ) then  
                call normalize_vector( u(:,ms), v )
             end if

          else if ( mm <= NBLK1 ) then

             allocate( utmp1(NBLK1) ) ; utmp1=zero

             do m=ms,me
                n=min(m-1,ne)

                if ( n-ns+1 > 0 ) then

!                  call inner_product_vector( u(:,ns:n), u(:,m), utmp1, v )
                   call inner_product_vector( u2(:,ns:n), u(:,m), utmp1, v )

!                  call dgemv('N',ML0,n-ns+1,-one,u(:,ns) &
                   call dgemv('N',ML0,n-ns+1,-one,u2(:,ns) &
                              ,ML0,utmp1,1,one,u(:,m),1)

                end if

                if ( m==1 .or. ( n==m-1 .and. m/=ns ) ) then
                   call normalize_vector( u(:,m), v )
                end if

             end do ! m

             deallocate( utmp1 )

          else

             MBLKH=max(MBLK/2,NBLK1)

!            call d_gram_schmidt_t_sub( u, v, ms,me, ns,ne, MBLKH )
             call d_gram_schmidt_t_sub( u(:,ms:me),u(:,ns:ne),v,msR,meR,nsR,neR,MBLKH )

          end if

       end do ! ns

    end do ! ms

    call write_border( 1, " d_gram_schmidt_t_sub(end)" )

    return

  END SUBROUTINE d_gram_schmidt_t_sub

!-----------------------------------------------------------------------

END MODULE gram_schmidt_pbp_module
