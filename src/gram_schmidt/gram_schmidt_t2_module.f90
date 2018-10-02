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

MODULE gram_schmidt_t2_module

  use io_tools_module
  use vector_tools_module, only: normalize_vector, inner_product_vector, vinfo
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none

  PRIVATE
  PUBLIC :: read_gram_schmidt_t2
  PUBLIC :: gram_schmidt_t2

  integer :: NBLK=0
  integer :: NBLK1=0

  INTERFACE gram_schmidt_t2
     MODULE PROCEDURE d_gram_schmidt_t, z_gram_schmidt_t
  END INTERFACE

CONTAINS


  SUBROUTINE read_gram_schmidt_t2
    implicit none
    call IOTools_readIntegerKeyword( "NBLK", NBLK )
    call IOTools_readIntegerKeyword( "NBLK1", NBLK1 )
  END SUBROUTINE read_gram_schmidt_t2


  SUBROUTINE d_gram_schmidt_t( u, v )
    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: irank_b,np_band,myrank_b,ns,ne,ms,me
    integer :: nbss,ib,NBAND_BLK,ncycle,icycle,MB

    call write_border( 1, " d_gram_schmidt_t2(start)" )

    MB       = size( u, 2 )
    np_band  = v(2)%pinfo%np
    myrank_b = v(2)%pinfo%me

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK+1

    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ns = NBAND_BLK*(icycle-1) + 1
       ne = min( ns+NBAND_BLK-1, MB )

       if ( irank_b == myrank_b ) then
#ifdef _TIMER_
!         call MPI_Barrier(v(1)%pinfo%comm, ierr)
          call timer_sta(2)
#endif
          call d_gram_schmidt_t_sub( u, v(1), ns,ne, ns,ne, NBLK )
#ifdef _TIMER_
!         call MPI_Barrier(v(1)%pinfo%comm, ierr)
          call timer_end(2)
#endif
       end if

       call rsdft_bcast( u(:,ns:ne), irank_b, v(2)%pinfo%comm )

       if ( ns <= MB-NBAND_BLK ) then

          do ib=1,(ncycle-1)/np_band+1

             nbss = (ib-1)*np_band + myrank_b + 1

             if ( nbss <= ncycle .and. nbss >= icycle+1 ) then

                ms = NBAND_BLK*(nbss-1)+1
                me = min(ms+NBAND_BLK-1,MB)

                if ( ms <= me ) then
#ifdef _TIMER_
                   call timer_sta(3)
#endif
                   call d_gram_schmidt_t_sub( u, v(1), ms,me, ns,ne, NBLK )
#ifdef _TIMER_
                   call timer_end(3)
#endif
                end if

             end if

          end do

       end if

    end do ! icycle

    call write_border( 1, " d_gram_schmidt_t2(end)" )

    return

  END SUBROUTINE d_gram_schmidt_t

  RECURSIVE SUBROUTINE d_gram_schmidt_t_sub( u, v, mm1,mm2, nn1,nn2, MBLK )
    implicit none
    real(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK
    integer :: n,ns,ne,m,ms,me,mm,nn,MBLKH,ML0
    real(8),allocatable :: utmp1(:), utmp2(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
!#ifndef _TIMER_
!    call write_border( 1, " d_gram_schmidt_t_sub(start)" )
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
             call d_gram_schmidt_t_sub( u, v, ms,me, ns,ne, MBLKH )

          end if

       end do ! ns

    end do ! ms
!#ifndef _TIMER_
!    call write_border( 1, " d_gram_schmidt_t_sub(end)" )
!#endif
    return

  END SUBROUTINE d_gram_schmidt_t_sub

!-----------------------------------------------------------------------

  SUBROUTINE z_gram_schmidt_t( u, v )
    implicit none
    complex(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: irank_b,np_band,myrank_b,ns,ne,ms,me
    integer :: nbss,ib,NBAND_BLK,ncycle,icycle,MB

    call write_border( 1, " z_gram_schmidt_t2(start)" )

    MB       = size( u, 2 )
    np_band  = v(2)%pinfo%np
    myrank_b = v(2)%pinfo%me

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK+1

    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ns = NBAND_BLK*(icycle-1) + 1
       ne = min( ns+NBAND_BLK-1, MB )

       if ( irank_b == myrank_b ) then

          call z_gram_schmidt_t_sub( u, v(1), ns,ne, ns,ne, NBLK )

       end if

       call rsdft_bcast( u(:,ns:ne), irank_b, v(2)%pinfo%comm )

       if ( ns <= MB-NBAND_BLK ) then

          do ib=1,(ncycle-1)/np_band+1

             nbss = (ib-1)*np_band + myrank_b + 1

             if ( nbss <= ncycle .and. nbss >= icycle+1 ) then

                ms = NBAND_BLK*(nbss-1)+1
                me = min(ms+NBAND_BLK-1,MB)

                if ( ms <= me ) then
                   call z_gram_schmidt_t_sub( u, v(1), ms,me, ns,ne, NBLK )
                end if

             end if

          end do

       end if

    end do ! icycle

    call write_border( 1, " z_gram_schmidt_t2(end)" )

    return

  END SUBROUTINE z_gram_schmidt_t

  RECURSIVE SUBROUTINE z_gram_schmidt_t_sub( u, v, mm1,mm2, nn1,nn2, MBLK )
    implicit none
    complex(8),intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK
    integer :: n,ns,ne,m,ms,me,mm,nn,MBLKH,ML0
    complex(8),allocatable :: utmp1(:), utmp2(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
!#ifndef _TIMER
!    call write_border( 1, " z_gram_schmidt_t_sub(start)" )
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

             call zgemm('N','N',ML0,mm,nn,-one,u(:,ns) &
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

                   call zgemv('N',ML0,n-ns+1,-one,u(:,ns) &
                              ,ML0,utmp1,1,one,u(:,m),1)

                end if

                if ( m==1 .or. ( n==m-1 .and. m/=ns ) ) then

                   call normalize_vector( u(:,m), v )

                end if

             end do ! m

             deallocate( utmp1 )

          else

             MBLKH=max(MBLK/2,NBLK1)
             call z_gram_schmidt_t_sub( u, v, ms,me, ns,ne, MBLKH )

          end if

       end do ! ns

    end do ! ms
!#ifndef _TIMER
!    call write_border( 1, " z_gram_schmidt_t_sub(end)" )
!#endif
    return

  END SUBROUTINE z_gram_schmidt_t_sub


END MODULE gram_schmidt_t2_module
