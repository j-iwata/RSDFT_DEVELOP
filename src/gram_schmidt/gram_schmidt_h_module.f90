MODULE gram_schmidt_h_module

  use io_tools_module, only: IOTools_readIntegerKeyword
  use vector_tools_module, only: normalize_vector,inner_product_vector,vinfo
  use rsdft_mpi_module, only: rsdft_bcast

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt_h
  PUBLIC :: read_gram_schmidt_h

  integer :: NBLK=0
  integer :: NBLK1=0

  INTERFACE gram_schmidt_h
     MODULE PROCEDURE d_gram_schmidt_h
  END INTERFACE

CONTAINS


  SUBROUTINE read_gram_schmidt_h
    implicit none
    call IOTools_readIntegerKeyword( "NBLK", NBLK )
    call IOTools_readIntegerKeyword( "NBLK1", NBLK1 )
  END SUBROUTINE read_gram_schmidt_h


  SUBROUTINE d_gram_schmidt_h( u, v )

    implicit none
    real(8),intent(INOUT)  :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    integer :: NBAND_BLK, np_band, myrank_b, comm_b
    integer :: icycle, irank_b, ib, MB, ncycle, ncycle2, nxx
    integer :: nn, ns, ne, ns_next, ne_next, ms, me, nbss
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE), ierr
    integer :: ireq, rep_rank, tru_rank, tri_cal, irank_t, iwait

    call write_border( 1, " d_gram_schmidt_h(start)" )

    MB = size( u, 2 )

    np_band  = v(2)%pinfo%np
    myrank_b = v(2)%pinfo%me
    comm_b   = v(2)%pinfo%comm

    if ( NBLK  == 0 ) NBLK  = MB
    if ( NBLK1 == 0 ) NBLK1 = 4

    NBAND_BLK = NBLK

    ncycle  = int( (MB-1)/NBAND_BLK ) + 1
    ncycle2 = int( (ncycle-1)/np_band ) + 1

    nxx = mod( ncycle, np_band )
    rep_rank = -1
    tru_rank = -1
    tri_cal  =  0
    iwait    =  0

    do icycle=1,ncycle

       irank_b = mod( icycle-1, np_band )

       ns = NBAND_BLK*(icycle-1)+1
       ne = min(ns+NBAND_BLK-1,MB)

       ns_next = NBAND_BLK*icycle+1
       ne_next = min( ns_next+NBAND_BLK-1, MB )

       if ( icycle==1 .or. ncycle-icycle+1<=np_band .or. irank_b==nxx ) then
          irank_t = irank_b
          if ( irank_b == myrank_b ) tri_cal=1
       else
          irank_t = rep_rank
       end if

       if ( iwait == 1 ) then
          call mpi_wait( ireq, istatus, ierr )
          iwait = 0
       end if

       if ( tri_cal == 1 ) then
#ifdef _TIMER_
!         call MPI_Barrier(v(1)%pinfo%comm, ierr)
          call timer_sta(2)
#endif
          call d_gram_schmidt_h_sub( u,v(1), ns,ne, ns,ne, NBLK )
#ifdef _TIMER_
!         call MPI_Barrier(v(1)%pinfo%comm, ierr)
          call timer_end(2)
#endif
          tri_cal=0
       end if

       call mpi_barrier( comm_b, ierr )

       call rsdft_bcast( u(:,ns:ne), irank_t, comm_b )

       rep_rank=irank_b
       if ( mod(icycle,np_band)==nxx .or. ncycle-icycle<=np_band ) rep_rank=-1

       if ( rep_rank>=0 .and. rep_rank==myrank_b ) then
          tru_rank=rep_rank+1
          if ( tru_rank == np_band ) tru_rank=0
          nn=size(u,1)*NBAND_BLK
          !nn=size(u,1)*(ne_next-ns_next+1)
          call mpi_irecv( u(:,ns_next:ne_next),nn,MPI_REAL8, &
                          tru_rank,1000,comm_b,ireq,ierr )
          tri_cal=1
          iwait=1
       end if

       if ( ns <= MB-NBAND_BLK ) then

          do ib=1,ncycle2

             nbss=(ib-1)*np_band+myrank_b+1

             if ( nbss<=ncycle .and. nbss>= icycle+1 ) then

                ms=NBAND_BLK*(nbss-1)+1
                me=min(ms+NBAND_BLK-1,MB)

                if ( ms <= me ) then
#ifdef _TIMER_
                   call timer_sta(3)
#endif
                   call d_gram_schmidt_h_sub( u,v(1), ms,me, ns,ne, NBLK )
#ifdef _TIMER_
                   call timer_end(3)
#endif
                end if

                if ( nbss == icycle+1 .and. rep_rank >= 0 ) then
                   nn=size(u,1)*(me-ms+1)
                   call mpi_isend( u(:,ms:me),nn,MPI_REAL8, &
                                   rep_rank,1000,comm_b,ireq,ierr )
                   iwait = 1
                end if

             end if

          end do ! ib

       end if

    end do ! icycle

    call write_border( 1, " d_gram_schmidt_h(end)" )

    return
  END SUBROUTINE d_gram_schmidt_h


  RECURSIVE SUBROUTINE d_gram_schmidt_h_sub( u,v, mm1,mm2, nn1,nn2, MBLK )

    implicit none
    real(8),intent(INOUT)  :: u(:,:)
    type(vinfo),intent(IN) :: v
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr
    integer :: ML0,i
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: utmp1(:), utmp2(:,:)

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

          call dgemm('N','N',ML0,mm,nn,-one,u(:,ns:ne),ML0 &
               ,utmp2,nn,one,u(:,ms:me),ML0)

          deallocate( utmp2 )

          if ( ms == ne+1 ) call normalize_vector( u(:,ms), v )

       else if ( mm <= NBLK1 ) then

          allocate( utmp1(NBLK1) ) ; utmp1=zero

          do m=ms,me

             n=min(m-1,ne)

             if ( n-ns+1 > 0 ) then

                call inner_product_vector( u(:,ns:n), u(:,m), utmp1, v )

                call dgemv('N',ML0,n-ns+1,-one,u(:,ns),ML0,utmp1,1,one,u(:,m),1)

             end if

             if ( m == 1 .or. ( n==m-1 .and. m /= ns ) ) then
                call normalize_vector( u(:,m), v )
             end if

          end do ! m

          deallocate( utmp1 )

       else

          MBLKH = max( MBLK/2, NBLK1 )
          call d_gram_schmidt_h_sub( u,v, ms,me, ns,ne, MBLKH )

       end if

    end do ! ns
    end do ! ms

    return
  END SUBROUTINE d_gram_schmidt_h_sub


END MODULE gram_schmidt_h_module
