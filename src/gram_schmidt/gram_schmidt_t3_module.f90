module gram_schmidt_t3_module

  use rgrid_variables, only: dV

  implicit none

  private
  public :: gram_schmidt_t3

  interface gram_schmidt_t3
    module procedure d_gram_schmidt_t3, z_gram_schmidt_t3
  end interface

  integer :: NBLK =0
  integer :: NBLK1=0
  integer :: comm_grid, comm_band
  real(8),allocatable :: dtmp1(:), dtmp2(:,:)
  complex(8),allocatable :: ztmp1(:), ztmp2(:,:)

contains

  subroutine d_gram_schmidt_t3( u, comm_grid_in, comm_band_in, u2 )
    use rsdft_bcast_module, only: d_rsdft_bcast_tmp
    implicit none
    real(8), intent(inout) :: u(:,:)
    integer, intent(in) :: comm_grid_in, comm_band_in
    real(8), optional, intent(inout) :: u2(:,:)
    integer :: nbss,k1,ib,NBAND_BLK,ncycle
    integer :: irank_b,myrank_b,np_band
    integer :: ML0,MB,ns,ne,ms,me,n

    call write_border( 1, " d_gram_schmidt_t3(start)" )

    comm_grid = comm_grid_in
    comm_band = comm_band_in
    call prep_mpi( myrank_b, np_band, comm_band )

    ML0 = size( u, 1 )
    MB  = size( u, 2 )

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK + 1

    do k1=1,ncycle

      ns = NBAND_BLK*(k1-1) + 1
      ne = min( ns+NBAND_BLK-1, MB )

      irank_b = mod( k1-1, np_band )

      if ( myrank_b == irank_b ) then

        if ( present(u2) ) then
          call d_gram_schmidt_sub(u,ns,ne,ns,ne,NBLK,u2)
        else
          call d_gram_schmidt_sub(u,ns,ne,ns,ne,NBLK)
        end if

      end if

      if ( np_band > 1 ) then
        n=ML0*(ne-ns+1)
        call d_rsdft_bcast_tmp( u(:,ns:ne), n, irank_b, comm_chr='band' )
        if ( present(u2) ) call d_rsdft_bcast_tmp( u2(:,ns:ne), n, irank_b, comm_chr='band' )
      end if

      if ( ns <= MB-NBAND_BLK ) then

        do ib=1,(ncycle-1)/np_band+1

          nbss=(ib-1)*np_band+myrank_b+1

          if ( nbss <= ncycle .and. nbss >= k1+1 ) then

            ms=NBAND_BLK*(nbss-1)+1
            me=min(ms+NBAND_BLK-1,MB)

            if ( ms <= me ) then
              if ( present(u2) ) then
                call d_gram_schmidt_sub(u,ms,me,ns,ne,NBLK,u2)
              else
                call d_gram_schmidt_sub(u,ms,me,ns,ne,NBLK)
              end if
            end if

            end if

        end do

      end if

    end do ! k1

    call write_border( 1, " d_gram_schmidt_t3(end)" )

  end subroutine d_gram_schmidt_t3

  recursive subroutine d_gram_schmidt_sub(u,mm1,mm2,nn1,nn2,MBLK,u2)
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    real(8),intent(inout) :: u(:,:)
    integer,intent(in) :: mm1,mm2,nn1,nn2,MBLK
    real(8),optional,intent(inout) :: u2(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8) :: a
    real(8) :: c
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr,ML0,i

    ML0 = size( u, 1 )
    a = -dV

    do ms=mm1,mm2,MBLK

      me=min(ms+MBLK-1,mm2)
      mm=me-ms+1

      do ns=nn1,nn2,MBLK

        ne=min(ns+MBLK-1,nn2)
        ne=min(ne,me-1)
        nn=ne-ns+1

        if ( nn <= 0 ) cycle

        if ( ms >= ne+1 ) then

          allocate( dtmp2(ns:ne,ms:me) )

          call DGEMM('T','N',nn,mm,ML0,a,u(1,ns),ML0,u(1,ms),ML0,zero,dtmp2,nn)
          call rsdft_allreduce( dtmp2, 'grid' )

          call DGEMM('N','N',ML0,mm,nn,one,u(1,ns),ML0,dtmp2,nn,one,u(1,ms),ML0)
          if ( present(u2) ) then
            call DGEMM('N','N',ML0,mm,nn,one,u2(1,ns),ML0,dtmp2,nn,one,u2(1,ms),ML0)
          end if

          deallocate( dtmp2 )

          if ( ms == ne+1 ) then

            c=sum( u(:,ms)*u(:,ms) )
            call rsdft_allreduce( c, 'grid' )

            c=1.0d0/sqrt(c*dV)
            u(:,ms) = c*u(:,ms)
            if ( present(u2) ) u2(:,ms)=c*u2(:,ms)

          end if

        else if ( mm <= NBLK1 ) then

          allocate( dtmp1(NBLK1) )

          do m=ms,me

            n=min(m-1,ne)

            if ( n-ns+1 > 0 ) then

              call DGEMV('T',ML0,n-ns+1,a,u(1,ns),ML0,u(1,m),1,zero,dtmp1,1)
              call rsdft_allreduce( dtmp1(1:n-ns+1), 'grid' )

              call DGEMV('N',ML0,n-ns+1,one,u(1,ns),ML0,dtmp1,1,one,u(1,m),1)
              if ( present(u2) ) then
                call DGEMV('N',ML0,n-ns+1,one,u2(1,ns),ML0,dtmp1,1,one,u2(1,m),1)
              end if

            end if

            if ( m==1 .or. (n==m-1 .and. m/=ns) ) then

              c=sum( u(:,m)*u(:,m) )
              call rsdft_allreduce( c, 'grid' )

              c=1.0d0/sqrt(c*dV)

              u(:,m) = c*u(:,m)
              if ( present(u2) ) u2(:,m) = c*u2(:,m)

            end if

          end do ! m

          deallocate( dtmp1 )

        else

          MBLKH=max(MBLK/2,NBLK1)
          if ( present(u2) ) then
            call d_gram_schmidt_sub(u,ms,me,ns,ne,MBLKH,u2)
          else
            call d_gram_schmidt_sub(u,ms,me,ns,ne,MBLKH)
          end if

        end if

      end do ! ns

    end do ! ms

  end subroutine d_gram_schmidt_sub


  subroutine prep_mpi( myrank, nprocs, comm )
    implicit none
    integer, intent(out) :: myrank, nprocs
    integer, intent(in)  :: comm
    integer :: ierr
    call MPI_Comm_rank( comm, myrank, ierr ) 
    call MPI_Comm_size( comm, nprocs, ierr ) 
  end subroutine prep_mpi


!-----------------------------------------------------------------------

  subroutine z_gram_schmidt_t3( u, comm_grid_in, comm_band_in, u2 )
    use rsdft_bcast_module, only: z_rsdft_bcast_tmp
    implicit none
    complex(8), intent(inout) :: u(:,:)
    integer, intent(in) :: comm_grid_in, comm_band_in
    complex(8), optional, intent(inout) :: u2(:,:)
    integer :: irank_b,myrank_b,np_band
    integer :: nbss,k1,ib,NBAND_BLK,ncycle
    integer :: ML0,MB,ns,ne,ms,me,n

    call write_border( 1, " z_gram_schmidt_t3(start)" )

    comm_grid = comm_grid_in
    comm_band = comm_band_in
    call prep_mpi( myrank_b, np_band, comm_band )

    ML0 = size( u, 1 )
    MB  = size( u, 2 )

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    NBAND_BLK = NBLK
    ncycle    = (MB-1)/NBAND_BLK + 1

    do k1=1,ncycle

      ns = NBAND_BLK*(k1-1) + 1
      ne = min( ns+NBAND_BLK-1, MB )

      irank_b = mod( k1-1, np_band )

      if ( myrank_b == irank_b ) then

        if ( present(u2) ) then
          call z_gram_schmidt_sub(u,ns,ne,ns,ne,NBLK,u2)
        else
          call z_gram_schmidt_sub(u,ns,ne,ns,ne,NBLK)
        end if

      end if

      if ( np_band > 1 ) then
        n=ML0*(ne-ns+1)
        call z_rsdft_bcast_tmp( u(:,ns:ne), n, irank_b, comm_chr='band' )
        if ( present(u2) ) call z_rsdft_bcast_tmp( u2(:,ns:ne), n, irank_b, comm_chr='band' )
      end if

      if ( ns <= MB-NBAND_BLK ) then

        do ib=1,(ncycle-1)/np_band+1

          nbss=(ib-1)*np_band+myrank_b+1

          if ( nbss <= ncycle .and. nbss >= k1+1 ) then

            ms=NBAND_BLK*(nbss-1)+1
            me=min(ms+NBAND_BLK-1,MB)

            if ( ms <= me ) then
              if ( present(u2) ) then
                call z_gram_schmidt_sub(u,ms,me,ns,ne,NBLK,u2)
              else
                call z_gram_schmidt_sub(u,ms,me,ns,ne,NBLK)
              end if
            end if

          end if

        end do

      end if

    end do ! k1

    call write_border( 1, " z_gram_schmidt_t3(end)" )

  end subroutine z_gram_schmidt_t3

  recursive subroutine z_gram_schmidt_sub(u,mm1,mm2,nn1,nn2,MBLK,u2)
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    complex(8),intent(inout) :: u(:,:)
    integer,intent(in) :: mm1,mm2,nn1,nn2,MBLK
    complex(8),optional,intent(inout) :: u2(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    complex(8) :: a
    real(8) :: c
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr,ML0,i

    ML0 = size( u, 1 )
    a = -dV

    do ms=mm1,mm2,MBLK

      me=min(ms+MBLK-1,mm2)
      mm=me-ms+1

      do ns=nn1,nn2,MBLK

        ne=min(ns+MBLK-1,nn2)
        ne=min(ne,me-1)
        nn=ne-ns+1

        if ( nn <= 0 ) cycle

        if ( ms >= ne+1 ) then

          allocate( ztmp2(ns:ne,ms:me) )

          call ZGEMM('C','N',nn,mm,ML0,a,u(1,ns),ML0,u(1,ms),ML0,zero,ztmp2,nn)
          call rsdft_allreduce( ztmp2, 'grid' )

          call ZGEMM('N','N',ML0,mm,nn,one,u(1,ns),ML0,ztmp2,nn,one,u(1,ms),ML0)
          if ( present(u2) ) then
            call ZGEMM('N','N',ML0,mm,nn,one,u2(1,ns),ML0,ztmp2,nn,one,u2(1,ms),ML0)
          end if

          deallocate( ztmp2 )

          if ( ms == ne+1 ) then

            c=sum( abs(u(:,ms))**2 )
            call rsdft_allreduce( c, 'grid' )

            c=1.0d0/sqrt(c*dV)
            u(:,ms) = c*u(:,ms)
            if ( present(u2) ) u2(:,ms)=c*u2(:,ms)

          end if

        else if ( mm <= NBLK1 ) then

          allocate( ztmp1(NBLK1) )

          do m=ms,me

            n=min(m-1,ne)

            if ( n-ns+1 > 0 ) then

              call ZGEMV('C',ML0,n-ns+1,a,u(1,ns),ML0,u(1,m),1,zero,ztmp1,1)
              call rsdft_allreduce( ztmp1(1:n-ns+1), 'grid' )

              call ZGEMV('N',ML0,n-ns+1,one,u(1,ns),ML0,ztmp1,1,one,u(1,m),1)
              if ( present(u2) ) then
                call ZGEMV('N',ML0,n-ns+1,one,u2(1,ns),ML0,ztmp1,1,one,u2(1,m),1)
              end if

            end if

            if ( m==1 .or. (n==m-1 .and. m/=ns) ) then

              c=sum( abs(u(:,m))**2 )
              call rsdft_allreduce( c, 'grid' )

              c=1.0d0/sqrt(c*dV)

              u(:,m) = c*u(:,m)
              if ( present(u2) ) u2(:,m) = c*u2(:,m)

            end if

          end do ! m

          deallocate( ztmp1 )

        else

          MBLKH=max(MBLK/2,NBLK1)
          if ( present(u2) ) then
            call z_gram_schmidt_sub(u,ms,me,ns,ne,MBLKH,u2)
          else
            call z_gram_schmidt_sub(u,ms,me,ns,ne,MBLKH)
          end if

        end if

      end do ! ns

    end do ! ms

  end subroutine z_gram_schmidt_sub

end module gram_schmidt_t3_module
