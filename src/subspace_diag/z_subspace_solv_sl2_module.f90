module z_subspace_solv_sl2_module

  implicit none

  private
  public :: z_subspace_solv_sl2

contains

  subroutine z_subspace_solv_sl2( sl, eig )
    use sl_variables, only: slinfo2, sl2, Hsub
    use sl_tools_module, only: slinfo, gather_matrix, distribute_matrix
    use rsdft_allreduce_module, only: rsdft_allreduce
    use rsdft_bcast_module, only: d_rsdft_bcast
    implicit none
    type(slinfo2),intent(inout) :: sl
    real(8),intent(inout) :: eig(:)
    integer :: ierr,nb,i
    integer :: LZWORK,LRWORK,LIWORK,NP0,NQ0,NPX,NQX
    integer :: itmp(1)
    integer,allocatable :: iwork(:)
    real(8) :: rtmp(1)
    real(8),allocatable :: rwork(:)
    complex(8) :: ztmp(1)
    complex(8),allocatable :: Vsub(:,:), zwork(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    character(100) :: msg
    character(1) :: ul
    logical :: disp_on
#ifndef _LAPACK_
    integer,external :: NUMROC
    type(slinfo) :: sl_

    call write_border( 1, " z_subspace_solv_sl2(start)" )
    call check_disp_switch( disp_on, 0 )

    ierr = 0
    nb = sl%nband
    ul = sl%uplo

    if ( sl%uplo /= '' ) then

      select case( sl%idiag )
      case( 'PZHEEVD' )

        allocate( Vsub(size(Hsub,1),size(Hsub,2)) ); Vsub=z0

        if ( sl%LWORK == 0 ) then
          call PZHEEVD('V',ul,nb,Hsub,1,1,sl%desca,eig,Vsub,1,1 &
                      ,sl%descz,ztmp,-1,rtmp,-1,itmp,-1,ierr)
          LZWORK = nint( real(ztmp(1)) )
          LRWORK = nint( rtmp(1) )
          LIWORK = itmp(1)
          NP0 = NUMROC(nb,sl%mbsize,0,0,sl%nprow)
          NQ0 = NUMROC(nb,sl%nbsize,0,0,sl%npcol)
          NPX = NUMROC(nb,sl%mbsize,sl%myrow,0,sl%nprow)
          NQX = NUMROC(nb,sl%nbsize,sl%mycol,0,sl%npcol)
          if ( disp_on ) then
            write(*,'("(PZHEEVD) Work-array sizes (by query)   : LWORK,LRWORK,LIWORK=",3i8)') LZWORK,LRWORK,LIWORK
            write(*,'("(PZHEEVD) Work-array sizes (in document): LWORK,LRWORK,LIWORK=",3i8)') &
              nb+(NP0+NQ0+sl%mbsize)*sl%mbsize, 1+9*nb+3*NPX*NQX, 7*nb+8*sl%npcol+2
          end if
          LZWORK = max( LZWORK, nb+(NP0+NQ0+sl%mbsize)*sl%mbsize )
          LRWORK = max( LRWORK, (1+9*nb+3*NPX*NQX) )
          LIWORK = max( LIWORK, 7*nb+8*sl%npcol+2 )
          if ( disp_on ) then
            write(*,'("(PZHEEVD) Work-array sizes (allocate)   : LWORK,LRWORK,LIWORK=",3i8)') LZWORK,LRWORK,LIWORK
          end if
          sl%lwork = LZWORK
          sl%lrwork = LRWORK
          sl%liwork = LIWORK
        else
          LZWORK = sl%lwork
          LRWORK = sl%lrwork
          LIWORK = sl%liwork
        end if

        allocate( zwork(LZWORK) ); zwork=z0
        allocate( rwork(LRWORK) ); rwork=z0
        allocate( iwork(LIWORK) ); iwork=0

        call PZHEEVD('V',ul,nb,Hsub,1,1,sl%desca,eig,Vsub,1,1 &
                    ,sl%descz,zwork,LZWORK,rwork,LRWORK,iwork,LIWORK,ierr)

        deallocate( iwork )
        deallocate( rwork )
        deallocate( zwork )

        Hsub = Vsub

        deallocate( Vsub )

      end select

    end if

    if ( ierr /= 0 ) then
      write(msg,'("ierr=",i4,"(stop@z_subspace_solv_sl2)")') ierr
      call stop_program( msg )
    end if

    ! ---
    !
    call d_rsdft_bcast( eig, comm_chr='grid' )
    call d_rsdft_bcast( eig, comm_chr='band' )
    !
    ! ---
    !
    if ( allocated(Hsub) ) then

      sl_%nprow  = sl%nprow
      sl_%npcol  = sl%npcol
      sl_%myrow  = sl%myrow
      sl_%mycol  = sl%mycol
      sl_%mbsize = sl%mbsize
      sl_%nbsize = sl%nbsize

      allocate( Vsub(nb,nb) ); Vsub=z0
      call gather_matrix( sl_, Hsub, Vsub )

      deallocate( Hsub )
    
    end if
    !
    ! ---

    if ( .not.allocated(Vsub) ) then
      allocate( Vsub(nb,nb) ); Vsub=z0
    end if

    call rsdft_allreduce( Vsub, 'g' )
    call rsdft_allreduce( Vsub, 'b' )

    ! ---
    !
    allocate( Hsub(sl2%ldr,sl2%ldc) ); Hsub=z0

    sl_%nprow  = sl2%nprow
    sl_%npcol  = sl2%npcol
    sl_%myrow  = sl2%myrow
    sl_%mycol  = sl2%mycol
    sl_%mbsize = sl2%mbsize
    sl_%nbsize = sl2%nbsize

    call distribute_matrix( sl_, Vsub, Hsub )
    !
    ! ---

    call write_border( 1, " z_subspace_solv_sl2(end)" )

    return

#endif

  end subroutine z_subspace_solv_sl2

end module z_subspace_solv_sl2_module
