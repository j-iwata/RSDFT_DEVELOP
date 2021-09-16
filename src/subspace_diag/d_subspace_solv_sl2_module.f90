module d_subspace_solv_sl2_module

  implicit none

  private
  public :: d_subspace_solv_sl2

contains

  subroutine d_subspace_solv_sl2( sl, eig )
    use sl_variables, only: slinfo2, sl2, Dsub
    use sl_tools_module, only: slinfo, gather_matrix, distribute_matrix
    use rsdft_allreduce_module, only: rsdft_allreduce
    use rsdft_bcast_module, only: d_rsdft_bcast
    implicit none
    type(slinfo2),intent(inout) :: sl
    real(8),intent(inout) :: eig(:)
    integer :: ierr,nb,i
    integer :: LRWORK,LIWORK,TRILWMIN,NP,NQ
    integer :: itmp(1)
    integer,allocatable :: iwork(:)
    real(8) :: rtmp(1)
    real(8),allocatable :: rwork(:)
    real(8),allocatable :: Vsub(:,:)
    real(8),parameter :: z0=0.0d0
    character(100) :: msg
    character(1) :: ul
    logical :: disp_on
#ifndef _LAPACK_
    integer,external :: NUMROC
    type(slinfo) :: sl_

    call write_border( 1, " d_subspace_solv_sl2(start)" )
    call check_disp_switch( disp_on, 0 )

    ierr = 0
    nb = sl%nband
    ul = sl%uplo

    select case( sl%idiag )
    case( 'PDSYEVD' )

      allocate( Vsub(size(Dsub,1),size(Dsub,2)) ); Vsub=z0

      if ( sl%LWORK == 0 ) then
        call PDSYEVD('V',ul,nb,Dsub,1,1,sl%desca,eig,Vsub,1,1,sl%descz,rtmp,-1,itmp,-1,ierr)
        LRWORK = nint( rtmp(1) )
        LIWORK = itmp(1)
        NP = NUMROC(nb,sl%mbsize,sl%myrow,0,sl%nprow)
        NQ = NUMROC(nb,sl%nbsize,sl%mycol,0,sl%npcol)
        TRILWMIN = 3*nb + max( sl%mbsize*(NP+1), 3*sl%mbsize )
        if ( disp_on ) then
          write(*,'("(PDSYEVD) Work-array sizes (by query)   : LRWORK,LIWORK=",3i8)') LRWORK,LIWORK
          write(*,'("(PDSYEVD) Work-array sizes (in document): LRWORK,LIWORK=",3i8)') &
          max( 1+6*nb+2*NP*NQ, TRILWMIN ) + 2*nb, 7*nb + 8*sl%npcol + 2
        end if
        LRWORK = max( LRWORK, max(1+6*nb+2*NP*NQ,TRILWMIN)+2*nb )
        LIWORK = max( LIWORK, 7*nb+8*sl%npcol+2 )
        if ( disp_on ) then
          write(*,'("(PDSYEVD) Work-array sizes (allocate)   : LRWORK,LIWORK=",3i8)') LRWORK,LIWORK
        end if
        sl%lwork = LRWORK
        sl%lrwork = LRWORK
        sl%liwork = LIWORK
      else
        LRWORK = sl%lrwork
        LIWORK = sl%liwork
      end if

      allocate( rwork(LRWORK) ); rwork=z0
      allocate( iwork(LIWORK) ); iwork=0

      call PDSYEVD('V',ul,nb,Dsub,1,1,sl%desca,eig,Vsub,1,1,sl%descz,rwork,LRWORK,iwork,LIWORK,ierr)

      deallocate( iwork )
      deallocate( rwork )

      Dsub = Vsub

      deallocate( Vsub )

    end select

    if ( ierr /= 0 ) then
      write(msg,'("ierr=",i4,"(stop@d_subspace_solv_sl2)")') ierr
      call stop_program( msg )
    end if

    ! ---
    !
    call d_rsdft_bcast( eig, size(eig), 0, comm_label='grid' )
    call d_rsdft_bcast( eig, size(eig), 0, comm_label='band' )
    !
    ! ---
    !
    if ( allocated(Dsub) ) then

      sl_%nprow  = sl%nprow
      sl_%npcol  = sl%npcol
      sl_%myrow  = sl%myrow
      sl_%mycol  = sl%mycol
      sl_%mbsize = sl%mbsize
      sl_%nbsize = sl%nbsize

      allocate( Vsub(nb,nb) ); Vsub=z0
      call gather_matrix( sl_, Dsub, Vsub )

      deallocate( Dsub )
    
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
    allocate( Dsub(sl2%ldr,sl2%ldc) ); Dsub=z0

    sl_%nprow  = sl2%nprow
    sl_%npcol  = sl2%npcol
    sl_%myrow  = sl2%myrow
    sl_%mycol  = sl2%mycol
    sl_%mbsize = sl2%mbsize
    sl_%nbsize = sl2%nbsize

    call distribute_matrix( sl_, Vsub, Dsub )
    !
    ! ---

    call write_border( 1, " d_subspace_solv_sl2(end)" )

    return

#endif

  end subroutine d_subspace_solv_sl2

end module d_subspace_solv_sl2_module
