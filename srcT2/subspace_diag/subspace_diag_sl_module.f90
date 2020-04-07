MODULE subspace_diag_sl_module

  use scalapack_module, only: LLD_R,LLD_C,imate,prep_scalapack &
                             ,LLD_R_e,LLD_C_e,MBSIZE,NBSIZE,irotv,idiag &
                             ,d_gather_matrix_scalapack
  use subspace_mate_sl_module
  use subspace_mate_sl_bp0_module, only: subspace_mate_sl_bp0
  use subspace_solv_sl_module
  use subspace_rotv_sl_module
  use subspace_rotv_sl_bp0_module, only: subspace_rotv_sl_bp0
  use subspace_diag_variables, only: MB_diag,Hsub,zero,Vsub,Hsub_e,Vsub_e
  use watch_module
  use wf_module, only: unk
  use vector_tools_module, only: vinfo
!  use reconst_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_sl

CONTAINS

  SUBROUTINE subspace_diag_sl( k, s, v )
    implicit none
    integer,intent(IN) :: k,s
    type(vinfo),intent(IN) :: v(2)
    type(time) :: t
    integer :: ml0,mb0,mbt,mblk,np_band,ierr,myrank
    real(8),allocatable :: u(:,:)
    include 'mpif.h'
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    call write_border( 1, " subspace_diag_sl(start)" )
    call start_timer( t )

    np_band = v(2)%pinfo%np

    call prep_scalapack( MB_diag, v )

    allocate( Hsub(LLD_R,LLD_C) )
    Hsub=zero

    if ( imate >= 0 ) then

       select case( imate )
       case( 0 )
          call subspace_mate_sl_bp0( k, s, unk(:,:,k,s), v )
       case( 1 )
          call subspace_mate_sl_bp1( k, s, unk(:,:,k,s), v )
       case( 2 )
          call subspace_mate_sl_bp2( k, s, unk(:,:,k,s), v )
       case( 3 )
          call subspace_mate_sl_bp3( k, s, unk(:,:,k,s), v )
       end select

    else

       if ( np_band > 1 ) goto 900
       call subspace_mate_sl(k,s)

    end if

    allocate( Vsub(LLD_R,LLD_C) )
    Vsub=zero

    if ( idiag == "eigen_s" ) call alloc_dealloc_eigen( "A" )

    call subspace_solv_sl(k,s)

!    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if ( idiag == "eigen_s" ) call alloc_dealloc_eigen( "D" )

    deallocate( Hsub )

    if ( irotv >= 0 ) then

       !call reconstruct_init( v(2)%pinfo%comm )

       select case( irotv )
       case( 0 )
          call subspace_rotv_sl_bp0( unk(:,:,k,s), v )
       case( 1 )
          call subspace_rotv_sl_bp1( unk(:,:,k,s), v )
       case( 2 )
          call subspace_rotv_sl_bp2( unk(:,:,k,s), v )
       end select

    else

       if ( np_band > 1 ) goto 900
       call subspace_rotv_sl( unk(:,:,k,s), v )

    end if

    deallocate( Vsub )

    call result_timer( "sd", t )
    call write_border( 1, " subspace_diag_sl(end)" )

    return

900 call stop_program("This routine is incompatible with memory-band-parallel")

  END SUBROUTINE subspace_diag_sl


  SUBROUTINE alloc_dealloc_eigen( ctrl )
    implicit none
    character(1),intent(IN) :: ctrl
    if ( ctrl == "a" .or. ctrl == "A" ) then
       if ( .not.allocated(Hsub_e) ) then
          allocate( Hsub_e(LLD_R_e,LLD_C_e) )
          Hsub_e=zero
       end if
       if ( .not.allocated(Vsub_e) ) then
          allocate( Vsub_e(LLD_R_e,LLD_C_e) )
          Vsub_e=zero
       end if
    else if ( ctrl == "d" .or. ctrl == "D" ) then
       if ( allocated(Hsub_e) ) deallocate(Hsub_e)
       if ( allocated(Vsub_e) ) deallocate(Vsub_e)
    end if
  END SUBROUTINE alloc_dealloc_eigen


END MODULE subspace_diag_sl_module
