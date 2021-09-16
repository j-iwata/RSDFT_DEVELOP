MODULE subspace_diag_sl_module

  use scalapack_module
  use subspace_mate_sl_module
  use subspace_solv_sl_module
  use subspace_rotv_sl_module
  use subspace_diag_variables, only: MB_diag,Hsub,zero,Vsub
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_sl

CONTAINS

  SUBROUTINE subspace_diag_sl( k, s )
    implicit none
    integer,intent(IN) :: k,s
    type(time) :: t

    call write_border( 1, " subspace_diag_sl(start)" )
    call start_timer( t_out=t )

    call prep_scalapack( MB_diag )

    allocate( Hsub(LLD_R,LLD_C) )
    Hsub=zero

    call subspace_mate_sl(k,s)

    allocate( Vsub(LLD_R,LLD_C) )
    Vsub=zero

    call subspace_solv_sl(k,s)

    deallocate( Hsub )

    call subspace_rotv_sl(k,s)

    deallocate( Vsub )

    call result_timer( "sd", t )
    call write_border( 1, " subspace_diag_sl(end)" )

  END SUBROUTINE subspace_diag_sl

END MODULE subspace_diag_sl_module
