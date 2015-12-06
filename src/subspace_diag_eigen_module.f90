MODULE subspace_diag_eigen_module

#ifdef _EIGEN_
  use eigen_libs
#endif
  use watch_module
  use wf_module, only: esp, MB_WF
  use subspace_diag_variables, only: Hsub, Vsub, zero
  use subspace_mate_eigen_module
  use subspace_solv_eigen_module
  use subspace_rotv_eigen_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_eigen

CONTAINS


  SUBROUTINE subspace_diag_eigen( k,s, disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    integer,intent(IN) :: k,s
#ifdef _EIGEN_
    integer :: LDR,LDC,n
    real(8) :: time_temp(2),time_eigen(2,8)
    character(4) :: time_index(8)

    call write_border( 1, " subspace_diag_eigen(start)" )

    time_eigen(:,:)=0.0d0
    time_index(:) = (/"init","alc1","mate","alc2" &
                     ,"solv","rotv","deal","free"/)

    call watchb( time_temp )

    call eigen_init()

    call watchb( time_temp, time_eigen(1,1) )

    call eigen_get_matdims( MB_WF, LDR, LDC )

    allocate( Hsub(LDR,LDC) ) ; Hsub=zero

    call watchb( time_temp, time_eigen(1,2) )

    call subspace_mate_eigen(k,s)

    call watchb( time_temp, time_eigen(1,3) )

    allocate( Vsub(LDR,LDC) ) ; Vsub=zero

    call watchb( time_temp, time_eigen(1,4) )

#ifdef _DRSDFT_
    call subspace_solv_eigen( MB_WF, Hsub, esp(:,k,s), Vsub )
#endif

    call watchb( time_temp, time_eigen(1,5) )

    call subspace_rotv_eigen(k,s)

    call watchb( time_temp, time_eigen(1,6) )

    deallocate( Vsub )
    deallocate( Hsub )

    call watchb( time_temp, time_eigen(1,7) )

    call eigen_free()

    call watchb( time_temp, time_eigen(1,8) )

    if ( disp_switch ) call write_watchb( time_eigen, 8, time_index )

    call write_border( 1, " subspace_diag_eigen(end)" )
#endif
  END SUBROUTINE subspace_diag_eigen


END MODULE subspace_diag_eigen_module
