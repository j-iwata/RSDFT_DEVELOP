MODULE velocity_scaling_ele_module

  use cpmd_variables, only: ekin1, ekinw

  implicit none

  PRIVATE
  PUBLIC :: velocity_scaling_ele

CONTAINS

  SUBROUTINE velocity_scaling_ele( fke, psi_v )
    implicit none
    real(8),intent(IN) :: fke
    real(8),intent(INOUT) :: psi_v(:,:,:,:)
    real(8) :: c
    if ( fke > ekin1 .and. fke /= 0.0d0 ) then
       c = sqrt( ekinw/fke )
       psi_v = c*psi_v
    end if
    return
  END SUBROUTINE velocity_scaling_ele

END MODULE velocity_scaling_ele_module
