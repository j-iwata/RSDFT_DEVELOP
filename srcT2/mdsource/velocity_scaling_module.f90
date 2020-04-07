MODULE velocity_scaling_module

  use cpmd_variables, only: trange
  use calc_kine_temp_module

  implicit none

  PRIVATE
  PUBLIC :: velocity_scaling

CONTAINS


  SUBROUTINE velocity_scaling( temp, Velocity )

    implicit none
    real(8),intent(IN) :: temp
    real(8),intent(INOUT) :: Velocity(:,:)
    real(8) :: temp_now, scale

    call calc_temp( Velocity, temp_now )

    if ( temp_now == 0.0d0 ) return

    if ( ( temp_now > temp+trange ) .or. ( temp_now < temp-trange ) ) then

       scale = sqrt( temp/temp_now )
       Velocity(:,:)=Velocity(:,:)*scale

    end if

  END SUBROUTINE velocity_scaling


END MODULE velocity_scaling_module
