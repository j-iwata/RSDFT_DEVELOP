!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!-------Berendsen thermostat for ions 08/16 KK,JI------------------------------

MODULE berendsen_module

  use calc_kine_temp_module

  implicit none

  PRIVATE
  PUBLIC :: berendsen

CONTAINS

  SUBROUTINE berendsen( temp, dt, Velocity )
    implicit none
    real(8),intent(INOUT) :: Velocity(:,:)
    real(8),intent(IN) :: temp, dt
    real(8),parameter :: thresh=1.d-6 ! don't scale if we are too "cold"
    real(8),parameter :: taubp =2.0d0 ! time constant
    real(8) :: temp_now,c,taubp2

    call calc_temp( Velocity, temp_now )

    if ( (temp_now/temp) > thresh ) then

       !taubp2 = taubp
       taubp2 = dt/0.0025d0

       c = sqrt( 1.0d0 + dt/taubp2*(temp/temp_now-1.0d0) )

       Velocity(:,:) = Velocity(:,:)*c

    end if

  END SUBROUTINE berendsen

END MODULE berendsen_module
