MODULE ps_init_module

  use ps_init_sol_module
  use ps_init_mol_module

  implicit none

  PRIVATE
  PUBLIC :: ps_init

CONTAINS

  SUBROUTINE ps_init( SYStype, Gcut, rho )
    implicit none
    integer,intent(IN) :: SYStype
    real(8),intent(IN) :: Gcut
    real(8),intent(INOUT) :: rho(:,:)

    select case( SYStype )
    case( 0 )
       call ps_init_sol( Gcut, rho )
    case( 1 )
       call ps_init_mol( Gcut, rho )
    end select
  END SUBROUTINE ps_init

END MODULE ps_init_module
