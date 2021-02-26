module hartree_sol_module

  use hartree_sol_1dffte_module, only: calc_hartree_sol_1dffte
  use hartree_sol_ffte_module, only: calc_hartree_sol_ffte
  use hartree_sol_fftw_module, only: calc_hartree_sol_fftw
  use hartree_sol_fft0_module, only: calc_hartree_sol_fft0

  implicit none

  private
  public :: calc_hartree_sol

contains

  subroutine calc_hartree_sol( rho, iswitch_fft )
    implicit none
    real(8),intent(in) :: rho(:,:)
    character(*),intent(in) :: iswitch_fft

    call write_border( 1, " calc_hartree_sol(start)" )

    select case( iswitch_fft )
    case( 'FFTE1' )

      call calc_hartree_sol_ffte( rho )

    case( 'FFTE2' )

      call calc_hartree_sol_1dffte( rho )

    case( 'FFTW1' )

      !call calc_hartree_sol_fftw( rho )

    case default

      call calc_hartree_sol_fft0( rho )

    end select

    call write_border( 1, " calc_hartree_sol(end)" )

  end subroutine calc_hartree_sol

end module hartree_sol_module
