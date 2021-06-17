module hartree_sol_module

  use fft_module, only: iswitch_fft
  use hartree_sol_1dffte_module, only: init_hartree_sol_1dffte, calc_hartree_sol_1dffte
  use hartree_sol_ffte_module  , only: init_hartree_sol_ffte  , calc_hartree_sol_ffte
  use hartree_sol_fftw_module  , only: init_hartree_sol_fftw  , calc_hartree_sol_fftw
  use hartree_sol_fft0_module  , only: calc_hartree_sol_fft0

  implicit none

  private
  public :: init_hartree_sol
  public :: calc_hartree_sol

contains

  subroutine init_hartree_sol( Ngrid, Igrid )
    implicit none
    integer,intent(in) :: Ngrid(3), Igrid(2,3)
    call write_border( 0, " init_hartree_sol(start)" )

    select case( iswitch_fft )                                                                                                          
    case( 'FFTE', 'FFTE1' )                                                                                                               

      call init_hartree_sol_ffte( Ngrid(1:3), Igrid(:,1:3) )                                                                            

    case( 'FFTE2' )                                                                                                                     

      call init_hartree_sol_1dffte( Ngrid(1:3), Igrid(:,1:3) )                                                                          

    case( 'FFTW', 'FFTW1' )

      call init_hartree_sol_fftw( Ngrid(1:3) )

    end select

    call write_border( 0, " init_hartree_sol(end)" )
  end subroutine init_hartree_sol

  subroutine calc_hartree_sol( rho )
    implicit none
    real(8),intent(in) :: rho(:,:)

    call write_border( 1, " calc_hartree_sol(start)" )

    select case( iswitch_fft )
    case( 'FFTE', 'FFTE1' )

      call calc_hartree_sol_ffte( rho )

    case( 'FFTE2' )

      call calc_hartree_sol_1dffte( rho )

    case( 'FFTW', 'FFTW1' )

      call calc_hartree_sol_fftw( rho )

    case default

      call calc_hartree_sol_fft0( rho )

    end select

    call write_border( 1, " calc_hartree_sol(end)" )

  end subroutine calc_hartree_sol

end module hartree_sol_module
