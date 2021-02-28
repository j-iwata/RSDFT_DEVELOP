module construct_vion_vh_floc_module

  use fft_module, only: iswitch_fft
  use construct_vion_vh_floc_ffte_module, only: construct_vion_vh_floc_ffte
  use construct_vion_vh_floc_1dffte_module, only: construct_vion_vh_floc_1dffte

  implicit none

  private
  public :: construct_vion_vh_floc

contains

  subroutine construct_vion_vh_floc( rho, Vion, Vh, force_local, Eh )
    implicit none
    real(8),intent(in)  :: rho(:,:)
    real(8),intent(out) :: Vion(:), Vh(:), force_local(:,:)
    real(8),intent(out) :: Eh
    select case( iswitch_fft )
    case( 'FFTE', 'FFTE1' )
      call construct_vion_vh_floc_ffte( rho, Vion, Vh, force_local, Eh )
    case( 'FFTE2' )
      call construct_vion_vh_floc_1dffte( rho, Vion, Vh, force_local, Eh )
    end select
  end subroutine construct_vion_vh_floc

end module construct_vion_vh_floc_module
