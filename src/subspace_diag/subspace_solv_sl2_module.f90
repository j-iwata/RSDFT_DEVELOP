module subspace_solv_sl2_module

  use d_subspace_solv_sl2_module, only: d_subspace_solv_sl2
  use z_subspace_solv_sl2_module, only: z_subspace_solv_sl2

  implicit none

  private
  public :: subspace_solv_sl2

contains

  subroutine subspace_solv_sl2( sl, esp )
    use sl_variables, only: slinfo2, Dsub
    implicit none
    type(slinfo2),intent(inout) :: sl
    real(8),intent(inout) :: esp(:)
    if ( allocated(Dsub) ) then
      call d_subspace_solv_sl2( sl, esp )
    else
      call z_subspace_solv_sl2( sl, esp )
    end if
  end subroutine subspace_solv_sl2
  
end module subspace_solv_sl2_module