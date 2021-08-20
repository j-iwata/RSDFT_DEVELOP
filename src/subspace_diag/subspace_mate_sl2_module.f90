module subspace_mate_sl2_module

  use d_subspace_mate_sl2_module, only: d_subspace_mate_sl2
  use z_subspace_mate_sl2_module, only: z_subspace_mate_sl2

  implicit none

  private
  public :: subspace_mate_sl2

  interface subspace_mate_sl2
    module procedure d_subspace_mate_sl2, z_subspace_mate_sl2
  end interface

end module subspace_mate_sl2_module