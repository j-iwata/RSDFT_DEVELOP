module subspace_rotv_sl2_module

  use d_subspace_rotv_sl2_module, only: d_subspace_rotv_sl2
  use z_subspace_rotv_sl2_module, only: z_subspace_rotv_sl2

  private
  public :: subspace_rotv_sl2
  public :: d_subspace_rotv_sl2
  public :: z_subspace_rotv_sl2

  interface subsapce_rotv_sl2
    module procedure d_subspace_rotv_sl2, z_subspace_rotv_sl2
  end interface

end module subspace_rotv_sl2_module
