module kinetic_sol_module

  use d_kinetic_sol_module, only: d_op_kinetic_sol, d_construct_matrix_kinetic_sol
  use z_kinetic_sol_module, only: z_op_kinetic_sol, z_construct_matrix_kinetic_sol

  implicit none

  private
  public :: op_kinetic_sol
  public :: construct_matrix_kinetic_sol

  interface op_kinetic_sol
    module procedure d_op_kinetic_sol, z_op_kinetic_sol
  end interface

  interface construct_matrix_kinetic_sol
    module procedure d_construct_matrix_kinetic_sol, z_construct_matrix_kinetic_sol
  end interface

end module kinetic_sol_module
