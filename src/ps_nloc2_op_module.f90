module ps_nloc2_op_module

  implicit none

  private
  public :: reset_init_flag_ps_nloc2_op

contains

  subroutine reset_init_flag_ps_nloc2_op
    use d_ps_nloc2_op_module, only: d_reset_init_flag_ps_nloc2_op
    use z_ps_nloc2_op_module, only: z_reset_init_flag_ps_nloc2_op
    implicit none
    call d_reset_init_flag_ps_nloc2_op
    call z_reset_init_flag_ps_nloc2_op
  end subroutine reset_init_flag_ps_nloc2_op

end module ps_nloc2_op_module
