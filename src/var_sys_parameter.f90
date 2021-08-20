module var_sys_parameter

  implicit none

  private
  public :: use_real8_wf
 
  character(4), public :: pp_kind = "NCPP"

  logical :: flag_real8 = .false.

contains

  logical function use_real8_wf( flag )
    implicit none
    logical,optional,intent(in) :: flag
    if ( present(flag) ) flag_real8 = flag
    use_real8_wf = flag_real8
  end function use_real8_wf

end module var_sys_parameter
