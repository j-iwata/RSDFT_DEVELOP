module var_sys_parameter

  implicit none

  private
  public :: use_real8_wf
  public :: systype_query

  character(4), public :: pp_kind = "NCPP"

  logical :: flag_real8 = .false.
  integer :: SYStype = 0

contains

  logical function use_real8_wf( flag )
    implicit none
    logical,optional,intent(in) :: flag
    if ( present(flag) ) flag_real8 = flag
    use_real8_wf = flag_real8
  end function use_real8_wf

  integer function systype_query( set_value )
    implicit none
    integer,optional,intent(in) :: set_value
    if ( present(set_value) ) SYStype = set_value
    systype_query = SYStype
  end function systype_query

end module var_sys_parameter
