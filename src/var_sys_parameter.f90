module var_sys_parameter

  implicit none

  private
  public :: use_real8_wf
  public :: systype_query
  public :: xctype_query

  character(4), public :: pp_kind = "NCPP"

  logical :: flag_real8 = .false.
  integer :: SYStype = 0
  character(9) :: XCtype = "GGA"

contains

  logical function use_real8_wf( set_flag )
    implicit none
    logical,optional,intent(in) :: set_flag
    if ( present(set_flag) ) flag_real8 = set_flag
    use_real8_wf = flag_real8
  end function use_real8_wf

  integer function systype_query( set_value )
    implicit none
    integer,optional,intent(in) :: set_value
    if ( present(set_value) ) SYStype = set_value
    systype_query = SYStype
  end function systype_query

  function xctype_query( set_str )
    implicit none
    character(9) :: xctype_query
    character(9),optional,intent(in) :: set_str
    if ( present(set_str) ) XCtype = set_str
    xctype_query = XCtype
  end function xctype_query

end module var_sys_parameter
