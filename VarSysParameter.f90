MODULE VarSysParameter
  implicit none
#ifdef _USPP_
  character(4) :: pp_kind='USPP'
#else
  character(4) :: pp_kind='NCPP'
#endif

END MODULE VarSysParameter
