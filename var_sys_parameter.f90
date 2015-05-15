MODULE VarSysParameter
  implicit none
#ifdef _USPP_
  character(4) :: pp_kind='USPP'
#else
  character(4) :: pp_kind='NCPP'
#endif

  logical :: isRootRank
  logical :: isTestRank
  integer :: testRank=1 ! this will show testRank status
#ifdef _ParallelTest_
  logical :: isParallelTest=.true.
#else
  logical :: isParallelTest=.false.
#endif

CONTAINS
  SUBROUTINE setDispSwitch(myrank,nprocs)
    implicit none
    integer,intent(IN) :: myrank
    integer,intent(IN) :: nprocs
    isRootRank = (myrank==0)
    isTestRank = (myrank==testRank)
    if (isParallelTest .and. nprocs>9) then
      stop 'nprocs>9 is not suitable for parallel test'
    endif
    return
  END SUBROUTINE setDispSwitch
    

END MODULE VarSysParameter
