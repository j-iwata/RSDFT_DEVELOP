module io_ctrl_parameters

  implicit none

  integer :: IO_ctrl=3
  integer :: IC=0
  integer :: OC=3
  integer :: OC2=100
  integer :: MBwr1=0
  integer :: MBwr2=0
  integer :: icount=0

  integer,allocatable :: lat_new(:,:)
  integer,allocatable :: lat_old(:,:)

  logical :: flag_overwrite_io = .true.
  logical,private :: FLAG_WF_AVAILABLE = .false.

  character(30) :: file_wf0   = "wf.dat1"
  character(30) :: file_vrho0 = "vrho.dat1"
  character(30) :: file_wf1   = "wf.dat1"
  character(30) :: file_vrho1 = "vrho.dat1"
  character(30) :: file_wf2   = "wf.dat1"
  character(30) :: file_vrho2 = "vrho.dat1"

contains

  subroutine read_io_ctrl_parameters
    use io_tools_module, only: IOTools_readIntegerKeyword
    implicit none
    integer :: itmp(2)
    itmp = (/ MBwr1, MBwr2 /)
    call IOTools_readIntegerKeyword( "IC"    , IC  )
    call IOTools_readIntegerKeyword( "OC"    , OC  )
    call IOTools_readIntegerKeyword( "OC2"   , OC2 )
    call IOTools_readIntegerKeyword( "IOCTRL", IO_ctrl )
    call IOTools_readIntegerKeyword( "MBWR"  , itmp )
    MBwr1=itmp(1)
    MBwr2=itmp(2)
    if ( IO_ctrl >= 10  ) then
      flag_overwrite_io = .false.
      IO_ctrl = IO_ctrl - 10
    end if
    if ( IC == 1 .or. IC == 3 ) FLAG_WF_AVAILABLE = .true.
  end subroutine read_io_ctrl_parameters

  subroutine deallocate_io_ctrl_parameters
    implicit none
    if ( allocated(lat_new) ) deallocate(lat_new)
    if ( allocated(lat_old) ) deallocate(lat_old)
  end subroutine deallocate_io_ctrl_parameters

  logical function wf_available()
    implicit none
    wf_available = FLAG_WF_AVAILABLE
  end function wf_available

end module io_ctrl_parameters
