MODULE io_ctrl_parameters

  use io_tools_module

  implicit none

  integer :: IO_ctrl=0
  integer :: IC=0
  integer :: OC=3
  integer :: OC2=100
  integer :: MBwr1=0
  integer :: MBwr2=0
  integer :: icount=0

  integer,allocatable :: lat_new(:,:)
  integer,allocatable :: lat_old(:,:)

  logical :: flag_overwrite=.true. 

CONTAINS

  SUBROUTINE read_io_ctrl_parameters
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
       flag_overwrite = .false.
       IO_ctrl = IO_ctrl - 10
    end if
  END SUBROUTINE read_io_ctrl_parameters

  SUBROUTINE deallocate_io_ctrl_parameters
    implicit none
    if ( allocated(lat_new) ) deallocate(lat_new)
    if ( allocated(lat_old) ) deallocate(lat_old)
  END SUBROUTINE deallocate_io_ctrl_parameters

END MODULE io_ctrl_parameters
