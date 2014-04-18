MODULE info_module

  implicit none

  PRIVATE
  PUBLIC :: open_info, write_info, close_info

  integer,parameter :: unit_info = 100
  character(30) :: file_name = "RSDFT_INFO"
  integer :: myrank

CONTAINS

  SUBROUTINE open_info(rank)
    implicit none
    integer,intent(IN) :: rank
    myrank = rank
    if ( rank == 0 ) then
       open(unit_info,file=file_name)
    end if
  END SUBROUTINE open_info

  SUBROUTINE write_info(info)
    implicit none
    character(*),intent(IN) :: info
    if ( info == "" ) return
    if ( myrank == 0 ) then
       write(*,*) "INFO: ",info
       write(unit_info,*) "INFO: ",info
    end if
  END SUBROUTINE write_info

  SUBROUTINE close_info
    if ( myrank == 0 ) close(unit_info)
  END SUBROUTINE close_info

END MODULE info_module
