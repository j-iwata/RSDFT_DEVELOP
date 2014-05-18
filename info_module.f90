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
    call header_info
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

  SUBROUTINE header_info
    implicit none
    character(8)  :: date
    character(10) :: time
    date=""
    time=""
    call date_and_time(DATE=date,TIME=time)
    if ( myrank == 0 ) then
       write(*,*) "(Subversion Revision) $Rev$"
       write(*,*) "(Subversion Date    ) $Date$"
       write(*,*) "(Compile date/time) ",date,"/",time
       write(unit_info,*) "(Subversion Revision) $Rev$"
       write(unit_info,*) "(Subversion Date    ) $Date$"
       write(unit_info,*) "(Compile date/time) ",date,"/",time
#ifdef _DRSDFT_
       write(*,*) "DRSDFT(REAL8)"
       write(unit_info,*) "DRSDFT(REAL8)"
#else
       write(*,*) "ZRSDFT(COMPLEX16)"
       write(unit_info,*) "ZRSDFT(COMPLEX16)"
#endif
#ifdef _SPLINE_
       write(*,*) "SPLINE(ps_nloc2_module)"
       write(unit_info,*) "SPLINE(ps_nloc2_module)"
#endif
#ifdef _FFTE_
       write(*,*) "FFTE(ps_local,hartree)"
       write(unit_info,*) "FFTE(ps_local,hartree)"
#endif
#ifdef _LAPACK_
       write(*,*) "LAPACK"
       write(unit_info,*) "LAPACK"
#endif
    end if
  END SUBROUTINE header_info

END MODULE info_module
