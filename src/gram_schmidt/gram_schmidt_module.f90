MODULE gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_t_module
  use gram_schmidt_t3_module
  use gram_schmidt_lusl_module
  use gram_schmidt_luslbp_module
 !use gram_schmidt_u_module
  use gram_schmidt_g_module
  use gram_schmidt_ncol_module, only: gram_schmidt_ncol, flag_noncollinear
  use var_sys_parameter
  use io_tools_module
  use wf_module, only: gather_b_wf, unk
  use watch_module
  use parallel_module, only: comm_grid, comm_band

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt

  integer :: iswitch_algorithm = 0
  include 'mpif.h'

  logical :: flag_init_read = .true.

  integer :: iparam_gs(9)=1

CONTAINS

  SUBROUTINE read_gram_schmidt
    implicit none
    call IOTools_readIntegerKeyword( "GS", iswitch_algorithm )
    call IOTools_readIntegerKeyword( "GSPARAM", iparam_gs )
    flag_init_read = .false.
  END SUBROUTINE read_gram_schmidt

  SUBROUTINE gram_schmidt(n0,n1,k,s)

    implicit none
    integer,intent(IN) :: n0,n1,k,s
    type(time) :: t

    if ( flag_noncollinear ) then
       call gram_schmidt_ncol( n0,n1,k,unk )
       return
    end if

!    call write_border( 1, " gram_schmidt(start)" )

!    call start_timer( t )

    if ( flag_init_read ) call read_gram_schmidt

    call gather_b_wf( k, s )

    if ( pp_kind == "USPP" ) then

       call gram_schmidt_g( n0,n1,k,s )

    else

       select case( iswitch_algorithm )
       case default
          call gram_schmidt_t(n0,n1,k,s)
       case( 1 )
          call gram_schmidt_m(n0,n1,k,s)
      !case( 2 )
      !   call gram_schmidt_u(n0,n1,k,s)
      case( 3 )
         call gram_schmidt_t3( unk(:,:,k,s), comm_grid, comm_band )
      case( 4 )
      !   do i=1,4
      !      if ( iparam_gs(i) /= 0 ) iparam_gs_lusl(i)=iparam_gs(i)
      !   end do
         call gram_schmidt_lusl( unk(:,:,k,s) )
      case( 5 )
#ifdef _DRSDFT_
         call gram_schmidt_luslbp( unk(:,:,k,s) )
#else
         write(*,*) "z_gram_schmidt_lusl is not implemented yet"
         call stop_program('gram_schmidt')
#endif
      end select

    end if

!    call result_timer( "gs", t )

!    call write_border( 1, " gram_schmidt(end)" )

  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module
