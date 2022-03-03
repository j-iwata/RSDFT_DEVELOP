module gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_t_module
  use gram_schmidt_t3_module
  use gram_schmidt_lusl_module
  use gram_schmidt_luslbp_module
  use gram_schmidt_ncol_module, only: gram_schmidt_ncol, flag_noncollinear
  use var_sys_parameter
  use io_tools_module
  use watch_module
  use parallel_module, only: comm_grid, comm_band

  implicit none

  private
  public :: gram_schmidt

  integer :: iswitch_algorithm = 0
  include 'mpif.h'

  logical :: flag_init_read = .true.

  integer :: iparam_gs(9)=1

contains

  subroutine read_gram_schmidt
    implicit none
    call IOTools_readIntegerKeyword( "GS", iswitch_algorithm )
    call IOTools_readIntegerKeyword( "GSPARAM", iparam_gs )
    flag_init_read = .false.
  end subroutine read_gram_schmidt

  subroutine gram_schmidt(n0,n1,k,s)
    use rgrid_variables, only: dV
    use wf_module, only: gather_b_wf, unk
    implicit none
    integer,intent(in) :: n0,n1,k,s
    type(time) :: t

    if ( flag_noncollinear ) then
       call gram_schmidt_ncol( n0,n1,k,unk )
       return
    end if

!    call write_border( 1, " gram_schmidt(start)" )

!    call start_timer( t )

    if ( flag_init_read ) call read_gram_schmidt

    call gather_b_wf( k, s )

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
      call gram_schmidt_luslbp( unk(:,:,k,s), dV )
    end select

    ! call gather_b_wf( k, s )

    ! call result_timer( "gs", t )

!    call write_border( 1, " gram_schmidt(end)" )

  end subroutine gram_schmidt

end module gram_schmidt_module
