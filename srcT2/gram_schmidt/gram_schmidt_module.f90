MODULE gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_m_bp_module
  use gram_schmidt_t_module
 !use gram_schmidt_u_module
  use gram_schmidt_g_module
  use gram_schmidt_ncol_module, only: gram_schmidt_ncol, flag_noncollinear
  use var_sys_parameter, only: pp_kind
  use io_tools_module
  use wf_module, only: gather_b_wf, unk
  use watch_module
  use gram_schmidt_t2_module, only: gram_schmidt_t2, read_gram_schmidt_t2
  use gram_schmidt_h_module
  use gram_schmidt_hbp_module
  use gram_schmidt_pbp_module     ! diag first
  use gram_schmidt_p_module
  use vector_tools_module
  use rgrid_module, only: dV
  use reconst_module
  use lugs_sl_test_module
  use lugs_sl_bp_test_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt

  integer :: iswitch_algorithm = 10
  logical :: flag_init_read = .true.

CONTAINS

  SUBROUTINE read_gram_schmidt
    implicit none
    call IOTools_readIntegerKeyword( "GS", iswitch_algorithm )
    flag_init_read = .false.
  END SUBROUTINE read_gram_schmidt

  SUBROUTINE gram_schmidt(n0,n1,k,s,v)

    implicit none
    integer,intent(IN) :: n0,n1,k,s
    type(vinfo),intent(IN) :: v(2)
    integer :: nblk,mbt,mb0,ml0,np_band
    type(time) :: t
    real(8),allocatable :: u(:,:)

!    if ( first_time ) then
!       v(1)%factor=dV
!       v(1)%pinfo%comm=comm_grid
!       v(1)%pinfo%np=np_grid
!       v(1)%pinfo%me=myrank_g
!       allocate( v(1)%pinfo%ir(0:np_grid-1) ) ; v(1)%pinfo%ir=ir_grid
!       allocate( v(1)%pinfo%id(0:np_grid-1) ) ; v(1)%pinfo%id=id_grid
!       v(2)%factor=1.0d0
!       v(2)%pinfo%comm=comm_band
!       v(2)%pinfo%np=np_band
!       v(2)%pinfo%me=myrank_b
!       allocate( v(2)%pinfo%ir(0:np_band-1) ) ; v(2)%pinfo%ir=ir_band
!       allocate( v(2)%pinfo%id(0:np_band-1) ) ; v(2)%pinfo%id=id_band
!       first_time=.false.
!    end if

    if ( flag_noncollinear ) then
       call gram_schmidt_ncol( n0,n1,k,unk )
       return
    end if

    call write_border( 1, " gram_schmidt(start)" )
    call start_timer( t )

    if ( flag_init_read ) call read_gram_schmidt

    np_band = v(2)%pinfo%np

    if ( pp_kind == "USPP" ) then

       call gram_schmidt_g( n0,n1,k,s )

    else

       select case( iswitch_algorithm )
       case default
          if ( np_band > 1 ) goto 900
          call gram_schmidt_m(n0,n1,k,s)
       case( 10 )
          call gram_schmidt_m_bp( unk(:,:,k,s), v )
       case( 1 )
          if ( np_band > 1 ) goto 900
          call gram_schmidt_t(n0,n1,k,s)
       case( 2 )
          if ( np_band > 1 ) goto 900
          call read_gram_schmidt_t2
          call gram_schmidt_t2( unk(:,:,k,s), v )
       case( 3 )
          if ( np_band > 1 ) goto 900
          call read_gram_schmidt_h
          call gram_schmidt_h( unk(:,:,k,s), v )
!      case( 11,12,13 )
       case( 11,12 )
          call read_gram_schmidt_hbp( nblk )
          call gram_schmidt_hbp( unk(:,:,k,s), v )
       case( 13 )
          call read_gram_schmidt_pbp( nblk )
          call gram_schmidt_pbp( unk(:,:,k,s), v )
       case( 4 )
          if ( np_band > 1 ) goto 900
          call LUGS_sl_test( unk(:,:,k,s), v )
       case( 14 )
          call LUGS_sl_bp_test( unk(:,:,k,s), v )
       case( 5 )
          if ( np_band > 1 ) goto 900
          call read_gram_schmidt_p
          call gram_schmidt_p( unk(:,:,k,s), v )
       case( 6,16 )
          call stop_program("QRGS is not available")
          !call QRGS_sl_test()
       end select

    end if

    call result_timer( "gs", t )

    call write_border( 1, " gram_schmidt(end)" )

    return

900 call stop_program("This routine is incompatible with memory-band-parallel")

  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module
