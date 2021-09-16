MODULE ps_nloc_initiate_module

  use ps_nloc2_init_module, only: ps_nloc2_init,rcfac,qcfac,etafac
  use ps_nloc2_module, only: prep_ps_nloc2
  use ps_nloc2_variables, only: read_fmax_conv_ps_nloc2
  use ps_nloc3_module, only: init_ps_nloc3, prep_ps_nloc3
  use ps_nloc_mr_module, only: prep_ps_nloc_mr
  use pseudopot_module, only: pselect
  use var_sys_parameter, only: pp_kind
  use var_ps_member, only: ps_type

  implicit none

  PRIVATE
  PUBLIC :: ps_nloc_initiate

CONTAINS

  SUBROUTINE ps_nloc_initiate( Gcut )

    implicit none
    real(8),intent(IN) :: Gcut

    call write_border( 0, " ps_nloc_initiate(start)" )

    select case( pselect )
    case default

       pp_kind="NCPP"

       if ( pselect == 1 ) then
          call ps_nloc2_init( Gcut, no_mask=.true. )
          pselect = 2
       else
          call ps_nloc2_init( Gcut )
       end if

       call read_fmax_conv_ps_nloc2

       if ( ps_type == 0 ) then
          call prep_ps_nloc2
       else if ( ps_type == 1 ) then
          call prep_ps_nloc_mr
       end if

    case( 3 )

       pp_kind="NCPP"

       call init_ps_nloc3
       call prep_ps_nloc3

    end select

    call write_border( 0, " ps_nloc_initiate(end)" )

  END SUBROUTINE ps_nloc_initiate

END MODULE ps_nloc_initiate_module
