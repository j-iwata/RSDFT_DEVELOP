MODULE pseudopotentials

  use strfac_module, only: construct_strfac,destruct_strfac
  use ps_local_module, only: init_ps_local,construct_ps_local,construct_ps_local_ffte
  use ps_pcc_module, only: init_ps_pcc,construct_ps_pcc
  use ps_initrho_module, only: init_ps_initrho,construct_ps_initrho
  use density_module, only: normalize_density
  use ps_nloc2_init_module, only: ps_nloc2_init,rcfac,qcfac,etafac
  use ps_nloc2_module, only: prep_ps_nloc2
  use ps_nloc3_module, only: init_ps_nloc3,prep_ps_nloc3
  use ps_nloc_mr_module, only: prep_ps_nloc_mr
  use pseudopot_module, only: pselect
#ifdef _USPP_
  use PSQInit, only: initKtoKPSQ,ps_Q_init
  use PSQRijPrep, only: prepQRijp102
  use PSnonLocPrepG, only: prepNzqr
#endif
  use parallel_module, only: myrank
  use VarSysParameter, only: pp_kind
  use VarPSMember, only: ippform

  implicit none

CONTAINS

  SUBROUTINE initiatePS(gcut)
    implicit none
    real(8),intent(IN) :: gcut

!    call init_ps_local
!    call init_ps_pcc
!    call init_ps_initrho
!    call construct_strfac
!#ifndef _FFTE_
!    call construct_ps_local
!#else
!    call construct_ps_local_ffte
!#endif
!    if (pselect/=4 .and. pselect/=5) then
!      call construct_ps_pcc
!      call construct_ps_initrho
!      call normalize_density
!    end if
!    call destruct_strfac

!----

    if ( all(ippform < 100) ) then
       pp_kind="NCPP"
    else
       pp_kind="USPP"
    end if

!----------------------------------------- PSELECT
    select case( pselect )

!    case( 2 )
!       call ps_nloc2_init(gcut)
!       call prep_ps_nloc2
!    case( 3 )
!       call init_ps_nloc3
!       call prep_ps_nloc3
!    case( 5 )
!       call prep_ps_nloc_mr

#ifdef _USPP_
    case( 102 )

       if ( myrank == 0 ) write(200,*) "----- init PS nonlocal for pselect==102 -----"
       if ( myrank == 0 ) write(200,*) "----- initKtoKPSQ -----"
       call initKtoKPSQ
       if ( myrank == 0 ) write(200,*) "----- ps_nloc2_init -----"
       call ps_nloc2_init(gcut)
       if ( myrank == 0 ) write(200,*) "----- ps_Q_init -----"
       call ps_Q_init(gcut,rcfac,qcfac,etafac)
       if ( myrank == 0 ) write(200,*) "----- prep_ps_nloc2 -----"
       call prep_ps_nloc2
       if ( myrank == 0 ) write(200,*) "----- prepNzqr -----"
       call prepNzqr
       if ( myrank == 0 ) write(200,*) "----- prepQRijp102 -----"
       call prepQRijp102
#endif

    end select
!========================================= PSELECT

    return
  END SUBROUTINE initiatePS

END MODULE pseudopotentials
