MODULE pseudopotentials
    use strfac_module, only: construct_strfac,destruct_strfac
    use ps_local_module, only: init_ps_local,construct_ps_local,construct_ps_local_ffte
    use ps_pcc_module, only: init_ps_pcc,construct_ps_pcc
    use ps_initrho_module, only: init_ps_initrho,construct_ps_initrho
    use density_module, only: normalize_density
    use ps_nloc2_init_module, only: ps_nloc2_init
    use ps_nloc2_module, only: prep_ps_nloc2
    use ps_nloc3_module, only: init_ps_nloc3,prep_ps_nloc3
    use ps_nloc_mr_module, only: prep_ps_nloc_mr
    use pseudopot_module, only: pselect
    implicit none

CONTAINS

    SUBROUTINE initiationPS(gcut)
        implicit none
        real(8),intent(IN) :: gcut
        
        call init_ps_local
        call init_ps_pcc
        call init_ps_initrho
        call construct_strfac

#ifndef _FFTE_
        call construct_ps_local
#else
        call construct_ps_local_ffte
#endif
        
        if (pselect/=4 .and. pselect/=5) then
            call construct_ps_pcc
            call construct_ps_initrho
            call normalize_density
        end if

        call destruct_strfac

        select case( pselect )
        case( 2 )
            call ps_nloc2_init(gcut)
            call prep_ps_nloc2
        case( 3 )
            call init_ps_nloc3
            call prep_ps_nloc3
        case( 5 )
            call prep_ps_nloc_mr

#ifdef _USPP_
        case( 102 )
            call initKtoKPSQ
            call ps_nloc2_init(gcut)
            call ps_Q_init(gcut)
            call prep_ps_nloc2
#endif

        end select

        return
    END SUBROUTINE initiationPS

END MODULE pseudopotentials
