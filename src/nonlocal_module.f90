MODULE nonlocal_module

  use pseudopot_module, only: pselect, ps_type, flag_so
  use ps_nloc2_op_module, only: op_ps_nloc2, op_ps_nloc2_hp
  use ps_nloc2_module, only: calc_force_ps_nloc2, prep_uvk_ps_nloc2 &
                            ,prep_rvk_ps_nloc2
  use force_ps_nloc2_module, only: calc_force_ps_nloc2_hp
  use ps_nloc3_module
  use ps_nloc_mr_module
  use PSnonLocOpG2
  use parallel_module, only: myrank, ir_grid, id_grid, myrank_g

  implicit none

  PRIVATE
  PUBLIC :: op_nonlocal
  PUBLIC :: calc_force_nonlocal
  PUBLIC :: update_k_dependence_nonlocal

CONTAINS

  SUBROUTINE op_nonlocal( tpsi, htpsi, k_in, s_in, htpsi00 )

    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8),intent(INOUT),optional :: htpsi00(:,:)
#else
    complex(8),intent(IN)  :: tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8),intent(INOUT),optional :: htpsi00(:,:)
#endif
    integer,optional,intent(IN) :: k_in, s_in
    integer :: k,s,n1,n2
    k=1 ; if ( present(k_in) ) k=k_in
    s=1 ; if ( present(s_in) ) s=s_in
    select case( pselect )
    case( 2 )
       if ( ps_type == 1 ) then
          call op_ps_nloc_mr( tpsi, htpsi, k )
       else
          call op_ps_nloc2_hp( tpsi, htpsi, k )
          !n1=id_grid(myrank_g)+1
          !n2=n1+ir_grid(myrank_g)-1
          !call op_ps_nloc2( k, tpsi, htpsi, n1, n2, 1, size(tpsi,2) )
       end if
    case(3)
       call op_ps_nloc3( tpsi, htpsi, k )
    case(102)
      call op_ps_nloc2_uspp( tpsi, htpsi, htpsi00, k, s )
    end select

  END SUBROUTINE op_nonlocal


  SUBROUTINE calc_force_nonlocal( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    select case( pselect )
    case( 2 )
       !call calc_force_ps_nloc2( MI, force )
       call calc_force_ps_nloc2_hp( MI, force )
    case( 3 )
       call calc_force_ps_nloc3( MI, force )
    case( 102 )
       call stop_program("stop@force_nonlocal(2)")
    end select
  END SUBROUTINE calc_force_nonlocal


  SUBROUTINE update_k_dependence_nonlocal( MBZ_0, MBZ_1, kbb, flag_momentum )
    implicit none
    integer,intent(IN) :: MBZ_0, MBZ_1
    real(8),intent(IN) :: kbb(:,:)
    logical,optional,intent(IN) :: flag_momentum
    select case( pselect )
    case( 2 )
       if ( ps_type == 1 ) then
          call prep_uvk_ps_nloc_mr(MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1))
          if ( flag_so ) &
               call prep_uvkso_ps_nloc_mr(MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1))
       else
          call prep_uvk_ps_nloc2(MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1))
       end if
    case( 3 )
       call init_ps_nloc3
       call prep_ps_nloc3
    end select
    if ( present(flag_momentum) .and. flag_momentum ) then
       select case( pselect )
       case( 2 )
          if ( ps_type == 1 ) then
             call prep_rvk_ps_nloc_mr( MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1) )
          else
             call prep_rvk_ps_nloc2( MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1) )
          end if
       case( 3 )
          !call prep_rvk_ps_nloc3(MBZ_0,MBZ_1,kbb(1,MBZ_0:MBZ_1) )
       end select
    end if
  END SUBROUTINE update_k_dependence_nonlocal


END MODULE nonlocal_module
