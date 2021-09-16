module nonlocal_module

  use pseudopot_module, only: pselect, ps_type

  implicit none

  private
  public :: op_nonlocal
  public :: calc_force_nonlocal
  public :: update_k_dependence_nonlocal

  interface op_nonlocal
    module procedure d_op_nonlocal, z_op_nonlocal
  end interface

contains

  subroutine d_op_nonlocal( tpsi, htpsi, n,k,s )
    use pseudopot_module, only: pselect
    use d_ps_nloc2_op_module, only: d_op_ps_nloc2_hp
    use ps_nloc3_module, only: d_op_ps_nloc3
    use ps_nloc_mr_module, only: d_op_ps_nloc_mr
    implicit none
    real(8),intent(in) :: tpsi(:,:)
    real(8),intent(inout) :: htpsi(:,:)
    integer,optional,intent(in) :: n, k, s
    select case( pselect )
    case( 2 )
      if ( ps_type == 1 ) then
        call d_op_ps_nloc_mr( tpsi, htpsi, n,k,s )
      else
        call d_op_ps_nloc2_hp( tpsi, htpsi, n,k,s )
      end if
    case(3)
      call d_op_ps_nloc3( tpsi, htpsi, n,k,s )
    end select
  end subroutine d_op_nonlocal

  subroutine z_op_nonlocal( tpsi, htpsi, n,k,s )
    use pseudopot_module, only: pselect
    use z_ps_nloc2_op_module, only: z_op_ps_nloc2_hp
    use ps_nloc3_module, only: z_op_ps_nloc3
    use ps_nloc_mr_module, only: z_op_ps_nloc_mr
    implicit none
    complex(8),intent(in) :: tpsi(:,:)
    complex(8),intent(inout) :: htpsi(:,:)
    integer,optional,intent(in) :: n,k,s
    select case( pselect )
    case( 2 )
      if ( ps_type == 1 ) then
        call z_op_ps_nloc_mr( tpsi, htpsi, n,k,s )
      else
        call z_op_ps_nloc2_hp( tpsi, htpsi, n,k,s )
      end if
    case(3)
      call z_op_ps_nloc3( tpsi, htpsi, n,k,s )
    end select
  end subroutine z_op_nonlocal


  subroutine calc_force_nonlocal( force )
    use force_ps_nloc2_module, only: calc_force_ps_nloc2
    use force_ps_nloc2_hp_module, only: calc_force_ps_nloc2_hp
    use ps_nloc3_module, only: calc_force_ps_nloc3
    use ps_nloc_mr_module, only: calc_force_ps_nloc_mr
    implicit none
    real(8),intent(inout) :: force(:,:)
    select case( pselect )
    case( 2 )
      if ( ps_type == 1 ) then
        call calc_force_ps_nloc_mr( force )
      else
        call calc_force_ps_nloc2( force )
        ! call calc_force_ps_nloc2_hp( force )
      end if
    case( 3 )
      call calc_force_ps_nloc3( force )
    case default
      call stop_program("stop@force_nonlocal(2)")
    end select
  end subroutine calc_force_nonlocal


  subroutine update_k_dependence_nonlocal( MBZ_0, MBZ_1, kbb, flag_momentum )
    use ps_nloc2_module, only: prep_uvk_ps_nloc2, prep_rvk_ps_nloc2
    use ps_nloc_mr_module, only: prep_uvk_ps_nloc_mr, prep_rvk_ps_nloc_mr
    use ps_nloc3_module, only: init_ps_nloc3, prep_ps_nloc3
    implicit none
    integer,intent(in) :: MBZ_0, MBZ_1
    real(8),intent(in) :: kbb(:,:)
    logical,optional,intent(in) :: flag_momentum
    select case( pselect )
    case( 2 )
      if ( ps_type == 1 ) then
        call prep_uvk_ps_nloc_mr(MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1))
        ! if ( flag_so ) call prep_uvkso_ps_nloc_mr(MBZ_0,MBZ_1,kbb(:,MBZ_0:MBZ_1))
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
  end subroutine update_k_dependence_nonlocal


end module nonlocal_module
