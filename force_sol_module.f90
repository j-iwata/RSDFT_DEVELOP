MODULE force_sol_module

  use ps_local_module, only: calc_force_ps_local
  use ps_nloc2_module, only: calc_force_ps_nloc2
  use ps_pcc_module, only: flag_pcc_0
  use ps_pcc_force_module, only: calc_ps_pcc_force
  use force_ewald_module, only: calc_force_ewald
  use watch_module
  use parallel_module, only: disp_switch_parallel
  use vdw_grimme_module, only: calc_F_vdw_grimme
  use atom_module, only: aa_atom

  implicit none

  PRIVATE
  PUBLIC :: calc_force_sol

CONTAINS

  SUBROUTINE calc_force_sol( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),allocatable :: work(:,:)
    real(8) :: ctt(0:4),ett(0:4)

    force(:,:) = 0.d0

    ctt(:)=0.d0
    ett(:)=0.d0

    allocate( work(3,MI) )

    call watch(ctt(0),ett(0))

#ifdef _FFTE_
    call calc_force_ps_local_ffte(MI,work)
#else
    call calc_force_ps_local(MI,work)
#endif
    force = force + work

    if ( flag_pcc_0 ) then
       call calc_ps_pcc_force( MI, work )
       force = force + work
    end if

    call watch(ctt(1),ett(1))

    call calc_force_ps_nloc2(MI,work)
    force = force + work

    call watch(ctt(2),ett(2))

    call calc_force_ewald(MI,work)
    force = force + work

    call watch(ctt(3),ett(3))

    call calc_F_vdw_grimme( MI, aa_atom, force )

    call watch(ctt(4),ett(4))

    deallocate( work )

    if ( disp_switch_parallel ) then
       write(*,*) "time(force1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(force4)",ctt(4)-ctt(3),ett(4)-ett(3)
    end if

  END SUBROUTINE calc_force_sol

END MODULE force_sol_module
