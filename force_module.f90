MODULE force_module

  use ps_local_module
  use ps_nloc2_module
  use ps_pcc_module
  use force_ewald_module
  use watch_module
  use parallel_module, only: disp_switch_parallel

  implicit none

  PRIVATE
  PUBLIC :: calc_force

CONTAINS

  SUBROUTINE calc_force(MI,force)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),allocatable :: work(:,:)
    real(8) :: ctt(0:3),ett(0:3)

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
       call calc_force_ps_pcc(MI,work)
       force = force + work
    end if

    call watch(ctt(1),ett(1))

    call calc_force_ps_nloc2(MI,work)
    force = force + work

    call watch(ctt(2),ett(2))

    call calc_force_ewald(MI,work)
    force = force + work

    call watch(ctt(3),ett(3))

! --- constraint & symmetry ---

!    call_symforce
!    call_constraint

    deallocate( work )

    if ( disp_switch_parallel ) then
       write(*,*) "time(force1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE calc_force

END MODULE force_module
