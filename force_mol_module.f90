MODULE force_mol_module

  use force_local_mol_module
  use force_nloc2_mol_module
  use force_ion_mol_module
  use watch_module
  use parallel_module, only: disp_switch_parallel

  implicit none

  PRIVATE
  PUBLIC :: calc_force_mol

CONTAINS

  SUBROUTINE calc_force_mol(MI,force)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),allocatable :: work(:,:)
    real(8) :: ctt(0:3),ett(0:3)

    force(:,:) = 0.0d0

    ctt(:)=0.d0
    ett(:)=0.d0

    allocate( work(3,MI) ) ; work=0.0d0

    call watch(ctt(0),ett(0))

    call calc_force_local_mol(work)
    force = force + work

!    if ( flag_pcc_0 ) then
!       write(*,*) "RSMOL: force with pcc is not available yet"
!       stop "stop@force_mol"
!    end if

    call watch(ctt(1),ett(1))

    call calc_force_nloc2_mol(work)
    force = force + work

    call watch(ctt(2),ett(2))

    call calc_force_ion_mol(work)
    force = force + work

    call watch(ctt(3),ett(3))

! --- constraint & symmetry ---

!    call_symforce
!    call_constraint

    deallocate( work )

    if ( disp_switch_parallel ) then
       write(*,*) "time(force_mol1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force_mol2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force_mol3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE calc_force_mol

END MODULE force_mol_module
