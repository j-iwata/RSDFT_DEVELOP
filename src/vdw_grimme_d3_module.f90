MODULE vdw_grimme_D3_module

#ifdef _DFTD3_
  use dftd3_api
#endif

  implicit none

  PRIVATE
  PUBLIC :: calc_E_vdw_grimme_D3
  PUBLIC :: calc_F_vdw_grimme_D3

CONTAINS


  SUBROUTINE calc_E_vdw_grimme_D3( aa, aa_atom, zatom, Edisp )

    implicit none
    real(8),intent(IN)  :: aa(3,3), aa_atom(:,:)
    integer,intent(IN)  :: zatom(:)
    real(8),intent(OUT) :: Edisp
    real(8),allocatable :: coords(:,:)

#ifdef _DFTD3_

    type(dftd3_input) :: input
    type(dftd3_calc) :: dftd3

    call write_border( 1, " calc_vdw_E_grimme_D3(start)" )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! You can set input variables if you like, or just leave them on their
  ! defaults, which are the same as the dftd3 program uses.
  
  !! Threebody interactions (default: .false.)
  !  input%threebody = .true.
  !
  !! Numerical gradients (default: .false.)
  !  input%numgrad = .false.
  !
  !! Cutoffs (below you find the defaults)
  !  input%cutoff = sqrt(9000.0_wp)
  !  input%cutoff_cn = sqrt(1600.0_wp)

    call dftd3_init(dftd3, input)

  !call dftd3_set_functional(dftd3, func='pbe', version=2, tz=.false.) ! -old
  !call dftd3_set_functional(dftd3, func='pbe', version=4, tz=.false.) ! -bj
   call dftd3_set_functional(dftd3, func='pbe', version=6, tz=.false.) ! -bjm

    allocate( coords(3,size(aa_atom,2)) ) ; coords=0.0d0
    coords=matmul( aa, aa_atom )

    call dftd3_pbc_dispersion( dftd3, coords, zatom, aa, Edisp )

    deallocate( coords )

    call write_border( 1, " calc_vdw_E_grimme_D3(end)" )

#else
    Edisp=0.0d0
    call stop_program("You need DFTD3 libary. (stop@calc_E_vdw_grimme_D3)") 
#endif

  END SUBROUTINE calc_E_vdw_grimme_D3


  SUBROUTINE calc_F_vdw_grimme_D3( aa, aa_atom, zatom, force )

    implicit none
    real(8),intent(IN)  :: aa(3,3), aa_atom(:,:)
    integer,intent(IN)  :: zatom(:)
    real(8),intent(INOUT) :: force(:,:)
    real(8),allocatable :: coords(:,:), grads(:,:)
    real(8) :: Ed, stress(3,3)

#ifdef _DFTD3_

    type(dftd3_input) :: input
    type(dftd3_calc) :: dftd3

    call write_border( 1, " calc_vdw_F_grimme_D3(start)" )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! You can set input variables if you like, or just leave them on their
  ! defaults, which are the same as the dftd3 program uses.
  
  !! Threebody interactions (default: .false.)
  !  input%threebody = .true.
  !
  !! Numerical gradients (default: .false.)
  !  input%numgrad = .false.
  !
  !! Cutoffs (below you find the defaults)
  !  input%cutoff = sqrt(9000.0_wp)
  !  input%cutoff_cn = sqrt(1600.0_wp)

  ! Initialize dftd3
    call dftd3_init(dftd3, input)

  ! call dftd3_set_functional(dftd3, func='pbe', version=4, tz=.false.) ! -bj
    call dftd3_set_functional(dftd3, func='pbe', version=6, tz=.false.) ! -bjm

    allocate( coords(3,size(aa_atom,2)) ) ; coords=0.0d0
    coords=matmul( aa, aa_atom )

    allocate( grads(3,size(aa_atom,2)) ) ; grads=0.0d0

    call dftd3_pbc_dispersion( dftd3, coords, zatom, aa, Ed, grads, stress )

    force=force-grads

    deallocate( grads  )
    deallocate( coords )

    call write_border( 1, " calc_vdw_F_grimme_D3(end)" )

#else
    call stop_program("You need DFTD3 libary. (stop@calc_F_vdw_grimme_D3)") 
#endif

  END SUBROUTINE calc_F_vdw_grimme_D3


END MODULE vdw_grimme_D3_module
