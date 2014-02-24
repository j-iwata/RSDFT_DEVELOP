!-----------------------------------------------------------------------
!     Evaluate Force BOMD
!-----------------------------------------------------------------------
SUBROUTINE getforce
  use bb_module
  use atom_module, only: Natom,aa_atom
  use strfac_module
  use ps_local_module
  use ps_pcc_module
  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use ewald_module
  use force_module
  use cpmd_variables, only: Rion,Force,pmass,iatom,AMU,disp_switch &
       ,MB_0_CPMD,MB_1_CPMD,MB_0_SCF,MB_1_SCF
  use scf_module
  use array_bound_module, only: MB_0,MB_1
  implicit none
  real(8) :: c
  integer :: a,Diter1,iter_final

  Diter1       = 100
  c            = 1.d0/(2.d0*acos(-1.d0))
  aa_atom(:,:) = matmul(transpose(bb),Rion)*c
  Force(:,:)   = 0.0d0

  call init_ion

  call calc_ewald(Eewald,disp_switch)

  call construct_strfac
#ifndef _FFTE_
  call construct_ps_local
#else
  call construct_ps_local_ffte
#endif
  call construct_ps_pcc
  call destruct_strfac

  select case(pselect)
  case(2)
     call prep_ps_nloc2
  case(3)
     call prep_ps_nloc3
  case(5)
     call prep_ps_nloc_mr
  end select

  MB_0=MB_0_SCF
  MB_1=MB_1_SCF

  call calc_scf(Diter1,0,iter_final,disp_switch)

  MB_0=MB_0_CPMD
  MB_1=MB_1_CPMD

  call calc_force(Natom,Force)

  do a=1,Natom
     Force(:,a)=Force(:,a)/(pmass(iatom(a))*AMU)
  enddo

  return
END SUBROUTINE getforce
