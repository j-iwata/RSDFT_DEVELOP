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
  use eion_module, only: calc_eion
  use force_module
  use cpmd_variables, only: Rion,Force,pmass,iatom,AMU,disp_switch &
       ,MB_0_CPMD,MB_1_CPMD,MB_0_SCF,MB_1_SCF
  use scf_module
  use array_bound_module, only: MB_0,MB_1
  use parallel_module, only: end_mpi_parallel, disp_switch_parallel

  use ps_prepNzqr_g_module
  use ps_qrij_prep_module

  implicit none
  real(8) :: c
  integer :: a,Diter1,ierr

  Diter1       = 100
  c            = 1.d0/(2.d0*acos(-1.d0))
  aa_atom(:,:) = matmul(transpose(bb),Rion)*c
  Force(:,:)   = 0.0d0

  call init_ion

  call calc_eion

  call construct_strfac
  call construct_ps_local
  call construct_ps_pcc
  call destruct_strfac

  select case(pselect)
  case(2)
     call prep_ps_nloc2
  case(3)
     call prep_ps_nloc3
  case(5)
     call prep_ps_nloc_mr
  case(102)
     call prep_ps_nloc2
     call prepNzqr
     call prepQRijp102
  end select

  MB_0=MB_0_SCF
  MB_1=MB_1_SCF

  call calc_scf( disp_switch_parallel, ierr, Diter1 )
  if ( ierr == -1 ) then
     if ( disp_switch ) write(*,*) "time limit !!!"
     call end_mpi_parallel
     stop "stop@getforce"
  end if
  if ( ierr == -2 ) then
     if ( disp_switch ) write(*,*) "SCF is not converged"
  end if

  MB_0=MB_0_CPMD
  MB_1=MB_1_CPMD

  call calc_force(Natom,Force)

  do a=1,Natom
     Force(:,a)=Force(:,a)/(pmass(iatom(a))*AMU)
  enddo

  return
END SUBROUTINE getforce
