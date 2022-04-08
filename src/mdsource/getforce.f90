!-----------------------------------------------------------------------
!     Evaluate Force BOMD
!-----------------------------------------------------------------------
SUBROUTINE getforce
  use bb_module
  use atom_module, only: Natom,aa_atom,ki_atom,zn_atom,shift_aa_coordinates_atom
  use strfac_module
  use ps_local_module
  use ps_pcc_module
  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use eion_module, only: calc_eion
  use force_module
  use cpmd_variables, only: Rion,Force,pmass,AMU,disp_switch &
       ,MB_0_CPMD,MB_1_CPMD,MB_0_SCF,MB_1_SCF,lcpmd
  use scf_module
  use array_bound_module, only: MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1
  use parallel_module, only: end_mpi_parallel, disp_switch_parallel
  use ps_nloc2_variables, only: prep_backup_uVunk_ps_nloc2
  use io_tools_module

  implicit none
  real(8) :: c
  integer :: a,ierr
  integer,save :: Diter1=0

  if ( Diter1 == 0 ) then
    Diter1 = 100
    call IOTools_readIntegerKeyword( "DITER", Diter1 )
  end if
  c            = 1.d0/(2.d0*acos(-1.d0))
  aa_atom(:,:) = matmul(transpose(bb),Rion)*c
  call shift_aa_coordinates_atom( aa_atom )

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
  end select

  MB_0=MB_0_SCF
  MB_1=MB_1_SCF

  call calc_scf( ierr, Diter1 )
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

  if ( lcpmd ) call prep_backup_uVunk_ps_nloc2( MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1 )

  Force(:,:)=0.0d0
  call calc_force(Natom,Force)

  do a=1,Natom
     Force(:,a)=Force(:,a)/(pmass(zn_atom(ki_atom(a)))*AMU)
  enddo

  return
END SUBROUTINE getforce
