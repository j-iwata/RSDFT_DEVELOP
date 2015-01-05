!-----------------------------------------------------------------------
!     Evaluate Force for CPMD
!-----------------------------------------------------------------------
SUBROUTINE getforce_cpmd(lewald,ltime)

  use eion_module, only: calc_eion
  use atom_module, only: Natom,aa_atom
  use bb_module
  use parallel_module, only: myrank
  use strfac_module
  use ps_local_module
  use ps_pcc_module
  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use localpot_module, only: Vloc
  use array_bound_module, only: MSP_0,MSP_1,MB_0,MB_1,ML_0,ML_1
  use density_module
  use xc_module
  use hartree_variables, only: Vh
  use hartree_module, only: calc_hartree
  use cpmd_variables, only: Force,Rion,disp_switch,AMU,pmass,iatom
  use watch_module
  use force_module

  use PSnonLocDij
  use PSnonLocPrepG
  use PSQRijPrep

  implicit none

  integer :: i,s,MB_0_BAK,MB_1_BAK
  real(8) :: ctime_force(0:9),etime_force(0:9),c
  logical :: lewald,ltime
  logical,save :: first_time=.true.
  integer,save :: icount=0

  c=1.d0/(2.d0*acos(-1.d0))
  aa_atom = matmul(transpose(bb),Rion)*c
  Force   = 0.0d0

  call init_ion

  if ( ltime ) call watch(ctime_force(0),etime_force(0))

  if ( lewald ) then
     call calc_eion
  end if

  if ( ltime ) call watch(ctime_force(1),etime_force(1))

  call construct_strfac
  call construct_ps_local
  call construct_ps_pcc
  call destruct_strfac

  if ( ltime ) call watch(ctime_force(2),etime_force(2))

  select case(pselect)
  case(2)
     call prep_ps_nloc2
  case(3)
     call prep_ps_nloc3
  case(5)
     call prep_ps_nloc_mr
#ifdef _USPP_
  case(102)
     call prep_ps_nloc2
     call prepNzqr
     call prepQRijp102
#endif
  end select

  if ( ltime ) call watch(ctime_force(3),etime_force(3))

  call calc_density

  if ( ltime ) call watch(ctime_force(4),etime_force(4))

  call calc_hartree(ML_0,ML_1,MSP_1-MSP_0+1,rho(ML_0,MSP_0))

  if ( ltime ) call watch(ctime_force(5),etime_force(5))

  call calc_xc

  if ( ltime ) call watch(ctime_force(6),etime_force(6))

  do s=MSP_0,MSP_1
     Vloc(:,s)=Vh(:)+Vxc(:,s)+Vion(:)
  enddo

#ifdef _USPP_
       call getDij
#endif
  if ( ltime ) call watch(ctime_force(7),etime_force(7))

  call calc_force(Natom,Force)

  if ( ltime ) call watch(ctime_force(8),etime_force(8))

  if ( myrank == 0 ) then
     if ( first_time ) open(33,file="force.dat")
     icount=icount+1
     write(33,*) "step=",icount
     do i=1,Natom
        write(33,'(1x,i8,3f20.15)') i,Force(1:3,i)
     end do
  end if

  do i=1,Natom
     Force(:,i)=Force(:,i)/(pmass(iatom(i))*AMU)
  enddo

  if ( ltime ) call watch(ctime_force(9),etime_force(9))

  if ( ltime .and. myrank==0 ) then
     write(17,'(9f10.5)') (etime_force(i+1)-etime_force(i),i=0,8)
  endif

  return

END SUBROUTINE getforce_cpmd
