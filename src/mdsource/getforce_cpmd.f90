!-----------------------------------------------------------------------
!     Evaluate Force for CPMD
!-----------------------------------------------------------------------
SUBROUTINE getforce_cpmd( ltime )

  use eion_module, only: calc_eion
  use atom_module, only: Natom,aa_atom,ki_atom,zn_atom
  use bb_module
  use parallel_module, only: myrank, np_band, myrank_b, comm_band
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
  use cpmd_variables, only: Force,Rion,disp_switch,AMU,pmass,ir_band_cpmd,id_band_cpmd
  use watch_module
  use force_module

  use ps_getDij_module
  use ps_prepNzqr_g_module
  use ps_qrij_prep_module

  use wf_module, only: unk
  use rsdft_mpi_module

  implicit none

  logical,intent(IN) :: ltime
  integer :: i,s
  real(8) :: ctime_force(0:9),etime_force(0:9),c

  c=1.0d0/(2.0d0*acos(-1.0d0))
  aa_atom = matmul(transpose(bb),Rion)*c
  Force   = 0.0d0

  if ( ltime ) call watch(ctime_force(0),etime_force(0))

  call calc_eion

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
  case(102)
     call prep_ps_nloc2
     call prepNzqr
     call prepQRijp102
  end select

  if ( ltime ) call watch(ctime_force(3),etime_force(3))

  if ( np_band > 1 ) call wf_gather_sub( unk )
  call calc_density

  if ( ltime ) call watch(ctime_force(4),etime_force(4))

  call calc_hartree(ML_0,ML_1,MSP_1-MSP_0+1,rho(ML_0,MSP_0))

  if ( ltime ) call watch(ctime_force(5),etime_force(5))

  call calc_xc

  if ( ltime ) call watch(ctime_force(6),etime_force(6))

  do s=MSP_0,MSP_1
     Vloc(:,s)=Vh(:)+Vxc(:,s)+Vion(:)
  enddo

  call getDij

  if ( ltime ) call watch(ctime_force(7),etime_force(7))

  call calc_force(Natom,Force)

  if ( ltime ) call watch(ctime_force(8),etime_force(8))

  do i=1,Natom
     Force(:,i)=Force(:,i)/(pmass(zn_atom(ki_atom(i)))*AMU)
  enddo

  if ( ltime ) call watch(ctime_force(9),etime_force(9))

  if ( ltime .and. myrank==0 ) then
     write(17,'(9f10.5)') (etime_force(i+1)-etime_force(i),i=0,8)
  endif

  return

CONTAINS

  SUBROUTINE wf_gather_sub( u )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: u(:,:,:,:)
#else
    complex(8),intent(INOUT) :: u(:,:,:,:)
#endif
    integer,save,allocatable :: ir(:),id(:)
    integer,save :: ib1,ib2
    integer :: k,s
    if ( .not.allocated(ir) ) then
       allocate( ir(0:np_band-1) ) ; ir=0
       allocate( id(0:np_band-1) ) ; id=0
       ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*size(u,1)
       id(0:np_band-1)=id_band_cpmd(0:np_band-1)*size(u,1)
       ib1=id_band_cpmd(myrank_b)+1
       ib2=id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)
    end if
    do s=1,size(u,4)
    do k=1,size(u,3)
       call rsdft_allgatherv( u(:,ib1:ib2,k,s),u(:,:,k,s),ir,id,comm_band )
    end do
    end do
  END SUBROUTINE wf_gather_sub

END SUBROUTINE getforce_cpmd
