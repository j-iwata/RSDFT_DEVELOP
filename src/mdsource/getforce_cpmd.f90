!-----------------------------------------------------------------------
!     Evaluate Force for CPMD
!-----------------------------------------------------------------------
SUBROUTINE getforce_cpmd( ltime )

  use eion_module, only: calc_eion
  use atom_module, only: Natom,aa_atom,ki_atom,zn_atom,shift_aa_coordinates_atom
  use bb_module
  use parallel_module, only: myrank, np_band, myrank_b, comm_band
  use strfac_module
  use ps_local_module
  use ps_pcc_module, only: flag_pcc_0, construct_ps_pcc
  use pseudopot_module
  use ps_nloc2_module
  use ps_nloc3_module
  use ps_nloc_mr_module
  use localpot_module, only: Vloc
  use array_bound_module, only: MSP_0,MSP_1,MB_0,MB_1,ML_0,ML_1,MBZ_0,MBZ_1
  use density_module, only: rho, calc_density_2
  use xc_module
  use hartree_variables, only: Vh, E_hartree
  use hartree_module, only: calc_hartree
  use cpmd_variables, only: Force,Rion,disp_switch,AMU,pmass,ir_band_cpmd,id_band_cpmd
  use watch_module
  use force_module
  use wf_module, only: unk, occ
  use rsdft_mpi_module
  use construct_vion_vh_floc_module, only: construct_vion_vh_floc
  use nonlocal_module, only: calc_force_nonlocal
  use force_ewald_module, only: calc_force_ewald
  use ps_nloc2_variables, only: prep_backup_uVunk_ps_nloc2
  use fft_module, only: iswitch_fft

  implicit none

  logical,intent(IN) :: ltime
  integer :: i,s
  real(8) :: ctime_force(0:9),etime_force(0:9),c,ttmp(2),ttt(2,2)
  real(8),allocatable :: work2(:,:)

  c=1.0d0/(2.0d0*acos(-1.0d0))
  aa_atom = matmul(transpose(bb),Rion)*c
  call shift_aa_coordinates_atom( aa_atom )
  !Force   = 0.0d0

  if ( ltime ) call watch(ctime_force(0),etime_force(0))

  call calc_eion

  if ( ltime ) call watch(ctime_force(1),etime_force(1))

  select case( iswitch_fft )
  case( 'FFT0','FFTW' )
    call construct_strfac
    call construct_ps_local
    call construct_ps_pcc
    call destruct_strfac
  case( 'FFTE','FFTE1','FFTE2' )
    if ( flag_pcc_0 ) then
      call construct_strfac
      call construct_ps_pcc
      call destruct_strfac
    end if
  end select

  if ( ltime ) call watch(ctime_force(2),etime_force(2))

  select case(pselect)
  case(2)
    call prep_ps_nloc2
  case(3)
    call prep_ps_nloc3
  case(5)
    call prep_ps_nloc_mr
  end select

  if ( ltime ) call watch(ctime_force(3),etime_force(3))

  if ( np_band > 1 ) call wf_gather_sub( unk ) ! This lien should be removed for memory-band-parallel CPMD.

  call calc_density_2( unk(:,MB_0:MB_1,:,:), occ(MB_0:MB_1,:,:) )

  if ( ltime ) call watch(ctime_force(4),etime_force(4))

  select case( iswitch_fft )
  case( 'FFT0', 'FFTW' )
    call calc_hartree( rho )
  case( 'FFTE','FFTE1','FFTE2' )
    call construct_vion_vh_floc( rho, Vion, Vh, Force, E_hartree )
  end select

  if ( ltime ) call watch(ctime_force(5),etime_force(5))

  call calc_xc

  if ( ltime ) call watch(ctime_force(6),etime_force(6))

  do s=MSP_0,MSP_1
    Vloc(:,s)=Vh(:)+Vxc(:,s)+Vion(:)
  enddo

  if ( ltime ) call watch(ctime_force(7),etime_force(7))

  call prep_backup_uVunk_ps_nloc2( MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1 )

  select case( iswitch_fft )
  case( 'FFTE','FFTE1','FFTE2' )
    allocate( work2(3,Natom) ); work2=0.0d0
    !call watchb( ttmp, barrier='on' ); ttt=0.0d0
    call calc_force_nonlocal( work2 )
    !call watchb( ttmp, ttt(:,1), barrier='on' )
    Force=Force+work2
    !call watchb( ttmp, barrier='on' )
    call calc_force_ewald( Natom, work2 )
    !call watchb( ttmp, ttt(:,2), barrier='on' )
    Force=Force+work2
    deallocate( work2 )
  case( 'FFT0', 'FFTW' )
    call calc_force( Natom, Force )
  end select

  if ( ltime ) call watch(ctime_force(8),etime_force(8))

  do i=1,Natom
    Force(:,i)=Force(:,i)/(pmass(zn_atom(ki_atom(i)))*AMU)
  enddo

  if ( ltime ) call watch(ctime_force(9),etime_force(9))

  if ( ltime .and. myrank==0 ) then
    write(17,'(9f10.5)') (etime_force(i+1)-etime_force(i),i=0,8)
    !write(*,'("fnloc,fewld",4f10.5)') ttt
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
