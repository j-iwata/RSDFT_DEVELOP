PROGRAM Real_Space_Solid

  use global_variables
  use parameters_module
  use func2gp_module
  use band_module
  use band_sseig_module
  use hamiltonian_matrix_module
  use rtddft_mol_module
  use omp_variables, only: init_omp
  use test_rtsol_module

  use ps_nloc_initiate_module
  use ps_getDij_module

  use WFtest
  use io_tools_module, only: init_io_tools, IOTools_readIntegerKeyword
  use lattice_module

  implicit none
  integer,parameter :: unit_input_parameters = 1
  integer,parameter :: unit_atomic_coordinates = 970
  real(8) :: ct0,ct1,et0,et1,exc_tmp,eh_tmp,eion_tmp,tmp,shift_factor
  integer :: i,n,k,s,iter,m,ierr,i1,i2,i3,m1,m2,m3,j,mm1,mm2,mm3,info
  real(8),allocatable :: force(:,:),forcet(:,:),vtmp(:)
  real(8) :: totalScfTime
  logical :: isDebug=.false.

  type(lattice) :: aa_obj, bb_obj

! --- start MPI ---

  call start_mpi_parallel

! --- global time counter start ---

  call global_watch(.false.)

! --- info ---

  call open_info(myrank)

! --- DISP_SWITCH ---

! call setDispSwitch(myrank,nprocs)
! DISP_SWITCH = .true.
  DISP_SWITCH = (myrank==0)
  disp_switch_parallel = (myrank==0)

  call check_disp_switch( DISP_SWITCH, 1 )
     
! --- input parameters ---

  call init_io_tools( myrank, unit_input_parameters )

  call read_parameters

! ---  Type of System ( RS-SOL or RS-MOL ) ---

  call IOTools_readIntegerKeyword( "SYSTYPE", Systype )

! --- atomic coordinates & R-space Lattice ---

  call check_format_atom( myrank, unit_atomic_coordinates, aa_obj )

  call write_border( 80, " main( aa & bb )(start)" )

  call init_aa( aa_obj )

  if ( SYStype == 0 ) then
     call convert_to_aa_coordinates_atom( aa_obj, aa_atom )
  else if ( SYStype == 1 ) then
     call convert_to_xyz_coordinates_atom( aa_obj, aa_atom )
  end if

  call backup_aa_lattice( aa_obj )

  call write_border( 80, " main( aa & bb )(end)" )

! --- Reciprocal Lattice ---

  call construct_bb(aa)
  call get_reciprocal_lattice( aa_obj, bb_obj )

! --- R-space & G-space Grid ---

  call Init_Rgrid( SYStype, Md, unit=2 )

  call Init_Ggrid( Ngrid, bb, Hgrid, disp_switch )

! --- Test ( Egg-Box Effect ) ---

  if ( iswitch_test == 2 ) then
     aa_atom(1,1) = aa_atom(1,1) + Hgrid(1)*0.5d0/ax
     if ( disp_switch ) then
        write(*,*) "--- EGG BOX TEST !!! ---"
        write(*,*) aa_atom(1,1),aa_atom(1,1)*ax,Hgrid(1)*0.5d0
     end if
  end if

! --- Symmetry ---

  call init_symmetry( Ngrid,dV,aa,bb, Natom,ki_atom,aa_atom )

! --- Brillouin Zone sampling ---

  if ( isymmetry == 0 ) then
     call generate_bz( disp_switch )
  else
     call generate_bz_sym( nsym, rgb, disp_switch )
  end if

! --- initial set up for parallel computation ---

!  call test_bcast

  call init_scalapack( Nband )

  call init_parallel( Ngrid, Nband, Nbzsm, Nspin )

  call InitParallel_Rgrid

  call InitParallel_Ggrid( nprocs, myrank )

  call prep_symmetry( Igrid )

! --- initialization for thread-parallel computation ---

  call init_omp( Igrid(1,1),Igrid(2,1),Igrid(1,2),Igrid(2,2) &
                ,Igrid(1,3),Igrid(2,3),Igrid(1,0),Igrid(2,0) &
                ,SYStype, disp_switch )

!- FD boundary set -

  call init_bcset( Md, SYStype )

! --- kinetic energy oprator coefficients ---

  call init_kinetic( aa, bb, Nbzsm, kbb, Hgrid, Igrid, MB_d, DISP_SWITCH )

! --- ??? ---

  call set_array_bound

! --- Pseudopotential, initial density, and partial core correction ---

  call read_pseudopot( Nelement, myrank )

! ---

  call write_info_atom( Zps, file_ps )

!-------- init density 

  call count_electron

  call init_density(Nelectron,dV)

!----------------------- SOL sol -----

  if ( SYStype == 0 ) then

     call init_ps_local
     call init_ps_pcc

     call construct_strfac  !----- structure factor

     call construct_ps_local
     call construct_ps_pcc
     call construct_ps_initrho( rho )

     call destruct_strfac   !----- structure factor

     call ps_nloc_initiate( Gcut )

!----------------------- MOL mol -----

  else if ( SYStype == 1 ) then

     call init_ps_local_mol(Gcut)
     call init_ps_pcc_mol
     call init_ps_initrho_mol

     call Construct_RgridMol(Igrid)

     call construct_ps_local_mol
     call construct_ps_pcc_mol
     call construct_ps_initrho_mol
     call normalize_density

     call ps_nloc2_init(Gcut)
     call prep_ps_nloc2_mol

     call ConstructBoundary_RgridMol(Md,Igrid)

  end if

  call sym_rho( ML_0, ML_1, Nspin, MSP_0, MSP_1, rho )

! --- Initial setup for Hybrid XC functional ---

  call init_xc_hybrid( ML_0, ML_1, Nelectron, Nspin, Nband &
       , MMBZ, Nbzsm, MBZ_0, MBZ_1, MSP, MSP_0, MSP_1, MB_0, MB_1 &
       , kbb, bb, Va, SYStype, np_fkmb, disp_switch )

!-------------------- Hamiltonian Test

  if ( iswitch_test == 1 ) then
     call test_hpsi2( 10 )
     goto 900
  end if

!-------------------- BAND with SSEIG

  if ( iswitch_band == 2 ) then
     call init_localpot
     call read_localpot("vrho.dat1",myrank)
     call band_sseig(disp_switch)
     goto 900
  end if

!--------------------

! --- Ion-Ion ---

  call init_eion( SYStype, disp_switch )

! --- Initialization of subspace diagonalization ---

  call init_subspace_diag( Nband )

! --- Initial wave functions ---

  call read_wf( myrank, 1 )
  call init_wf( SYStype )

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

! --- Initial occupation ---

  call init_occ_electron(Nelectron,Ndspin,Nbzsm,weight_bz,occ)

  if ( DISP_SWITCH ) then
     write(*,'(a60," main")') repeat("-",60)
     write(*,*) "Natom    =",Natom
     write(*,*) "Nelement =",Nelement
     write(*,*) "Nband =",Nband
     write(*,*) "Nspin =",Nspin
     write(*,*) "Nbzsm =",Nbzsm,MBZ
     write(*,*) "Nelectron =",Nelectron
     write(*,*) "Next_electron =",Next_electron
     write(*,*) "Ndspin,Nfixed =",Ndspin,Nfixed
     write(*,*) "Zps   =",Zps(1:Nelement)
     write(*,*) "sum(occ)=",sum(occ)
     if ( Nspin == 2 ) then
        write(*,*) "sum(occ(up))  =",sum(occ(:,:,1))
        write(*,*) "sum(occ(down))=",sum(occ(:,:,Nspin))
     endif
     do n=max(1,nint(Nelectron/2)-20),min(nint(Nelectron/2)+80,Nband)
        do k=1,Nbzsm
           write(*,*) n,k,(occ(n,k,s),s=1,Nspin)
        end do
     end do
  end if

! --- Initial Potential ---

  call init_hartree( Igrid, Nspin, Md, SYStype )
  call calc_hartree( ML_0, ML_1, MSP, rho )

  call calc_xc

  call init_localpot

  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

!-------------------- Real-Time Test

  if ( iswitch_test == 3 ) then
     call test_rtsol
     goto 900
  end if

! --- Read previous w.f. , density , potentials ---

  call read_data(disp_switch)

  call getDij

! the following GS should be performed when MB1_tmp is smaller than Nband,
! otherwise not necessary

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

! ---

  if ( GetParam_IO(1) == 1 .or. GetParam_IO(1) >= 3 ) then
     call control_xc_hybrid(1)
     if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid
  end if

! ---
! The followings are just to get H and XC energies,
! but potentials are also recalculated with new rho.

  call calc_hartree(ML_0,ML_1,MSP,rho)
  call calc_xc
  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

! --- Init vdW ---

  call read_vdw_grimme(myrank,1)
  call init_vdw_grimme(XCtype,aa,Natom,nprocs,myrank,ki_atom,zn_atom)
  call calc_E_vdw_grimme( Natom, aa_atom )

! ---

  call getDij

  call calc_with_rhoIN_total_energy(disp_switch)
  call calc_total_energy(.true.,disp_switch,999)

! ---

  if ( Nsweep > 0 ) then
!     call init_sweep( 2, Nband, 1.d-7 )
     call calc_sweep( Nsweep, ierr, disp_switch )
     if ( ierr == -1 ) goto 900
  end if

! ---

  select case( iswitch_scf )
  case( 1 )
     call calc_scf( Diter_scf, ierr, disp_switch, feps )
     if ( ierr < 0 ) goto 900
  case( 2 )
     call calc_scf_chefsi( Diter_scf_chefsi, ierr, disp_switch )
     if ( ierr < 0 ) goto 900
  case( -1 )
     if ( nprocs == 1 ) then
        call construct_hamiltonian_matrix( Ngrid(0) )
     end if
     goto 900
  end select

! ---

!
! --- BAND ---
!
#ifndef _DRSDFT_
  if ( iswitch_band == 1 ) then
     call control_xc_hybrid(1)
     call band(nint(Nelectron*0.5d0),disp_switch)
  end if
#endif

!
! --- Force test, atomopt, CPMD ---
!

  select case( iswitch_opt )
  case( -1 )

     call test_force(SYStype)

  case( 1,2 ) ! --- atomopt ---

     call atomopt(iswitch_opt,disp_switch)

  case( 3 ) ! --- CPMD ---

#ifdef _DRSDFT_
     call bomd
#else
     write(*,*) "RS-CPMD is not available for COMPLEX16"
     write(*,*) "Please re-compile the program"
#endif

  end select

!
! --- TDDFT ---
!
  select case( iswitch_tddft )
  case( 0 )

  case( 1,2 )

     select case( SYStype )
     case( 1 )
        call init_rtddft_mol( 1, myrank )
        call rtddft_mol( iswitch_tddft )
        goto 900
     case default
        write(*,*) "real-time tddft is available only for rsmol"
        goto 900
     end select

  case default

     write(*,'(1x,"iswitch_tddft=",i2," is not available")') iswitch_tddft

  end select

! --- finalize ---

  if ( DISP_SWITCH ) write(*,*) "END_PROGRAM : MAIN" 

900 continue
  if ( DISP_SWITCH ) write(*,*) 'intentional end'
  call global_watch(disp_switch)
  call close_info
  call end_mpi_parallel

END PROGRAM Real_Space_Solid
