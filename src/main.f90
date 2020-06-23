PROGRAM Real_Space_DFT

  use parallel_module
  use global_variables, only: iswitch_test,iswitch_scf,iswitch_tddft,iswitch_band,iswitch_opt,iswitch_dos,iswitch_latopt
  use kinetic_variables, only: SYStype, Md, kin_select
  use grid_module, only: grid_info
  use rgrid_module
  use ggrid_module
  use atom_module
  use pseudopot_module
  use bb_module
  use electron_module
  use bz_module
  use density_module
  use ps_local_module, only: Vion
  use array_bound_module, only: ML_0,ML_1,MSP_0,MSP_1,MBZ_0,MBZ_1,MBZ,MSP,MB_0,MB_1,set_array_bound
  use wf_module, only: occ, init_wf ,esp,unk
  use xc_module
  use localpot_module
  use hartree_variables
  use xc_hybrid_module
  use io_module
  use atomopt_module
  use scf_chefsi_module
  use scf_module
  use total_energy_module
  use sweep_module

  use watch_module
  use info_module
  use force_module
  use symmetry_module
  use scalapack_module
  use bc_module
  use kinetic_module
  use rgrid_mol_module
  use test_hpsi2_module
  use eion_module
  use subspace_diag_module
  use gram_schmidt_module
  use init_occ_electron_module
  use hartree_module
  use test_force_module

  use parameters_module, only: read_parameters
  use func2gp_module
  use band_module
  use dos_module, only: calc_dos
  use hamiltonian_matrix_module
  use rtddft_mol_module
  use omp_variables, only: init_omp
  use test_rtsol_module
  use ps_getDij_module, only: getDij
  use ps_init_module, only: ps_init
  use io_tools_module, only: init_io_tools, IOTools_readIntegerKeyword, IOTools_readStringKeyword, IOTools_findKeyword
  use lattice_module
  use ffte_sub_module, only: init_ffte_sub
  use fftw_module
  use vdw_grimme_module
  use efield_module
  use stress_module, only: test_stress ! MIZUHO-IR for cellopt
  use linear_response_module
  use kinetic_sym_ini_module
  use kinetic_allgatherv_module, only: init_kinetic_allgatherv
  use noncollinear_module, only: flag_noncollinear, io_read_noncollinear &
                                ,init_noncollinear, calc_xc_noncollinear
  use init_occ_electron_ncol_module
  use rtddft_sol_module
  use aa_module, only: init_aa
  use rsdft_mpi_module, only: init_rsdft_mpi
  use allel_module, only: init_allel
!  use bcast_module, only: test_bcast
!  use rsdft_sendrecv_module, only: test_sendrecv
!  use vector_tools_module, only: set_vinfo_vector_tools
!  use sl_tools_module, only: backup_param_sl_tools

  implicit none
  integer,parameter :: unit_input_parameters = 1
  integer,parameter :: unit_atomic_coordinates = 970
  real(8) :: ct0,ct1,et0,et1,exc_tmp,eh_tmp,eion_tmp,tmp,shift_factor
  integer :: i,n,k,s,iter,m,ierr,i1,i2,i3,m1,m2,m3,j,mm1,mm2,mm3,info
  real(8),allocatable :: force(:,:),forcet(:,:),vtmp(:)
  type(lattice) :: aa_obj, bb_obj
  type(grid_info) :: rgrid
  logical,parameter :: recalc_esp=.true.
  logical :: flag_read_ncol=.false., DISP_SWITCH
  logical :: skip_read_data=.false.
  real(8) :: Etot, Ehwf
  integer :: info_level=1
  character(32) :: lattice_index
  character(20) :: systype_in="SOL"
  integer :: nloop,itmp(3)

! --- start MPI ---

  call start_mpi_parallel

  call init_rsdft_mpi

! --- global time counter start ---

  call global_watch(.false.)

! --- info ---

  call open_info(myrank)

! --- init_io_tools ---

  if ( myrank == 0 ) then
     open(unit_input_parameters  ,file="fort.1" ,status="old")
     open(unit_atomic_coordinates,file="fort.970",status="old")
  end if

  call init_io_tools( myrank, unit_input_parameters )

! --- STDOUT and LOG control (ext0/write_info.f90) ---

  DISP_SWITCH = (myrank==0)
  disp_switch_parallel = (myrank==0)

  call check_disp_switch( DISP_SWITCH, 1 )
  call check_log_switch( DISP_SWITCH, 1 )

  call IOTools_readIntegerKeyword( "INFOLEVEL", info_level )

  call check_disp_length( info_level, 1 )
     
! ---  Type of System ( RS-SOL or RS-MOL ) ---

  call IOTools_readStringKeyword( "SYSTYPE", systype_in )
  call convert_to_capital( systype_in )
  SYStype=0
  select case( systype_in(1:3) )
  case( "SOL" )
     SYStype=0
  case( "MOL" )
     SYStype=1
  end select

! --- input parameters ---

  call read_parameters

  call read_atom( myrank, unit_atomic_coordinates, aa_obj )

  call read_lattice( aa_obj )

  call construct_lattice( aa_obj )
  call backup_aa_lattice( aa_obj )

! --- coordinate transformation of atomic positions ---

  if ( SYStype == 0 ) then
     call convert_to_aa_coordinates_atom( aa_obj, aa_atom )
     if ( myrank == 0 ) then
        call write_coordinates_atom( 96, 3, "atomic_coordinates_aa_ini" )
     end if
  else if ( SYStype == 1 ) then
     call convert_to_xyz_coordinates_atom( aa_obj, aa_atom )
     if ( myrank == 0 ) then
        call write_coordinates_atom( 96, 3, "atomic_coordinates_xyz_ini" )
     end if
  end if

! --- Real-Space Grid (MOL) ( lattice is defiend by MOLGRID ) ---

  if ( SYStype /= 0 ) then
     call Init_Rgrid( SYStype, Md, aa_obj )
     call construct_lattice( aa_obj )
  end if

! --- Lattice ---

  call backup_aa_lattice( aa_obj )

  call init_aa( aa_obj )

  call check_lattice( aa_obj%LatticeVector, lattice_index )
  if ( disp_switch ) write(*,*) "lattice_index: ",lattice_index

! --- Real-Space Grid (SOL) ( lattice is used to define the grid ) ---

  if ( SYStype == 0 ) call Init_Rgrid( SYStype, Md, aa_obj )

! --- Reciprocal Lattice ---

  call construct_bb( aa_obj%LatticeVector )

! --- G-space Grid ---

  call Init_Ggrid( Ngrid, bb )

! --- Pseudopotential ---

  call read_pseudopot( Nelement, myrank )

  call init_allel( Zps, abs(aa_obj%Volume) )

! --- info atoms ---

  call write_info_atom( Zps, file_ps )

! --- count the total # of electrons & check # of bands ---

  call count_electron

  call check_Nband_electron

! --- init_force ---

  call init_force( SYStype )

! --- Test ( Egg-Box Effect ) ---

  if ( iswitch_test == 10 ) then
     aa_atom(1,:) = aa_atom(1,:) + Hgrid(1)*0.5d0/aa_obj%LatticeConstant
     if ( disp_switch ) then
        write(*,*) "--- EGG BOX TEST !!! ---"
        do i=1,size(aa_atom,2)
           write(*,*) aa_atom(1,i),aa_atom(1,i)*aa_obj%LatticeConstant,Hgrid(1)*0.5d0
        end do
     end if
  end if

! --- Symmetry ---

  call init_symmetry( Ngrid,dV,aa_obj%LatticeVector,bb, Natom,ki_atom,aa_atom )

! --- Brillouin Zone sampling ---

  call generate_bz

  if ( myrank == 0 ) call write_info_bz( bb )

! --- initial set up for parallel computation ---

!  call test_bcast
!  call test_sendrecv( node_partition )
!  goto 900

  call init_scalapack( Nband )
!  call backup_param_sl_tools( NPROW, NPCOL, MBSIZE, NBSIZE )

  call init_parallel( Ngrid, Nband, Nbzsm, Nspin )

  call InitParallel_Rgrid

  call InitParallel_Ggrid( nprocs, myrank )

  call prep_symmetry( Igrid )

! --- initialization for thread-parallel computation ---

  call init_omp( Igrid(1,1),Igrid(2,1),Igrid(1,2),Igrid(2,2) &
                ,Igrid(1,3),Igrid(2,3),Igrid(1,0),Igrid(2,0) &
                ,SYStype, disp_switch )

! --- grid info ( real space ) ---

  call set_grid_info_rgrid( rgrid, aa_obj%LatticeVector, Md )

! --- Initialization for FFT ---

  if ( SYStype == 0 ) then

     itmp(:)=(/Igrid(1,1),Igrid(1,2),Igrid(1,3)/)
     call init_ffte_sub(itmp,Ngrid(1:3),node_partition(1:3),comm_grid)

     call init_fftw( Ngrid(1:3), node_partition(1:3), comm_grid, myrank_g )

  end if

!- FD boundary set -

  call IOTools_readIntegerKeyword( "MD", Md )

  call init_bcset( Md, SYStype )

! --- kinetic energy oprator coefficients ---

  call read_kinetic ! Md & kin_select

  call init_kinetic( aa_obj%LatticeVector, bb, Nbzsm, kbb, Hgrid, Igrid, MB_d )

  if ( kin_select == 2 ) then
     call init_kinetic_sym( lattice_index, aa_obj%LatticeVector, ierr )
     if ( ierr /= 0 ) kin_select=0
  else if ( kin_select == 3 ) then
     call init_kinetic_allgatherv( Igrid, comm_grid )
  end if

! --- ??? ---

  call set_array_bound

  if ( SYStype == 1 ) then ! MOL mol
     call Construct_RgridMol(Igrid)
     call ConstructBoundary_RgridMol(Md,Igrid)
  end if

! --- vector info ---

!  call set_vinfo_vector_tools( 1, dV   , comm_grid, np_grid, myrank_g, ir_grid, id_grid )
!  call set_vinfo_vector_tools( 2, 1.0d0, comm_band, np_band, myrank_b, ir_band, id_band )

! --- init density & potentials ---

  call init_density( Nelectron, dV )

  call ps_init( SYStype, Gcut, rho )

  call normalize_density( rho )

! --- External electric field (included in the local ionic potential) ---

  call sawtooth_efield( Vion )

! --- symmetrize density ---

  call sym_rho( ML_0, ML_1, Nspin, MSP_0, MSP_1, rho )

! --- init noncollinear ---

  flag_noncollinear = flag_so
  if ( flag_noncollinear ) then
     if ( Nspin /= 2 .or. np_spin == 2 ) then
        write(*,*) "Nspin, np_spin=",Nspin,np_spin
        write(*,*) "Nspin==2 and np_spin==1 is only available"
        goto 900
     end if
  end if

!-------------------- Hamiltonian Test

  if ( iswitch_test == 1 ) then
     nloop=10
     call IOTools_readIntegerKeyword( "NLOOP", nloop )
     call test_hpsi2( nloop )
     goto 900
  end if

!--------------------

! --- Ion-Ion ---

  call init_eion( SYStype, disp_switch )

! --- Initialization of subspace diagonalization ---

  call init_subspace_diag( Nband )

! --- Initial wave functions ---

  call init_wf( SYStype )

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

! --- Initial occupation ---

  if ( flag_noncollinear ) then
     call init_occ_electron_ncol(Nelectron,Ndspin,Nbzsm,weight_bz,occ)
  else
     call init_occ_electron(Nelectron,Ndspin,Nbzsm,weight_bz,occ)
  end if

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
     write(*,*) "sum(occ)=",sum(occ),count(occ(:,1,1)/=0.0d0)
     if ( Nspin == 2 ) then
        write(*,*) "sum(occ(up))  =",sum(occ(:,:,1)),count(occ(:,1,1)/=0.0d0)
        write(*,*) "sum(occ(down))=",sum(occ(:,:,Nspin)),count(occ(:,1,Nspin)/=0.0d0)
     endif
     !do n=max(1,nint(Nelectron/2)-10),min(nint(Nelectron/2)+10,Nband)
     !   do k=1,Nbzsm
     !      write(*,*) n,k,(occ(n,k,s),s=1,Nspin)
     !   end do
     !end do
  end if

! --- Initial setup for Hybrid XC functional ---

  call init_xc_hybrid( ML_0, ML_1, Nelectron, Nspin, Nband &
       , MMBZ, Nbzsm, MBZ_0, MBZ_1, MSP, MSP_0, MSP_1, MB_0, MB_1 &
       , kbb, bb, aa_obj%Volume, SYStype, np_fkmb, disp_switch )

! --- Initial Potential ---

  call init_hartree( Igrid, Nspin, Md, SYStype )
  call calc_hartree( ML_0, ML_1, MSP, rho )

  call calc_xc

  if ( flag_noncollinear ) then
     call init_noncollinear( rho, Vxc )
     Vxc=0.0d0
  end if

  call init_localpot( ML_0,ML_1, MSP_0,MSP_1 )

  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

!-------------------- Force Test

  if ( iswitch_test == 2 ) then
     call test_force( SYStype, .false. )
     goto 900
  end if

!-------------------- Real-Time Test

  if ( iswitch_test == 3 ) then
     call test_rtsol
     goto 900
  end if

! --- Read previous w.f. , density , potentials ---

  call IOTools_findKeyword( "SKIP_READ_DATA", skip_read_data, flag_bcast=.true. )

  if ( .not.skip_read_data ) call read_data(disp_switch)

  call io_read_noncollinear( myrank, flag_read_ncol )

  call getDij

! the following GS should be performed when MB1_tmp is smaller than Nband,
! otherwise not necessary

  if ( flag_read_ncol ) then
  else

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
     call gram_schmidt(1,Nband,k,s)
  end do
  end do

  end if

! ---

  if ( GetParam_IO(1) == 1 .or. GetParam_IO(1) >= 3 ) then
     call control_xc_hybrid(1)
     if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid
  end if

! ---
! The followings are just to get H and XC energies,
! but potentials are also recalculated with new rho.

  if ( flag_read_ncol ) then

     call calc_xc_noncollinear( unk, occ(:,:,1), rho, Vxc )
     call calc_hartree( ML_0,ML_1,MSP,rho )

  else

     call calc_hartree(ML_0,ML_1,MSP,rho)
     call calc_xc
     do s=MSP_0,MSP_1
        Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
     end do

  end if

! ---

  call getDij

! --- init_vdW_Grimme ---

  call init_vdw_grimme( aa_obj%LatticeVector, ki_atom, zn_atom )
  call calc_E_vdw_grimme( aa_atom )

! --- total energy ---

  call calc_with_rhoIN_total_energy( Ehwf, flag_ncol=flag_noncollinear )
  call calc_total_energy( recalc_esp, Etot, flag_ncol=flag_noncollinear )

! ---

  call write_border( 0, 'End Initialization' )
  call global_watch( disp_switch, indx='INIT')
  call write_border( 0, '' )



! ---

  if ( iswitch_dos == 1 ) then

     call control_xc_hybrid(1)
     call calc_dos( ierr )
     goto 900

  else

     call calc_sweep( ierr, flag_ncol_in=flag_noncollinear )
     if ( ierr < 0 ) goto 900

  end if

! ---

  select case( iswitch_scf )
  case( 1 )
     call calc_scf( ierr, tol_force_in=feps, Etot_out=Etot )
     if ( ierr < 0 ) goto 900
     call calc_total_energy( recalc_esp,Etot,unit_in=6,flag_ncol=flag_noncollinear )
  case( 2 )
     call calc_scf_chefsi( Diter_scf_chefsi, ierr, disp_switch )
     if ( ierr < 0 ) goto 900
     call calc_total_energy( recalc_esp,Etot,unit_in=6,flag_ncol=flag_noncollinear )
  case( -1 )
     if ( nprocs == 1 ) then
        call construct_hamiltonian_matrix( Ngrid(0) )
     end if
     goto 900
  end select

!
! --- Force test, atomopt, CPMD ---
!
  if( iswitch_opt == -1 ) then
     call test_force( SYStype, DISP_SWITCH )
  end if

  if( iswitch_latopt == -1 ) then
     call test_stress(SYStype)
  end if

  if( iswitch_opt >= 1 .or. iswitch_latopt >= 1 ) then

     if ( iswitch_opt /= 3 ) then
        call atomopt(iswitch_opt,iswitch_latopt)
        call calc_total_energy( recalc_esp, Etot, unit_in=6 )
     else
        if ( SYStype == 0 ) then
           !if ( flag_noncollinear ) then
           !   call init_occ_electron_ncol(Nelectron,Ndspin,Nbzsm,weight_bz,occ)
           !else
           !   call init_occ_electron(Nelectron,Ndspin,Nbzsm,weight_bz,occ)
           !end if
           call bomd
        else
           write(*,*) "MD for SYStype/=0 is not avaiable"
           goto 900
        end if
     end if

  end if

!
! --- BAND ---
!

  if ( iswitch_band > 0 ) then
     call control_xc_hybrid(1)
     call band(nint(Nelectron*0.5d0),disp_switch,iswitch_band)
  end if

!
! --- TDDFT ---
!
  select case( iswitch_tddft )
  case( 0 )

  case( 1,2 )

     select case( SYStype )
     case default
        call rtddft_sol( iswitch_tddft )
     case( 1 )
        call init_rtddft_mol( 1, myrank )
        call rtddft_mol( iswitch_tddft )
     end select

  case( 3 )

     call calc_dielectric_constant

  case default

     write(*,'(1x,"iswitch_tddft=",i2," is not available")') iswitch_tddft

  end select

! --- finalize ---

  if ( DISP_SWITCH ) write(*,*) "END_PROGRAM : MAIN" 

900 continue

  call global_watch(disp_switch)
  if ( SYStype == 0 ) call finalize_fftw
  call close_info
  call end_mpi_parallel

END PROGRAM Real_Space_DFT
