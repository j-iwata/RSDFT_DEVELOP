PROGRAM Real_Space_Solid

  use global_variables
  use parameters_module
  use omp_variables

  use esm_rgrid_module
  use esm_rshell_module
  use esm_cylindrical_test
  use ps_local_rs_module
  use esm_genpot_module
  use esm_kinetic_module

  use func2gp_module

  use rgrid_mol_module
  use ps_local_mol_module
  use eion_mol_module
  use ps_pcc_mol_module
  use ps_initrho_mol_module
  use ps_nloc2_mol_module
  use bc_mol_module
  use kinetic_mol_module
#ifndef _DRSDFT_
  use band_module
#endif

  use ps_gth_module
  use ps_nloc_mr_module

  use bcast_module

  use localpot2_variables
  use localpot2_ion_module
  use localpot2_density_module
  use localpot2_vh_module
  use localpot2_xc_module
  use localpot2_module
  use localpot2_te_module

  use ps_nloc3_module

  use test_hpsi2_module
  use test_force_module

  use info_module

  use init_occ_electron_module

  use esp_calc_module

  implicit none

  real(8) :: ct0,ct1,et0,et1,exc_tmp,eh_tmp,eion_tmp,tmp,shift_factor
  integer :: i,n,k,s,iter,m,ierr,i1,i2,i3,m1,m2,m3,j,mm1,mm2,mm3
  real(8),allocatable :: esp0(:,:,:),force(:,:),forcet(:,:),vtmp(:)
  real(8),allocatable :: work(:,:,:)
  logical :: flag_conv=.false.
  logical :: flag_exit=.false.
  logical :: flag_end =.false.
  logical :: flag_scf =.false.

! --- start MPI ---

  call start_mpi_parallel

! --- global time counter start ---

  call global_watch(flag_end)

! --- info ---

  call open_info(myrank)

! --- DISP_SWITCH ---

! DISP_SWITCH = .true.
  DISP_SWITCH = (myrank==0)
  disp_switch_parallel = (myrank==0)

  if ( disp_switch ) then
#ifdef _DRSDFT_
     write(*,*) "DRSDFT(REAL8)"
#else
     write(*,*) "RSDFT(COMPLEX16)"
#endif
#ifdef _SPLINE_
     write(*,*) "SPLINE(ps_nloc2_module)"
#endif
#ifdef _FFTE_
     write(*,*) "FFTE(ps_local,hartree)"
#endif
  end if
     
! --- input parameters ---

  if (DISP_SWITCH) write(*,'(a60," read_param")') repeat("-",60)

  call read_parameters
!  call read_oldformat_parameters

! --- initial preparetaion ---

  call construct_aa
  call construct_bb(aa)

! --- RSMOL ---

  if ( SYStype == 1 ) then
     call read_rgrid_mol(myrank,2)
     call init_rgrid_mol( Ngrid(1),Hgrid(1),aa,bb,disp_switch )
  end if

! --- R-grid and G-grid ---

  call construct_rgrid(aa)
  call get_cutoff_ggrid
  call construct_NMGL_ggrid

  if ( DISP_SWITCH ) then
     write(*,*) "Gcut,Ecut=",Gcut,Ecut
     write(*,*) "NGgrid=",NGgrid
     write(*,*) "NMGL=",NMGL
     write(*,'(1x,"Hgrid(1:3)=",3f20.10)') Hgrid(1:3)
  end if

! --- Brillouin Zone sampling ---

  call generate_bz(disp_switch)

! --- initial set up for parallel computation ---

  call test_bcast

  call prep_0_scalapack(Nband,disp_switch)

  call init_parallel(disp_switch,Ngrid(1),Nband,Nbzsm,Nspin)

  call parallel_rgrid(node_partition(1:3),myrank_g)
  call parallel_ggrid(nprocs,myrank)

  call init_bcset(Md)

! --- parallel computation for RSMOL ---

  if ( SYStype == 1 ) then

     select case(iswitch_eqdiv)
     case default
        call mesh_div_1(np_grid,node_partition,Ngrid(1),pinfo_grid,disp_switch)
     case(2)
        call mesh_div_2
     end select

     id_grid(:)=pinfo_grid(7,:)
     ir_grid(:)=pinfo_grid(8,:)
     idisp(myrank)=id_grid(myrank_g)
     ircnt(myrank)=ir_grid(myrank_g)
     call mpi_allgather(idisp(myrank),1,mpi_integer,idisp,1 &
                       ,mpi_integer,mpi_comm_world,ierr)
     call mpi_allgather(ircnt(myrank),1,mpi_integer,ircnt,1 &
                       ,mpi_integer,mpi_comm_world,ierr)
     Ngrid(0)=sum(ir_grid)
     Igrid(1,0)=pinfo_grid(7,myrank_g)+1
     Igrid(2,0)=pinfo_grid(7,myrank_g)+pinfo_grid(8,myrank_g)
     Igrid(1,1)=pinfo_grid(1,myrank_g)+1
     Igrid(2,1)=pinfo_grid(1,myrank_g)+pinfo_grid(2,myrank)
     Igrid(1,2)=pinfo_grid(3,myrank_g)+1
     Igrid(2,2)=pinfo_grid(3,myrank_g)+pinfo_grid(4,myrank)
     Igrid(1,3)=pinfo_grid(5,myrank_g)+1
     Igrid(2,3)=pinfo_grid(5,myrank_g)+pinfo_grid(6,myrank)

     call init_bcset_mol(Md,Ngrid(1),np_grid,myrank_g,comm_grid,pinfo_grid)

  end if

! --- ESM esm ---

  if ( SYStype == 3 ) then

     call read_esm_rgrid(myrank,1)
     call prep_esm_rgrid(Md)
     call construct0_esm_rgrid

     id_grid(:)=0
     ir_grid(:)=0
     id_grid(myrank_g) = ML0_ESM-1
     ir_grid(myrank_g) = ML1_ESM-ML0_ESM+1
     call mpi_allgather(id_grid(myrank_g),1,mpi_integer,id_grid,1,mpi_integer,comm_grid,ierr)
     call mpi_allgather(ir_grid(myrank_g),1,mpi_integer,ir_grid,1,mpi_integer,comm_grid,ierr)
     if ( myrank == 0 ) then
        do i=0,np_grid-1
           write(*,*) i,id_grid(i),id_grid(i)+ir_grid(i)-1,ir_grid(i)
        end do
     end if
     idisp(:)=0
     ircnt(:)=0
     idisp(myrank)=id_grid(myrank_g)
     ircnt(myrank)=ir_grid(myrank_g)
     call mpi_allgather(idisp(myrank),1,mpi_integer,idisp,1,mpi_integer,comm_grid,ierr)
     call mpi_allgather(ircnt(myrank),1,mpi_integer,ircnt,1,mpi_integer,comm_grid,ierr)
     if ( myrank == 0 ) then
        do i=0,nprocs-1
           write(*,*) i,idisp(i)+1,idisp(i)+ircnt(i),ircnt(i)
        end do
     end if
     Ngrid(0) = ML_ESM
     Igrid(1,0) = id_grid(myrank_g)+1
     Igrid(2,0) = id_grid(myrank_g)+ir_grid(myrank_g)
     call flush(6)
     call mpi_barrier(mpi_comm_world,ierr)
     write(*,'(1x,i6,2x,8i8)') myrank,Igrid(:,:)
     call flush(6)

  end if

! --- array bounds ---

  call set_array_bound

! --- OpenMP parallel ---

  call init_omp(Igrid(1,1),Igrid(2,1),Igrid(1,2),Igrid(2,2) &
               ,Igrid(1,3),Igrid(2,3),ML_0,ML_1,disp_switch)

! --- kinetic energy oprator coefficients ---

  call get_coef_kinetic(aa,bb,MBZ,kbb,DISP_SWITCH,SYStype_in=SYStype)
  if ( SYStype == 1 ) call get_coef_kinetic_mol(Md)
  if ( SYStype == 3 ) call get_coef_esm_kinetic(aa,bb,MBZ,kbb,Md)

! --- Pseudopotential, initial density, and partial core correction ---

  select case(pselect)
  case default
     call read_pseudopot(myrank)
  case(4,5)
     call read_ps_gth(myrank)
  end select

  call count_electron

!-------- init density 

  call init_density(Nelectron,dV)

!----------------------- SOL sol -----

  if ( SYStype == 0 ) then

     call init_ps_local
     call init_ps_pcc
     call init_ps_initrho
     call watcht(disp_switch,"strf",0)

     call construct_strfac !----- structure factor
     call watcht(disp_switch,"strf",1)

#ifndef _FFTE_
     call construct_ps_local
#else
     call construct_ps_local_ffte
#endif
     call watcht(disp_switch,"loc&pcc",1)
     call watcht(disp_switch,"loc",1)

     if ( pselect /= 4 .and. pselect /= 5 ) then
        call construct_ps_pcc
        call watcht(disp_switch,"pcc",1)
        call construct_ps_initrho
        call normalize_density
     end if

     call destruct_strfac !----- structure factor

if (myrank==0) write(*,*) "start initPS"
     select case( pselect )
     case( 2 )
        call ps_nloc2_init(Gcut)
        call prep_ps_nloc2
     case( 3 )
        call init_ps_nloc3
        call prep_ps_nloc3
     case( 5 )
        call prep_ps_nloc_mr
     end select

!----------------------- ESM esm -----

  else if ( SYStype == 3 ) then

     call ps_nloc2_init(Gcut)
     call prep_ps_nloc2_esm

     call init_ps_local_rs
     call read_rshell(myrank,1)

     if ( allocated(Vion) ) deallocate(Vion)
     if ( allocated(rho) ) deallocate(rho)
     allocate( Vion(ML0_ESM:ML1_ESM)      ) ; Vion=0.d0
     allocate( rho(ML0_ESM:ML1_ESM,Nspin) ) ; rho=0.d0

     call construct_ps_local_rs(Vion)

     call construct_ps_initrho_rs(ML0_ESM,ML1_ESM,Nspin,rho)
     call normalize_density
!     c0=sum(rho)
!     call mpi_allreduce(c0,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
!     rho=rho-c/ML_ESM+Nelectron/(ML_ESM*dV)
     write(*,'(1x,"sum(rho)*dV",3f15.10)') sum(rho)*dV,minval(rho),maxval(rho)

     call flush(6)

!----------------------- MOL mol -----

  else if ( SYStype == 1 ) then

     call init_ps_local_mol(Gcut)
     call init_ps_pcc_mol
     call init_ps_initrho_mol

     call construct_rgrid_mol
     call construct_ps_local_mol
     call construct_ps_pcc_mol
     call construct_ps_initrho_mol
     call normalize_density

     call ps_nloc2_init(Gcut)
     call prep_ps_nloc2_mol

     call construct_boundary_rgrid_mol(Md)

  end if

!-------------------- Hamiltonian Test

  if ( iswitch_test == 1 ) then
     call test_hpsi2( 10 )
     goto 900
  end if

!-------------------- BAND with SSEIG
#ifndef _DRSDFT_
  if ( iswitch_band == 2 ) then
     call init_localpot
     call read_localpot("vrho.dat1",myrank)
     call read_band(myrank,1)
     call band_sseig(disp_switch)
     goto 900
  end if
#endif
!--------------------


! --- ESM esm ---

  if ( SYStype == 3 ) then

     call construct_ps_density_longloc
     call read_esm_genpot(myrank,1)
     allocate( vtmp(ML0_ESM:ML1_ESM) )
     vtmp=0.d0
     call esm_genpot(vtmp)
     Vion(:) = Vion(:) + vtmp(:)
     deallocate( vtmp )

  end if

! --- Ewald sum ---

  select case(SYStype)
  case default

     call watcht(disp_switch,"",0)
     call test_ewald(Eewald,disp_switch)
     call watcht(disp_switch,"test_ewald",1)
     call calc_ewald(Eewald,disp_switch)
     call watcht(disp_switch,"calc_ewald",1)

  case( 1 )

     call watcht(disp_switch,"",0)
     call calc_eion_mol(Eewald)
     call watcht(disp_switch,"eion_mol",1)
     if ( disp_switch ) write(*,*) "Ewld(MOL)=",Eewald

  end select 

! --- preparing for subspace diagonalization ---

  call prep_subspace_diag(Nband,disp_switch)

! --- Initial wave functions ---

  call init_wf

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
!     call gram_schmidt_m(1,Nband,k,s)
     call gram_schmidt_t(1,Nband,k,s)
  end do
  end do
!  call test_on_wf(dV,myrank==0)

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

  call calc_hartree(ML_0,ML_1,MSP,rho,SYStype)

  call calc_xc

  call init_localpot

  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do

! --- Read previous w.f. , density , potentials ---

  call watcht(disp_switch,"read_data",0)
  call read_data(disp_switch)
  call watcht(disp_switch,"read_data",1)

!--------------------------------------------- LOCALLPOT2
!
  if ( DISP_SWITCH ) write(*,'(a60," LOCALPOT2")') repeat("-",60)
  call read_localpot2(1,myrank)
  if ( flag_localpot2 ) then
     call prep_parallel_localpot2
     call init_localpot2(Ngrid(1),Ngrid(2),Ngrid(3))
     m1=Ngrid_dense(1)
     m2=Ngrid_dense(2)
     m3=Ngrid_dense(3)
     mm1=Igrid_dense(2,1)-Igrid_dense(1,1)+1
     mm2=Igrid_dense(2,2)-Igrid_dense(1,2)+1
     mm3=Igrid_dense(2,3)-Igrid_dense(1,3)+1
     call test_localpot2
     call localpot2_ion(mm1,mm2,mm3,Nelement,Ecut,vion_nl)
     call localpot2_density(mm1,mm2,mm3,rho_nl)
     call localpot2_vh(mm1,mm2,mm3,Ecut,rho_nl,vh_nl,E_hartree)
     call localpot2_xc(mm1,mm2,mm3,rho_nl,vxc_nl,Exc)
     vloc_dense(:,:,:) = vion_nl(:,:,:)+vh_nl(:,:,:)+vxc_nl(:,:,:)
     vloc_dense_old(:,:,:) = vloc_dense(:,:,:)
     call test2_localpot2(mm1,mm2,mm3,vloc_dense)
  end if
!
!---------------------------------------------
! ---

  call calc_hartree(ML_0,ML_1,MSP,rho,SYStype)
  call calc_xc

! ---

  call calc_with_rhoIN_total_energy(disp_switch)
  call calc_total_energy(.true.,disp_switch)

  if ( mod(imix,2) == 0 ) then
     call init_mixing(ML_1-ML_0+1,MSP_1-MSP_0+1, rho(ML_0,MSP_0))
  else
     call init_mixing(ML_1-ML_0+1,MSP_1-MSP_0+1,Vloc(ML_0,MSP_0))
  end if

  allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.d0

  flag_exit = .false.
  flag_scf  = .false.

  do iter=1,Diter

     if ( disp_switch ) write(*,'(a40," iter=",i4)') repeat("-",40),iter

     if ( iter > Nsweep ) then
        if ( iswitch_scf == 1 ) then
           flag_scf = .true.
        else
           flag_exit = .true.
        end if
     end if

     call watch(ct0,et0)

     if ( .not. flag_exit ) then

     esp0=esp
     do s=MSP_0,MSP_1
     do k=MBZ_0,MBZ_1
        call watcht(disp_switch,"",0)
        if ( iter == 1 .or. flag_scf ) then
#ifdef _LAPACK_
           call subspace_diag_la(k,s)
#else
           call subspace_diag_sl(k,s,disp_switch)
#endif
        end if
        call watcht(disp_switch,"diag",1)
        call conjugate_gradient(ML_0,ML_1,Nband,k,s,Ncg,iswitch_gs &
                               ,unk(ML_0,1,k,s),esp(1,k,s),res(1,k,s))
        call watcht(disp_switch,"cg  ",1)
        call gram_schmidt_t(1,Nband,k,s)
        call watcht(disp_switch,"gs  ",1)
        if ( Ndiag /= 1 ) then
           !if ( .not.flag_scf ) then
#ifdef _LAPACK_
           call subspace_diag_la(k,s)
#else
           call subspace_diag_sl(k,s,disp_switch)
#endif
           !else
           !call esp_calc &
           !     (k,s,unk(ML_0,MB_0,k,s),ML_0,ML_1,MB_0,MB_1,esp(MB_0,k,s))
           !end if
        end if
        call watcht(disp_switch,"diag",1)
     end do
     end do

     call esp_gather(Nband,Nbzsm,Nspin,esp)
     call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                    ,esp,weight_bz,occ,disp_switch)

     if ( disp_switch ) then
        write(*,'(a4,a6,a20,2a13,1x)') &
             "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
        do k=1,Nbzsm
        do n=max(1,nint(Nelectron/2)-20),min(nint(Nelectron/2)+80,Nband)
           write(*,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
           ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
        end do
        end do
        write(*,*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
        write(*,'(1x,"flag_scf,iswitch_scf=",l2,i2)') flag_scf,iswitch_scf
     end if

     call calc_with_rhoIN_total_energy(disp_switch)

     if ( flag_scf ) then
        call calc_density ! n_out
        call watcht(disp_switch,"hartree",0)
        call calc_hartree(ML_0,ML_1,MSP,rho)
        call watcht(disp_switch,"hartree",1)
        call calc_xc
        call calc_total_energy(.false.,disp_switch)
        if ( mod(imix,2) == 0 ) then
           call perform_mixing(ML_1-ML_0+1,MSP_1-MSP_0+1,rho(ML_0,MSP_0),flag_conv,disp_switch)
           call normalize_density
           m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
           call mpi_allgather(rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8,comm_spin,ierr)
           call watcht(disp_switch,"hartree",0)
           call calc_hartree(ML_0,ML_1,MSP,rho)
           call watcht(disp_switch,"hartree",1)
           call calc_xc
           do s=MSP_0,MSP_1
              Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
           end do
        else if ( mod(imix,2) == 1 ) then
           do s=MSP_0,MSP_1
              Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
           end do
           call perform_mixing(ML_1-ML_0+1,MSP_1-MSP_0+1,Vloc(ML_0,MSP_0),flag_conv,disp_switch)
        end if

!---------------------------------- LPOT2
        if ( flag_localpot2 ) then
           call localpot2_density(mm1,mm2,mm3,rho_nl)
           call localpot2_calc_eion(mm1,mm2,mm3,vion_nl,rho_nl,eion_tmp)
           call localpot2_vh(mm1,mm2,mm3,Ecut,rho_nl,vh_nl,eh_tmp)
           call localpot2_xc(mm1,mm2,mm3,rho_nl,vxc_nl,exc_tmp)
           vloc_dense=vion_nl+vh_nl+vxc_nl
           vloc_dense=beta*vloc_dense+(1.d0-beta)*vloc_dense_old
           vloc_dense_old=vloc_dense
           call test2_localpot2(mm1,mm2,mm3,vloc_dense)
           call localpot2_te(eion_tmp,eh_tmp,exc_tmp,disp_switch)
        end if
!------------------------------------------

     end if ! flag_scf

     end if ! .not.flag_exit

     if ( abs(diff_etot) <= 1.d-14 ) then
        if ( iswitch_scf == 1 ) then
           flag_scf = .true.
        else
           flag_exit = .true.
        end if
     end if

     call watch(ct1,et1)
     if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0
     call global_watch(flag_end)
     flag_exit = (flag_exit.or.flag_conv.or.flag_end.or.(iter==Diter))

     if ( flag_scf .or. Nsweep > 0 ) then
     call watcht(disp_switch,"",0)
     call write_data(disp_switch,flag_exit)
     call watcht(disp_switch,"io",1)
     end if

     if ( flag_exit ) exit

  end do ! iter

  if ( disp_switch ) then
     write(*,*) "------------ SCF result ----------"
     write(*,'(a4,a6,a20,2a13,1x)') &
          "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
     do k=1,Nbzsm
     do n=1,Nband
        write(*,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
             ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
     end do
     end do
     write(*,*) "iter,sqerr=",iter,sqerr_out(1:Nspin)
     rewind 98
     write(98,'(a4,a6,a20,2a13,1x)') &
          "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
     do k=1,Nbzsm
     do n=1,Nband
        write(98,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
             ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
     end do
     end do
  end if

  call calc_total_energy(.true.,disp_switch)

  if ( flag_end ) then
     if ( disp_switch ) write(*,*) "flag_end=",flag_end
     call end_mpi_parallel
     stop
  end if

  if ( disp_switch ) write(*,'(a40," etot(with latest wf)")') repeat("-",40)
  call calc_density
  call calc_hartree(ML_0,ML_1,MSP,rho)
  call calc_xc
  do s=MSP_0,MSP_1
     Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
  end do
  call calc_total_energy(.true.,disp_switch)

!
! --- force calculation ---
!
  if ( disp_switch ) write(*,'(a40," Force")') repeat("-",40)

  if ( SYStype == 0 ) then
     select case( pselect )
     case( 2 )
        call ps_nloc2_init_derivative
     end select
  end if

  if ( iswitch_opt == -1 ) then

     call test_force(SYStype)

  end if

! --- BAND ---
#ifndef _DRSDFT_
  if ( iswitch_band == 1 ) then
     call read_band(myrank,1)
     call band(nint(Nelectron*0.5d0),disp_switch)
  end if
#endif

  select case(iswitch_opt)
  case( 1,2 )
     call atomopt(iswitch_opt,disp_switch)
  case( 3 )
#ifdef _DRSDFT_
! --- CPMD ---
     call bomd
#else
     write(*,*) "RS-CPMD is not available for COMPLEX16"
     write(*,*) "Please re-compile the program"
#endif
  end select

! --- finalize ---

  if ( DISP_SWITCH ) then
     write(*,*) "END_PROGRAM : MAIN" 
  end if
900 continue
  call close_info
  call end_mpi_parallel

END PROGRAM Real_Space_Solid
