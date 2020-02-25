MODULE scf_module

  use aa_module
  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module, only: init_mixing, perform_mixing, calc_sqerr_mixing &
                         , finalize_mixing, imix, beta
  use sqerr_module
  use xc_hybrid_module, only: control_xc_hybrid, get_flag_xc_hybrid
  use xc_module
  use hartree_variables, only: Vh
  use hartree_module, only: calc_hartree
  use ps_local_module
  use bz_module
  use wf_module
  use cg_module
  use array_bound_module
  use gram_schmidt_module
  use io_module
  use total_energy_module, only: calc_total_energy,calc_with_rhoIn_total_energy
  use fermi_module
  use fermi_ncol_module
  use subspace_diag_module
  use esp_gather_module
  use density_module
  use watch_module
  use ps_getDij_module
  use ggrid_module, only: Ecut
  use rgrid_module, only: dV, Ngrid
  use esp_calc_module
  use force_module, only: get_fmax_force
  use hamiltonian_module
  use io_tools_module
  use eigenvalues_module
  use rsdft_mpi_module
  use noncollinear_module
  use phase_control_module

  implicit none

  PRIVATE
  PUBLIC :: calc_scf, read_scf, Diter_scf

  integer :: Diter_scf   = 100
  integer :: Ndiag       = 1
  integer :: second_diag = 2
  real(8) :: scf_conv(2) = 0.0d0
  real(8) :: fmax_conv   = 0.0d0
  real(8) :: etot_conv   = 0.0d0
  logical :: errlog_scf  = .false.
  logical :: nodiag_scf  = .false.
  integer,parameter :: unit_err=40

  real(8),allocatable :: Vloc0(:,:)
  real(8),allocatable :: rho_0(:,:),vxc_0(:,:),vht_0(:)
  real(8) :: diff_vrho(7)
  real(8) :: time_scf(4)
  real(8) :: fmax,fmax0,fdiff

CONTAINS


  SUBROUTINE read_scf
    implicit none
    integer :: itmp(2)
    call IOTools_readReal8Keyword( "SCFCONV" , scf_conv )
    call IOTools_readReal8Keyword( "FMAXCONV", fmax_conv )
    call IOTools_readReal8Keyword( "ETOTCONV", etot_conv )
    call IOTools_readIntegerKeyword( "DITER", Diter_scf )
    call IOTools_findKeyword( "ERRLOG", errlog_scf )
    call IOTools_findKeyword( "NODIAG", nodiag_scf )
    itmp=-1
    call IOTools_readIntegerKeyword( "NDIAG", itmp )
    if ( itmp(1) > 0 ) Ndiag=itmp(1)
    if ( itmp(2) > 0 ) second_diag=itmp(2)
  END SUBROUTINE read_scf


  SUBROUTINE calc_scf(ierr_out,Diter_in,tol_force_in,outer_loop_info,Etot_out)
    implicit none
    integer,intent(OUT) :: ierr_out
    integer,optional,intent(IN) :: Diter_in
    real(8),optional,intent(IN) :: tol_force_in
    character(*),optional,intent(IN)  :: outer_loop_info
    real(8),optional,intent(OUT) :: Etot_out
    real(8),allocatable :: v(:,:,:)
    real(8) :: tol_force
    real(8) :: Etot, Ehwf, diff_etot
    real(8) :: Ntot(4), sqerr_out(4), t_out(2,14), t_ini(2), t_tmp(2)
    integer :: iter,s,k,n,m,ierr,idiag,i,j,Diter
    integer :: ML01,MSP01,ib1,ib2,iflag_hybrid,iflag_hybrid_0
    logical :: flag_exit,flag_conv,flag_conv_f,flag_conv_e
    logical :: flag_end, flag_end1, flag_end2
    logical :: flag_recalc_esp = .false.
    logical :: disp_switch
    logical,external :: exit_program
    character(40) :: chr_iter
    character(22) :: add_info
    type(eigv) :: eval
    type(time) :: etime, etime_tot, etime_lap(10)
    type(time) :: tt

    call write_border( 0, "" )
    call write_border( 0, " SCF START -----------" )

    call check_disp_switch( disp_switch, 0 )

    call init_time_watch( etime_tot )
    call init_time_watch( etime_lap(1) )

    flag_end    = .false.
    flag_end1   = .false.
    flag_end2   = .false.
    flag_exit   = .false.
    flag_conv   = .false.
    flag_conv_f = .false.
    flag_conv_e = .false.
    ierr_out    = 0
    fmax0       =-1.d10
    fdiff       =-1.d10
    tol_force   = 0.0d0 ; if ( present(tol_force_in) ) tol_force=tol_force_in
    Diter       = Diter_scf ; if ( present(Diter_in) ) Diter=Diter_in
    Etot        = 0.0d0
    time_scf(:) = 0.0d0
    t_out       = 0.0d0

    add_info ="" ; if ( present(outer_loop_info) ) add_info=outer_loop_info

    ML01        = ML_1-ML_0+1
    MSP01       = MSP_1-MSP_0+1
    ib1         = max(1,nint(Nelectron/2)-10)
    ib2         = min(nint(Nelectron/2)+10,Nband)

    if ( flag_noncollinear ) then
       !call calc_xc_noncollinear( unk, occ(:,:,1), vxc_out=Vxc )
       !occ(:,:,2) = occ(:,:,1)
!       Vxc=0.0d0
       E_exchange=0.0d0
       E_correlation=0.0d0
    end if

    do s=MSP_0,MSP_1
       Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
    end do

    call init_mixing(ML01,MSP,MSP_0,MSP_1,comm_grid,comm_spin &
                    ,dV,rho(ML_0,MSP_0),Vloc(ML_0,MSP_0) &
                    ,ir_grid,id_grid,myrank)

    call init_sqerr( ML01, MSP01, MSP, Vloc(ML_0,MSP_0), rho(ML_0,MSP_0) )

    rho_in = rho

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    if ( iflag_hunk == 1 ) then

       hunk(:,:,:,:)=0.0d0

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do m=MB_0,MB_1,MB_d
             n=min(m+MB_d-1,MB_1)
             call hamiltonian &
                  (k,s,unk(:,m:n,k,s),hunk(:,m:n,k,s),ML_0,ML_1,m,n)
          end do
       end do
       end do

       allocate( Vloc0(ML_0:ML_1,MSP_0:MSP_1) ) ; Vloc0=0.0d0
       Vloc0(:,:)=Vloc(:,:)

    end if

    call calc_time_watch( etime_lap(1) )

    do iter=1,Diter

       write(chr_iter,'(" scf_iter=",i4,1x,a)') iter, add_info
!       call write_border( 0, chr_iter(1:len_trim(chr_iter)) )

!       call start_timer( tt )
!       call init_time_watch( etime )
!       call init_time_watch( etime_lap(2) )
!       call watchb( t_ini(1) ) ; t_out=0.0d0

       esp0=esp

       call get_flag_xc_hybrid( iflag_hybrid_0 )

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1

!          call watchb( t_tmp )

          call control_xc_hybrid( iflag_hybrid_0 )

          if ( iflag_hunk == 1 ) then

             do n=MB_0,MB_1
                hunk(:,n,k,s) = hunk(:,n,k,s) &
                     + ( Vloc(:,s)-Vloc0(:,s) )*unk(:,n,k,s)
             end do
             Vloc0(:,s)=Vloc(:,s)

          else if ( iflag_hunk == 2 ) then

             call get_flag_xc_hybrid( iflag_hybrid )

             if ( iflag_hybrid == 2 ) then
                call control_xc_hybrid(0)
                allocate( workwf(ML_0:ML_1,MB_d) ) ; workwf=0.0d0
                do m=MB_0,MB_1,MB_d
                   n=min(m+MB_d-1,MB_1)
                   workwf(:,1:n-m+1)=hunk(:,m:n,k,s)
                   call hamiltonian &
                        (k,s,unk(:,m:n,k,s),hunk(:,m:n,k,s),ML_0,ML_1,m,n)
                   hunk(:,m:n,k,s)=hunk(:,m:n,k,s)+workwf(:,1:n-m+1)
                end do
                deallocate( workwf )
             end if

          end if

!          call watchb( t_tmp, t_out(:,1) )

          if ( .not.nodiag_scf ) call subspace_diag( k,s,ML_0,ML_1,unk,esp )

!          call watchb( t_tmp, t_out(:,2) )

#ifdef _DRSDFT_
          call mpi_bcast( unk,size(unk),MPI_REAL8,0,comm_fkmb,ierr )
#else
          call mpi_bcast( unk,size(unk),RSDFT_MPI_COMPLEX16,0,comm_fkmb,ierr )
#endif
          call mpi_bcast( esp,size(esp),MPI_REAL8, 0, comm_fkmb,ierr )

          call control_xc_hybrid(1)

!          call watchb( t_tmp, t_out(:,3) )

          do idiag=1,Ndiag

             if ( Ndiag > 1 .and. disp_switch ) then
                write(*,'(a5," idiag=",i4)') repeat("-",5),idiag
             end if

!             call watchb( t_tmp )

             call conjugate_gradient(ML_0,ML_1,Nband,k,s,unk,esp,res)

!             call watchb( t_tmp, t_out(:,4) )

             call gram_schmidt(1,Nband,k,s)

!             call watchb( t_tmp, t_out(:,5) )

             if ( second_diag == 1 .or. idiag < Ndiag ) then
                call subspace_diag( k,s,ML_0,ML_1,unk,esp )
!                call watchb( t_tmp, t_out(:,2) )
             else if ( second_diag == 2 .and. idiag == Ndiag ) then
                call esp_calc(k,s,ML_0,ML_1,MB_0,MB_1,unk,esp)
!                call watchb( t_tmp, t_out(:,6) )
             end if

          end do ! idiag

       end do ! k

       if ( flag_noncollinear ) exit

       end do ! s

       call phase_control( ML_0, ML_1, unk )

!       call calc_time_watch( etime_lap(2) )
!       call init_time_watch( etime_lap(3) )
!       call watchb( t_tmp )

       call esp_gather(Nband,Nbzsm,Nspin,esp)

#ifdef _DRSDFT_
       call mpi_bcast( unk,size(unk),MPI_REAL8,0,comm_fkmb,ierr )
#else
       call mpi_bcast( unk,size(unk),RSDFT_MPI_COMPLEX16,0,comm_fkmb,ierr )
#endif
       call mpi_bcast( esp,size(esp),MPI_REAL8, 0, comm_fkmb,ierr )

       if ( flag_noncollinear ) then
          call calc_fermi_ncol(iter,Nfixed,Nband,Nbzsm,Nspin &
               ,Nelectron,Ndspin,esp,weight_bz,occ )
       else
          call calc_fermi(iter,Nfixed,Nelectron,Ndspin,esp,weight_bz,occ)
       end if

!       call calc_time_watch( etime_lap(3) )
!       call init_time_watch( etime_lap(4) )
!       call watchb( t_tmp, t_out(:,7) )

! --- total energy ---

       call calc_with_rhoIN_total_energy( Ehwf, flag_noncollinear )

       if ( flag_noncollinear ) then
#ifndef _DRSDFT_
          call calc_xc_noncollinear( unk, occ(:,:,1), rho, Vxc )
#endif
          occ(:,:,2)=occ(:,:,1)
!          Vxc=0.0d0
          E_exchange=0.0d0
          E_correlation=0.0d0
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_spin_density( rho, Ntot )
       else
          call calc_density( Ntot )
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_xc
       end if
       call calc_total_energy( flag_recalc_esp, Etot, &
                               flag_ncol=flag_noncollinear )
       if ( present(Etot_out) ) Etot_out = Etot

! --- convergence check by total energy ---

       diff_etot = Etot - Ehwf

       if ( abs(diff_etot) < etot_conv ) flag_conv_e = .true.

! ---

       if ( disp_switch ) then
          write(*,*)
          write(*,'(1x,"Total Energy :",f26.8,2x,"(Hartree)")') Etot
          write(*,'(1x,"Harris Energy:",f26.8,2x,"(Hartree)")') Ehwf
          write(*,'(1x," ( diff/tol =",es13.5," /",es12.5," )",l5)') &
               diff_etot,etot_conv,flag_conv_e
       end if
       if ( disp_switch ) then
          write(*,*)
          write(*,'(1x,"Charge:",2f22.15)') (Ntot(s),s=1,Nspin)
          if ( Nspin == 2 ) then
             write(*,'(1x,"N_up - N_down        :",f22.15)') Ntot(3)
             write(*,'(1x,"sum|rho_up-rho_down| :",f22.15)') Ntot(4)
          end if
       end if

!       call calc_time_watch( etime_lap(4) )
!       call init_time_watch( etime_lap(5) )
!       call watchb( t_tmp, t_out(:,8) )

! --- convergence check by density & potential ---

       allocate( v(ML_0:ML_1,MSP_0:MSP_1,2) ) ; v=0.0d0
       do s=MSP_0,MSP_1
          v(:,s,1) = Vion(:) + Vh(:) + Vxc(:,s)
          v(:,s,2) = rho(:,s)
       end do
       call calc_sqerr( size(v,1), size(v,2), MSP, 2, v, sqerr_out )
       deallocate( v )

       if ( maxval( sqerr_out(1:MSP) ) <= scf_conv(1) .or. &
            maxval( sqerr_out(MSP+1:2*MSP) ) <= scf_conv(2) ) flag_conv=.true.

       if ( disp_switch ) then
          write(*,*)
          if ( Nspin == 1 ) then
             write(*,'(1x,"|Vloc-Vloc0|^2/tol =",es14.5,2x," /",es12.5,l5)') &
                  (sqerr_out(i),i=1,MSP),scf_conv(1),flag_conv
             write(*,'(1x,"|rho -rho0 |^2/tol =",es14.5,2x," /",es12.5)') &
                  (sqerr_out(i),i=MSP+1,2*MSP),scf_conv(2)
          else if ( Nspin == 2 ) then
             write(*,'(1x,"|Vloc-Vloc0|^2/tol =",2es14.5,2x," /",es12.5,l5)') &
                  (sqerr_out(i),i=1,MSP),scf_conv(1),flag_conv
             write(*,'(1x,"|rho -rho0 |^2/tol =",2es14.5,2x," /",es12.5)') &
                  (sqerr_out(i),i=MSP+1,2*MSP),scf_conv(2)
          end if
       end if

!       call calc_time_watch( etime_lap(5) )
!       call init_time_watch( etime_lap(6) )
!       call watchb( t_tmp, t_out(:,9) )

! --- convergence check by Fmax ---

       if ( fmax_conv > 0.0d0 ) then
          call get_fmax_force( fmax )
          fdiff=fmax-fmax0
          if ( fmax > tol_force .and. abs(fdiff) < fmax_conv ) then
             flag_conv_f=.true.
          else
             fmax0 = fmax
          end if
          if ( disp_switch ) then
             write(*,*)
             write(*,'(1x,"Fmax/fmax_tol:",es12.5," /",es12.5,4x,"(Hartree/bohr)")') &
                  fmax, tol_force
             write(*,'(1x," ( diff/tol =",es13.5," /",es12.5," )",l5)') &
                  fdiff,fmax_conv,flag_conv_f
             !write(11,*) iter,fmax
          end if
       end if

       flag_conv = ( flag_conv .or. flag_conv_f .or. flag_conv_e )

!       call calc_time_watch( etime_lap(6) )
!       call watchb( t_tmp, t_out(:,10) )

! ---

!       call global_watch( .false., flag_end1 )

       flag_end2 = exit_program()

       flag_end = ( flag_end1 .or. flag_end2 )

       flag_exit = (flag_conv.or.flag_end.or.(iter==Diter))

!       call init_time_watch( etime_lap(7) )
!       call watchb( t_tmp, t_out(:,11) )

! --- mixing ---

       if ( .not.flag_exit ) then

          do s=MSP_0,MSP_1
             Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
          end do

          if ( .not.flag_noncollinear ) then
             call perform_mixing( ML01, MSP_0, MSP_1, rho(ML_0,MSP_0) &
                                 ,Vloc(ML_0,MSP_0), disp_switch )
          end if

          if ( mod(imix,2) == 0 ) then
             if ( flag_noncollinear ) then
             else
                call normalize_density( rho )
                m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
                call rsdft_allgather( rho(:,MSP_0:MSP_1), rho, comm_spin )
                call calc_hartree(ML_0,ML_1,MSP,rho)
                call calc_xc
             end if
             do s=MSP_0,MSP_1
                Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
             end do
          end if

          call getDij

          rho_in = rho

       end if

!       call calc_time_watch( etime_lap(7) )
!       call watchb( t_tmp, t_out(:,12) )

! ---

       call write_esp_wf
       call construct_eigenvalues( Nband, Nbzsm, Nspin, esp, eval )
       if ( myrank == 0 ) call write_eigenvalues( eval )

!       call calc_time_watch( etime )
!       call watchb( t_tmp, t_out(:,13) )

       call write_data(disp_switch,flag_exit)
       if ( flag_noncollinear ) call io_write_noncollinear( myrank,flag_exit )
       call write_info_scf( (myrank==0) )
       call write_err_info(iter,sqerr_out(1:2*Nspin),diff_etot,flag_exit,(myrank==0) )

!       call result_timer( "scf", tt )

!       call watchb( t_ini, t_out(:,14) )

       if ( disp_switch ) then
       write(*,'(1x,"elapsed_time ",f8.3,"(scf)" &
                                   ,f8.3,"(sd) " &
                                   ,f8.3,"(esp)" &
                                   ,f8.3,"(cg) " &
                                   ,f8.3,"(gs) " &
                                   ,f8.3,"(for)" &
                                   ,f8.3,"(oth)")') &
                                   t_out(2,14) &
                                   ,t_out(2,2) &
                                   ,t_out(2,6) &
                                   ,t_out(2,4),t_out(2,5) &
                                   ,t_out(2,10) &
                                   ,t_out(2,1)+t_out(2,3)+sum( t_out(2,6:9) )+sum( t_out(2,11:13) )
       end if

       if ( flag_exit ) then
          call finalize_mixing
          exit
       end if

    end do ! iter

    if ( allocated(Vloc0) ) deallocate( Vloc0 )
    deallocate( esp0 )

    if ( .not.flag_conv ) then
       ierr_out = -2
    else
       ierr_out = iter
    end if

    if ( flag_end1 ) then
       ierr_out = -1
    else if ( flag_end2 ) then
       ierr_out = -3
    else
       call gather_wf
    end if

    call calc_time_watch( etime_tot )

    call write_border( 0, " SCF END -----" )
    if ( disp_switch ) then
       if ( ierr_out >= 0 ) then
          write(*,'(1x,"scf converged ( total # of iteration =",i4," )")') &
               ierr_out
       else if ( ierr_out == -1 ) then
          write(*,*) "Exit SCF loop due to the time limit"
       else if ( ierr_out == -2 ) then
          write(*,'(1x,"scf not converged !!!")')
       else if ( ierr_out == -3 ) then
          write(*,*) "'EXIT' file was found"
       end if
       write(*,'(1x,"Total SCF time :",f10.3,"(rank0)",f10.3,"(min)",f10.3,"(max)")') &
       etime_tot%t0, etime_tot%tmin, etime_tot%tmax
    end if
    call write_border( 0, "" )

  END SUBROUTINE calc_scf


  SUBROUTINE write_info_scf( flag )
    implicit none
    logical,intent(IN) :: flag
    integer :: s,k,n
    integer,parameter :: u=99
!    call write_border( 1, " write_info_scf(start)" )
    if ( flag ) then
       write(u,*) "Eigenvalues"
       write(u,'(a4,a6,a20,2a13,1x)') &
            "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
       do k=1,Nbzsm
       do n=1,Nband
          write(u,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
               ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
       end do
       end do
    end if
!    call write_border( 1, " write_info_scf(end)" )
  END SUBROUTINE write_info_scf


  SUBROUTINE write_err_info( iter, sqerr, dEtot, flag_finalize, flag )
    implicit none
    integer,intent(IN) :: iter
    real(8),intent(IN) :: sqerr(:), dEtot
    logical,intent(IN) :: flag_finalize, flag
    logical :: flag_new=.true.
    if ( .not.errlog_scf .or. .not.flag ) return
    if ( flag_new ) then
       open(unit_err,file="errlog_scf")
       flag_new=.false.
    else
       open(unit_err,file="errlog_scf",position="append")
    end if
    write(unit_err,'(1x,i6,5f18.10)') iter,log10(sqerr(:)),dEtot
    if ( flag_finalize ) write(unit_err,'(/)')
    close(unit_err)
  END SUBROUTINE write_err_info


!  SUBROUTINE sub_localpot2_scf( disp_switch )
!    implicit none
!    logical,intent(IN) :: disp_switch
!    integer :: mm1,mm2,mm3
!    real(8) :: eion_tmp,eh_tmp,exc_tmp
!    call localpot2_density( rho_nl )
!    call localpot2_calc_eion( vion_nl, rho_nl, eion_tmp )
!    call localpot2_vh( Ecut, rho_nl, vh_nl, eh_tmp )
!    call localpot2_xc( rho_nl, vxc_nl, exc_tmp )
!    vloc_dense=vion_nl+vh_nl+vxc_nl
!    vloc_dense=beta*vloc_dense+(1.d0-beta)*vloc_dense_old
!    vloc_dense_old=vloc_dense
!    call test2_localpot2( vloc_dense )
!    call localpot2_te( eion_tmp, eh_tmp, exc_tmp, disp_switch )
!  END SUBROUTINE sub_localpot2_scf


  SUBROUTINE init_diff_vrho_scf
    implicit none
    allocate( rho_0(ML_0:ML_1,MSP)         ) ; rho_0=0.0d0
    allocate( vxc_0(ML_0:ML_1,MSP_0:MSP_1) ) ; vxc_0=0.0d0
    allocate( vht_0(ML_0:ML_1)             ) ; vht_0=0.0d0
  END SUBROUTINE init_diff_vrho_scf

  SUBROUTINE diff_vrho_scf( disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    integer :: i,j,m,ierr
    real(8) :: dtmp_0,dtmp_1

    diff_vrho(:)=0.0d0

    do i=MSP_0,MSP_1
       diff_vrho(i+0) = sum( (rho(:,i)-rho_0(:,i))**2 )
       diff_vrho(i+2) = sum( (Vxc(:,i)-vxc_0(:,i))**2 )
    end do

    m=MSP_1-MSP_0+1
    call rsdft_allgather( diff_vrho(MSP_0:MSP_1), diff_vrho, comm_spin )
    call rsdft_allgather( diff_vrho(MSP_0+2:MSP_1+2),diff_vrho(3:4),comm_spin )

    diff_vrho(5) = sum( (Vh(:)-vht_0(:))**2 )

    do i=ML_0,ML_1

       dtmp_1=0.0d0
       dtmp_0=0.0d0
       do j=1,MSP
          dtmp_1 = dtmp_1 + rho(i,j)
          dtmp_0 = dtmp_0 + rho_0(i,j)
       end do
       diff_vrho(6) = diff_vrho(6) + (dtmp_1-dtmp_0)**2

       dtmp_1=0.0d0
       dtmp_0=0.0d0
       do j=1,MSP
          dtmp_1 = dtmp_1 - rho(i,j)
          dtmp_0 = dtmp_0 - rho_0(i,j)
       end do
       diff_vrho(7) = diff_vrho(7) + (dtmp_1-dtmp_0)**2

    end do ! i

    diff_vrho(:)=diff_vrho(:)/ML
    call rsdft_allreduce_sum( diff_vrho, comm_grid )

    if ( disp_switch ) then
       do i=1,MSP
          write(*,'(1x,"diff_rho(",i1,"/",i1,")=",g12.5)') i,MSP,diff_vrho(i)
       end do
       write(*,'(1x,"diff_tod     =",g12.5)') diff_vrho(6)
       write(*,'(1x,"diff_spd     =",g12.5)') diff_vrho(7)
       do i=1,MSP
          write(*,'(1x,"diff_vxc(",i1,"/",i1,")=",g12.5)') i,MSP,diff_vrho(2+i)
       end do
       write(*,'(1x,"diff_vht     =",g12.5)') diff_vrho(5)
    end if

    do i=MSP_0,MSP_1
       rho_0(:,i) = rho(:,i)
       vxc_0(:,i) = Vxc(:,i)
    end do
    vht_0(:) = Vh(:)

  END SUBROUTINE diff_vrho_scf

  SUBROUTINE end_diff_vrho_scf
    implicit none
    deallocate( vht_0 )
    deallocate( vxc_0 )
    deallocate( rho_0 )
  END SUBROUTINE end_diff_vrho_scf


END MODULE scf_module
