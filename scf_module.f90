MODULE scf_module

  use aa_module
  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module
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
  use total_energy_module
  use fermi_module
  use subspace_diag_module
  use esp_gather_module
  use density_module
  use watch_module
  use ggrid_module, only: Ecut
  use rgrid_module, only: dV, Ngrid
  use esp_calc_module

  use localpot2_variables, only: vloc_dense,vloc_dense_old,rho_nl &
                                ,vxc_nl,vh_nl,vion_nl
  use localpot2_module, only: flag_localpot2, test2_localpot2
  use localpot2_density_module, only: localpot2_density
  use localpot2_vh_module, only: localpot2_vh
  use localpot2_xc_module, only: localpot2_xc
  use localpot2_ion_module, only: localpot2_calc_eion
  use localpot2_te_module, only: localpot2_te, diff_Etot_lpot2

  use force_module, only: get_fmax_force
  use hamiltonian_module

  implicit none

  PRIVATE
  PUBLIC :: calc_scf, read_scf, Diter_scf

  integer :: Diter_scf   = 100
  integer :: Ndiag       = 1
  integer :: second_diag = 2
  real(8) :: scf_conv(2) = 0.0d0
  real(8) :: fmax_conv   = 0.d0
  real(8) :: etot_conv   = 0.d0

  real(8),allocatable :: esp0(:,:,:),Vloc0(:,:)
  real(8),allocatable :: rho_0(:,:),vxc_0(:,:),vht_0(:)
  real(8) :: diff_vrho(7)
  real(8) :: time_scf(4)
  real(8) :: fmax,fmax0,fdiff

CONTAINS


  SUBROUTINE read_scf( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(8) :: cbuf,ckey
    scf_conv(1) = 1.d-20
    scf_conv(2) = 0.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:7) == "SCFCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,scf_conv
          else if ( ckey(1:8) == "FMAXCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,fmax_conv
          else if ( ckey(1:8) == "ETOTCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,etot_conv
          else if ( ckey(1:5) == "DITER" ) then
             backspace(unit)
             read(unit,*) cbuf,Diter_scf
          else if ( ckey(1:5) == "NDIAG" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndiag,second_diag
          end if
       end do
999    continue
       write(*,*) "scf_conv    =",scf_conv
       write(*,*) "fmax_conv   =",fmax_conv
       write(*,*) "etot_conv   =",etot_conv
       write(*,*) "Diter_scf   =",Diter_scf
       write(*,*) "Ndiag       =",Ndiag
       write(*,*) "second_diag =",second_diag
    end if
    call mpi_bcast( scf_conv  ,2,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fmax_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(etot_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Diter_scf  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Ndiag      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(second_diag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE read_scf


  SUBROUTINE calc_scf( Diter, ierr_out, disp_switch, Etot_out )
    implicit none
    integer,intent(IN)  :: Diter
    integer,intent(OUT) :: ierr_out
    logical,intent(IN) :: disp_switch
    real(8),optional,intent(OUT) :: Etot_out
    integer :: iter,s,k,n,m,ierr,idiag
    integer :: ML01,MSP01,ib1,ib2,iflag_hybrid,iflag_hybrid_0
    real(8) :: ct0,et0,ct1,et1,ct(0:5),et(0:5),ctt(0:7),ett(0:7)
    logical :: flag_exit,flag_end,flag_conv,flag_conv_f

    if ( myrank == 0 ) write(*,*) "------------ SCF START ----------"

    flag_end    = .false.
    flag_exit   = .false.
    flag_conv   = .false.
    flag_conv_f = .false.
    ierr_out    = 0
    fmax0       =-1.d10
    fdiff       =-1.d10

    time_scf(:) = 0.0d0

    ML01        = ML_1-ML_0+1
    MSP01       = MSP_1-MSP_0+1
    ib1         = max(1,nint(Nelectron/2)-20)
    ib2         = min(nint(Nelectron/2)+80,Nband)

    do s=MSP_0,MSP_1
       Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
    end do

    call init_mixing(ML01,MSP,MSP_0,MSP_1,comm_grid,comm_spin &
                    ,dV,rho(ML_0,MSP_0),Vloc(ML_0,MSP_0),scf_conv &
                    ,ir_grid,id_grid,myrank)

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    if ( iflag_hunk == 1 ) then

       hunk(:,:,:,:)=0.0d0

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do m=MB_0,MB_1,MB_d
             n=min(m+MB_d-1,MB_1)
             call hamiltonian &
                  (k,s,unk(ML_0,m,k,s),hunk(ML_0,m,k,s),ML_0,ML_1,m,n)
          end do
       end do
       end do

       allocate( Vloc0(ML_0:ML_1,MSP_0:MSP_1) ) ; Vloc0=0.0d0
       Vloc0(:,:)=Vloc(:,:)

    end if

    do iter=1,Diter

       if ( disp_switch ) then
          write(*,'(a60," scf_iter=",i4)') repeat("-",60),iter
       end if

       call watch(ct0,et0)

       esp0=esp

       call get_flag_xc_hybrid( iflag_hybrid_0 )

       ct(:)=0.0d0
       et(:)=0.0d0

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1

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
                        (k,s,unk(ML_0,m,k,s),hunk(ML_0,m,k,s),ML_0,ML_1,m,n)
                   hunk(:,m:n,k,s)=hunk(:,m:n,k,s)+workwf(:,1:n-m+1)
                end do
                deallocate( workwf )
             end if

          end if

          call watchs(ct(0),et(0),0)

          call subspace_diag(k,s)

          call watchs(ct(0),et(0),1)

          call control_xc_hybrid(1)

          do idiag=1,Ndiag

             if ( Ndiag > 1 .and. disp_switch ) then
                write(*,'(a5," idiag=",i4)') repeat("-",5),idiag
             end if

             call watchs(ct(1),et(1),0)

             call conjugate_gradient(ML_0,ML_1,Nband,k,s,Ncg,iswitch_gs &
                                    ,unk(ML_0,1,k,s),esp(1,k,s),res(1,k,s))

             call watchs(ct(1),et(1),1)

             call gram_schmidt(1,Nband,k,s)

             call watchs(ct(2),et(2),1)

             if ( second_diag == 1 .or. idiag < Ndiag ) then
                call subspace_diag(k,s)
             else if ( second_diag == 2 .and. idiag == Ndiag ) then
                call esp_calc(k,s,unk(ML_0,MB_0,k,s) &
                             ,ML_0,ML_1,MB_0,MB_1,esp(MB_0,k,s))
             end if

             call watchs(ct(3),et(3),1)

          end do ! idiag

       end do ! k
       end do ! s


       ctt(0)=et(0)
       ctt(1)=et(1)
       ctt(2)=et(2)
       ctt(3)=et(3)
       call mpi_allreduce(ctt,ett(0),4,mpi_real8,mpi_min,mpi_comm_world,ierr)
       call mpi_allreduce(ctt,ett(4),4,mpi_real8,mpi_max,mpi_comm_world,ierr)

       if ( disp_switch ) then
          write(*,*) "time(diag)    =",ctt(0),ett(0),ett(4)
          write(*,*) "time(cg)      =",ctt(1),ett(1),ett(5)
          write(*,*) "time(gs)      =",ctt(2),ett(2),ett(6)
          write(*,*) "time(esp/diag)=",ctt(3),ett(3),ett(7)
       end if

       call watcht(disp_switch,"    ",0)

       call esp_gather(Nband,Nbzsm,Nspin,esp)

       call watcht(disp_switch,"esp_gather",1)

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call watcht(disp_switch,"fermi",1)

       call calc_with_rhoIN_total_energy(disp_switch)

       call watcht(disp_switch,"harris",1)

! ---
       call calc_density ! n_out
       call calc_hartree(ML_0,ML_1,MSP,rho)
       call calc_xc
       call calc_total_energy( .false., disp_switch, .true. )
       if ( present(Etot_out) ) Etot_out = Etot
! ---

       call watcht(disp_switch,"etot",1)

! --- convergence check by Fmax ---

       if ( fmax_conv > 0.0d0 ) then
          call get_fmax_force( fmax, ierr )
          fdiff=fmax-fmax0
          if ( ierr == 0 ) then
             if ( disp_switch ) then
                write(*,*) "fmax=",fmax,fdiff,flag_conv_f
             end if
             if ( abs(fdiff) < fmax_conv ) then
                flag_conv_f=.true.
             else
                fmax0 = fmax
             end if
          end if
       end if

! ---
       do s=MSP_0,MSP_1
          Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
       end do

       call perform_mixing( ML01, MSP_0, MSP_1, rho(ML_0,MSP_0) &
            ,Vloc(ML_0,MSP_0), flag_conv_f, flag_conv, disp_switch )

       if ( mod(imix,2) == 0 ) then
          call normalize_density
          m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
          call mpi_allgather &
               (rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8,comm_spin,ierr)
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_xc
          do s=MSP_0,MSP_1
             Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
          end do
       end if
! ---

       call watcht(disp_switch,"mixing",1)

       if ( flag_localpot2 ) then
          call sub_localpot2_scf( disp_switch )
          flag_conv=.false.
          if ( abs(diff_etot_lpot2) < 1.d-10 ) flag_conv=.true.
       end if

       call write_info_scf( ib1, ib2, iter, disp_switch, 0 )

       call watch(ct1,et1)
       time_scf(1)=ct1-ct0
       time_scf(2)=et1-et0
       time_scf(3)=time_scf(3)+time_scf(1)
       time_scf(4)=time_scf(4)+time_scf(2)
       if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0

       call global_watch(.false.,flag_end)

       flag_conv = (flag_conv.or.flag_conv_f)
       flag_exit = (flag_conv.or.flag_end.or.(iter==Diter))

       call watcht(disp_switch,"",0)
       call write_data(disp_switch,flag_exit)
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) then
          call finalize_mixing
          exit
       end if

    end do ! iter

    if ( myrank == 0 ) then
       write(*,*) "------------ SCF result ----------"
       call write_info_total_energy( myrank==0, .true. )
       call write_info_scf( 1, Nband, iter, myrank==0, 1 )
    end if

    if ( allocated(Vloc0) ) deallocate( Vloc0 )
    deallocate( esp0 )

    if ( .not.flag_conv ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "scf not converged"
    else
       ierr_out = iter
       if ( myrank == 0 ) then
          if ( flag_conv_f ) then
             write(*,'(A,2(E12.5,1X))') &
                  " exit SCF loop:  fmax-fmax0, fmaxconv= ",fdiff,fmax_conv
          else
             if ( nspin == 1 ) then
                write(*,'(2(A,1(E12.5,2X)),A,2E12.5)') &
                     " exit SCF loop:  rsqe=",sqerr_out(1:nspin), &
                     " vsqe=",sqerr_out(nspin+1:nspin*2)," scfconv= ",scf_conv
             else
                write(*,'(2(A,2(E12.5,1X),1X),1X,A,2E12.5)') &
                     " exit SCF loop:  rsqe=",sqerr_out(1:nspin), &
                     " vsqe=",sqerr_out(nspin+1:nspin*2),"scfconv= ",scf_conv
             end if
          end if
       end if
    end if

    if ( flag_end ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "exit SCF loop: time limit"
       return
    end if

    call gather_wf

    if ( myrank == 0 ) write(*,*) "------------ SCF END ----------"

  END SUBROUTINE calc_scf


  SUBROUTINE write_info_scf( ib1, ib2, iter, disp_switch, flag )
    implicit none
    integer,intent(IN) :: ib1, ib2, iter, flag
    logical,intent(IN) :: disp_switch
    integer :: s,k,n,nb1,nb2,i,u(2)
    u(:) = (/ 6, 99 /)
    do i=1,2
       if ( myrank /= 0 ) cycle
       if ( u(i) == 6 ) then
          nb1 = ib1
          nb2 = ib2
       else
          nb1 = 1
          nb2 = Nband
       end if
       if ( u(i) == 99 .or. disp_switch ) then
          write(u(i),'(a4,a6,a20,2a13,1x)') &
               "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
          do k=1,Nbzsm
          do n=nb1,nb2
             write(u(i),'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
                  ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
          end do
          end do
          write(u(i),*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
       end if
       if ( u(i) == 6 .and. flag == 0 ) then
          if ( NSPIN == 1 ) then
             write(u(i), '(A,I4,2(A,E12.5,2X),A,2(E12.5,2X),A,E12.5)') &
                  " iter= ",iter,"  rsqerr= ",sqerr_out(1),&
                  " vsqerr= ",sqerr_out(2), " fm,fm-fm0= ",fmax,fdiff, &
                  " time=",time_scf(2)
          else
             write(u(i), '(A,I4,4(A,2(E12.5,2x),A,E12.5))') &
                  " iter= ",iter,"  rsqerr= ",sqerr_out(1:2), &
                  " vsqerr= ",sqerr_out(3:4)," sum_dspin,|dspin|= ", &
                  sum_dspin(1:2)," fm,fm-fm0= ",fmax,fdiff, &
                  " time=",time_scf(2)
          end if
       else
          if ( NSPIN == 1 ) then
             write(u(i), '(1X,2(A,E12.5,2X),A,2(E14.5,2X))') &
                  "RSQERR= ",sqerr_out(1),  " VSQERR= ",sqerr_out(2), &
                  " FM,FM-FM0= ",fmax,fdiff
          else
             write(u(i), '(1X,4(A,2(E12.5,2x)))') &
                  "RSQERR= ",sqerr_out(1:2)," VSQERR= ",sqerr_out(3:4), &
                  " sum_DSPIN,|DSPIN|= ",sum_dspin(1:2), &
                  " FM,FM-FM0= ",fmax,fdiff
          end if
       end if
       if ( u(i) == 6 .and. .not.disp_switch ) cycle
       if ( u(i) == 99 ) then
          write(u(i),'("AX", f20.15)') ax
          write(u(i),'("A1",3f20.15)') aa(1:3,1)/ax
          write(u(i),'("A2",3f20.15)') aa(1:3,2)/ax
          write(u(i),'("A3",3f20.15)') aa(1:3,3)/ax
          write(u(i),'("VA", f30.15)') Va
          write(u(i),'("NGRID",3i5,i10)') Ngrid(1:3),Ngrid(0)
       end if
       if ( u(i) == 99 .or. flag == 1 ) then
          write(u(i),'(1X,A,f10.5,2x,A,A,2x,A,i4,2x,2(A,f11.5,2X))') &
               "ECUT=",Ecut,"XC=",XCtype,"ITER=",iter, &
               "CTIME=",time_scf(3),"ETIME=",time_scf(4)
       end if
    end do
  END SUBROUTINE write_info_scf


  SUBROUTINE sub_localpot2_scf( disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    integer :: mm1,mm2,mm3
    real(8) :: eion_tmp,eh_tmp,exc_tmp
    call localpot2_density( rho_nl )
    call localpot2_calc_eion( vion_nl, rho_nl, eion_tmp )
    call localpot2_vh( Ecut, rho_nl, vh_nl, eh_tmp )
    call localpot2_xc( rho_nl, vxc_nl, exc_tmp )
    vloc_dense=vion_nl+vh_nl+vxc_nl
    vloc_dense=beta*vloc_dense+(1.d0-beta)*vloc_dense_old
    vloc_dense_old=vloc_dense
    call test2_localpot2( vloc_dense )
    call localpot2_te( eion_tmp, eh_tmp, exc_tmp, disp_switch )
  END SUBROUTINE sub_localpot2_scf


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
    call mpi_allgather(diff_vrho(MSP_0),m,MPI_REAL8 &
                      ,diff_vrho,m,MPI_REAL8,comm_spin,ierr)
    call mpi_allgather(diff_vrho(MSP_0+2),m,MPI_REAL8 &
                      ,diff_vrho(3),m,MPI_REAL8,comm_spin,ierr)

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
    call mpi_allreduce( MPI_IN_PLACE,diff_vrho,7,MPI_REAL8 &
                       ,MPI_SUM,comm_grid,ierr )

    if ( disp_switch ) then
       do i=1,MSP
          write(*,'(1x,"diff_rho(",i1,"/",i1,")=",g12.5)') i,MSP,diff_vrho(i)
       end do
       write(*,'(1x,"diff_tod     =",g12.5)') diff_vrho(6)
       write(*,'(1x,"diff_spd     =",g12.5)') diff_vrho(7)
       do i=1,MSP
          write(*,'(1x,"diff_vxc(",i1,"/",i1,")=",g12.5)') i,MSP,diff_vrho(2+i)
       end do
       write(*,'(1x,"diff_vht     =",g12.5)'),diff_vrho(5)
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
