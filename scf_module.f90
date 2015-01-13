MODULE scf_module

  use aa_module
  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module
  use xc_hybrid_module
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

  implicit none

  PRIVATE
  PUBLIC :: calc_scf, read_scf, Diter_scf

  integer :: Diter_scf   = 100
  integer :: Ndiag       = 1
  integer :: second_diag = 2
  real(8) :: scf_conv    = 1.d-20
  real(8) :: fmax_conv   = 0.d0
  real(8) :: etot_conv   = 0.d0

  real(8),allocatable :: esp0(:,:,:)
  real(8),allocatable :: rho_0(:,:),vxc_0(:,:),vht_0(:)
  real(8) :: diff_vrho(7)

CONTAINS


  SUBROUTINE read_scf( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(8) :: cbuf,ckey
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
    call mpi_bcast( scf_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fmax_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(etot_conv  ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Diter_scf  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Ndiag      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(second_diag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE read_scf


  SUBROUTINE calc_scf( Diter, ierr_out, disp_switch )
    implicit none
    integer,intent(IN)  :: Diter
    integer,intent(OUT) :: ierr_out
    logical,intent(IN) :: disp_switch
    integer :: iter,s,k,n,m,ierr,idiag
    integer :: ML01,MSP01,ib1,ib2
    real(8) :: ct0,et0,ct1,et1
    real(8) :: fmax,fmax0
    logical :: flag_exit,flag_end,flag_conv,flag_conv_f

    flag_end    = .false.
    flag_exit   = .false.
    flag_conv   = .false.
    flag_conv_f = .false.
    ierr_out    = 0
    fmax0       =-1.d10

    ML01      = ML_1-ML_0+1
    MSP01     = MSP_1-MSP_0+1
    ib1       = max(1,nint(Nelectron/2)-20)
    ib2       = min(nint(Nelectron/2)+80,Nband)

    call init_mixing(ML01,MSP,MSP_0,MSP_1,comm_grid,comm_spin &
                    ,dV,rho(ML_0,MSP_0),Vloc(ML_0,MSP_0),scf_conv)

!    call init_diff_vrho_scf

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    do iter=1,Diter

       if ( disp_switch ) then
          write(*,'(a40," scf_iter=",i4)') repeat("-",40),iter
       end if

       call watch(ct0,et0)

       esp0=esp

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1

          call watcht(disp_switch,"",0)

          call subspace_diag(k,s)

          call watcht(disp_switch,"diag",1)

          call control_xc_hybrid(1)

          do idiag=1,Ndiag

             if ( disp_switch ) then
                write(*,'(a5," idiag=",i4)') repeat("-",5),idiag
             end if

             call watcht(disp_switch,"",0)

             call conjugate_gradient(ML_0,ML_1,Nband,k,s,Ncg,iswitch_gs &
                                    ,unk(ML_0,1,k,s),esp(1,k,s),res(1,k,s))

             call watcht(disp_switch,"cg  ",1)

             call gram_schmidt(1,Nband,k,s)

             call watcht(disp_switch,"gs  ",1)

             if ( second_diag == 1 .or. idiag < Ndiag ) then
                call subspace_diag(k,s)
                call watcht(disp_switch,"diag",1)
             else if ( second_diag == 2 .and. idiag == Ndiag ) then
                call esp_calc(k,s,unk(ML_0,MB_0,k,s) &
                             ,ML_0,ML_1,MB_0,MB_1,esp(MB_0,k,s))
                call watcht(disp_switch,"esp_calc",1)
             end if

          end do ! idiag

       end do ! k
       end do ! s

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
       call control_xc_hybrid(1)
       call calc_xc
       call control_xc_hybrid(2)
!       call diff_vrho_scf( disp_switch )
       call calc_total_energy( .false., disp_switch, .true. )
! ---

       call watcht(disp_switch,"etot",1)

! --- convergence check by Fmax ---

       if ( fmax_conv > 0.0d0 ) then
          call get_fmax_force( fmax, ierr )
          if ( ierr == 0 ) then
             if ( abs(fmax-fmax0) < fmax_conv ) flag_conv_f=.true. 
             if ( disp_switch ) then
                write(*,*) "fmax=",fmax,fmax-fmax0,flag_conv_f
             end if
             fmax0 = fmax
          end if
       end if

! ---
       do s=MSP_0,MSP_1
          Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
       end do

       call perform_mixing( ML01, MSP_0, MSP_1, rho(ML_0,MSP_0) &
                           ,Vloc(ML_0,MSP_0), flag_conv, disp_switch )

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

       call write_info_scf( ib1, ib2, iter, disp_switch )

       call watch(ct1,et1)
       if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0

       call global_watch(.false.,flag_end)

       flag_exit = (flag_conv_f.or.flag_conv.or.flag_end)

       call watcht(disp_switch,"",0)
       call write_data(disp_switch,flag_exit)
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) exit

    end do ! iter

    if ( myrank == 0 ) write(*,*) "------------ SCF result ----------"
    call write_info_scf( 1, Nband, iter, myrank==0 )

    deallocate( esp0 )

!    call end_diff_vrho_scf

    ierr_out = iter

    if ( flag_end ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "time limit !!"
       return
    end if

    if ( iter > Diter ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "scf not converged"
    end if

    call gather_wf

  END SUBROUTINE calc_scf


  SUBROUTINE write_info_scf( ib1, ib2, iter, disp_switch )
    implicit none
    integer,intent(IN) :: ib1, ib2, iter
    logical,intent(IN) :: disp_switch
    integer :: s,k,n,nb1,nb2,i,u(3)
    u(:) = (/ 6, 98, 99 /)
    do i=1,3
       if ( u(i) == 6 ) then
          nb1 = ib1
          nb2 = ib2
       else
          nb1 = 1
          nb2 = Nband
       end if
       if ( u(i) == 6  .and. .not.disp_switch ) cycle
       if ( u(i) /= 6  .and. myrank /= 0 ) cycle
       if ( u(i) == 98 .and. myrank == 0 ) rewind 98
       if ( u(i) == 99 .and. myrank == 0 ) then
          write(u(i),'("AX", f20.15)') ax
          write(u(i),'("A1",3f20.15)') aa(1:3,1)/ax
          write(u(i),'("A2",3f20.15)') aa(1:3,2)/ax
          write(u(i),'("A3",3f20.15)') aa(1:3,3)/ax
          write(u(i),'("VA", f30.15)') Va
          write(u(i),'("NGRID",3i5,i10)') Ngrid(1:3),Ngrid(0)
       end if
       write(u(i),'("ECUT ",f10.5)') Ecut
       write(u(i),'("XC",a10)') XCtype
       write(u(i),*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
       write(u(i),*) "iter,sqerr=",iter,sqerr_out(1:Nspin)
       write(u(i),'(a4,a6,a20,2a13,1x)') &
            "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
       do k=1,Nbzsm
       do n=nb1,nb2
          write(u(i),'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
               ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
       end do
       end do
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
