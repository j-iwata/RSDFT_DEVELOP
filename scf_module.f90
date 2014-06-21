MODULE scf_module

  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module
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
  use esp_calc_module

  use localpot2_variables, only: vloc_dense,vloc_dense_old,rho_nl &
                                ,vxc_nl,vh_nl,vion_nl
  use localpot2_module, only: flag_localpot2, test2_localpot2
  use localpot2_density_module, only: localpot2_density
  use localpot2_vh_module, only: localpot2_vh
  use localpot2_xc_module, only: localpot2_xc
  use localpot2_ion_module, only: localpot2_calc_eion
  use localpot2_te_module, only: localpot2_te

  implicit none

  PRIVATE
  PUBLIC :: calc_scf, init_scf

  integer :: Ndiag=1

  real(8),allocatable :: esp0(:,:,:)

!  logical :: second_diag=.false.
  logical :: second_diag=.true.

CONTAINS

  SUBROUTINE init_scf( Ndiag_in )
     implicit none
     integer,intent(IN) :: Ndiag_in
     if ( Ndiag_in > 0 ) Ndiag=Ndiag_in
  END SUBROUTINE init_scf

  SUBROUTINE calc_scf( Diter, ierr_out, disp_switch )
    implicit none
    integer,intent(IN)  :: Diter
    integer,intent(OUT) :: ierr_out
    logical,intent(IN) :: disp_switch
    integer :: iter,s,k,n,m,ierr,idiag
    integer :: ML01,MSP01,ib1,ib2
    real(8) :: ct0,et0,ct1,et1
    logical :: flag_exit,flag_end,flag_conv

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    ierr_out  = 0

    ML01      = ML_1-ML_0+1
    MSP01     = MSP_1-MSP_0+1
    ib1       = max(1,nint(Nelectron/2)-20)
    ib2       = min(nint(Nelectron/2)+80,Nband)

    if ( mod(imix,2) == 0 ) then
       call init_mixing( ML01, MSP01,  rho(ML_0,MSP_0) )
    else
       call init_mixing( ML01, MSP01, Vloc(ML_0,MSP_0) )
    end if

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

          if ( second_diag .or. idiag < Ndiag ) then
             call subspace_diag(k,s)
             call watcht(disp_switch,"diag",1)
          else if ( idiag == Ndiag ) then
             call esp_calc &
                  (k,s,unk(ML_0,MB_0,k,s),ML_0,ML_1,MB_0,MB_1,esp(MB_0,k,s))
             call watcht(disp_switch,"esp_calc",1)
          end if


          end do ! idiag

       end do
       end do

       call esp_gather(Nband,Nbzsm,Nspin,esp)

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call calc_with_rhoIN_total_energy(disp_switch)

       if ( disp_switch ) call write_info_scf( ib1, ib2 )

! ---
       call calc_density ! n_out
       call calc_hartree(ML_0,ML_1,MSP,rho)
       call calc_xc

       call calc_total_energy(.false.,disp_switch)

       if ( mod(imix,2) == 0 ) then

          call perform_mixing(ML01,MSP01,rho(ML_0,MSP_0),flag_conv,disp_switch)
          call normalize_density
          m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
          call mpi_allgather &
               (rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8,comm_spin,ierr)
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_xc
          do s=MSP_0,MSP_1
             Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
          end do

       else if ( mod(imix,2) == 1 ) then

          do s=MSP_0,MSP_1
             Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
          end do
          call perform_mixing(ML01,MSP01,Vloc(ML_0,MSP_0),flag_conv,disp_switch)

       end if

       if ( flag_localpot2 ) call sub_localpot2_scf( disp_switch )

       call watch(ct1,et1)
       if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0

       call global_watch(.false.,flag_end)

       flag_exit = (flag_conv.or.flag_end)

       call watcht(disp_switch,"",0)
       call write_data(disp_switch,flag_exit)
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) exit

    end do ! iter

    if ( myrank == 0 ) then
       write(*,*) "------------ SCF result ----------"
       call write_info_scf( 1, Nband )
       write(*,*) "iter,sqerr=",iter,sqerr_out(1:Nspin)
    end if

    deallocate( esp0 )

    ierr_out = iter

    if ( flag_end ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       return
    end if

    if ( iter > Diter ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "scf not converged"
    end if

    call gather_wf

  END SUBROUTINE calc_scf


  SUBROUTINE write_info_scf(ib1,ib2)
    implicit none
    integer,intent(IN) :: ib1,ib2
    integer :: s,k,n
    write(*,'(a4,a6,a20,2a13,1x)') &
         "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
    do k=1,Nbzsm
    do n=ib1,ib2
       write(*,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
            ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
    end do
    end do
    write(*,*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
  END SUBROUTINE write_info_scf


  SUBROUTINE sub_localpot2_scf( disp_switch )
    implicit none
    logical,intent(IN) :: disp_switch
    integer :: mm1,mm2,mm3
    real(8) :: eion_tmp,eh_tmp,exc_tmp
    call localpot2_density(mm1,mm2,mm3,rho_nl)
    call localpot2_calc_eion(mm1,mm2,mm3,vion_nl,rho_nl,eion_tmp)
    call localpot2_vh(mm1,mm2,mm3,Ecut,rho_nl,vh_nl,eh_tmp)
    call localpot2_xc(mm1,mm2,mm3,rho_nl,vxc_nl,exc_tmp)
    vloc_dense=vion_nl+vh_nl+vxc_nl
    vloc_dense=beta*vloc_dense+(1.d0-beta)*vloc_dense_old
    vloc_dense_old=vloc_dense
    call test2_localpot2(mm1,mm2,mm3,vloc_dense)
    call localpot2_te(eion_tmp,eh_tmp,exc_tmp,disp_switch)
  END SUBROUTINE sub_localpot2_scf


END MODULE scf_module
