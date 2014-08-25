MODULE scf_module

  use parallel_module
  use electron_module
  use localpot_module
  use mixing_module
  use xc_module
  use hartree_module
  use ps_local_module
  use bz_module
  use wf_module
  use cg_module
  use array_bound_module
  use gram_schmidt_t_module
  use io_module
  use total_energy_module
  use fermi_module
  use subspace_diag_la_module
  use subspace_diag_sl_module
  use esp_gather_module
  use density_module
  use watch_module
#ifdef _USPP_
  use PSnonLocDij
#endif

  implicit none

  PRIVATE
  PUBLIC :: calc_scf

  integer :: Ndiag=2

CONTAINS

  SUBROUTINE calc_scf(Diter,Nsweep,iter_final,disp_switch)
    implicit none
    integer,intent(IN)  :: Diter,Nsweep
    integer,intent(OUT) :: iter_final
    logical,intent(IN) :: disp_switch
    integer :: iter,s,k,n,m,ierr
    real(8) :: ct0,et0,ct1,et1
    real(8),allocatable :: esp0(:,:,:)
    logical :: flag_scf,flag_exit,flag_end,flag_conv

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    flag_scf  = .false. ; if ( Nsweep <= 0 ) flag_scf = .true.

    allocate( esp0(Nband,Nbzsm,Nspin) )

    do s=MSP_0,MSP_1
       Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
    end do
#ifdef _USPP_
        call getDij
#endif

    do iter=1,Diter

       if ( disp_switch ) write(*,'(a40," iter=",i4)') repeat("-",40),iter

       if ( iter > Nsweep ) flag_scf = .true.

       call watch(ct0,et0)

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
#ifdef _LAPACK_
             call subspace_diag_la(k,s)
#else
             call subspace_diag_sl(k,s,disp_switch)
#endif
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
          write(*,*) "flag_scf=",flag_scf
       end if

       call calc_with_rhoIN_total_energy(disp_switch)

       if ( flag_scf ) then
          call calc_density ! n_out
          call calc_hartree(ML_0,ML_1,MSP,rho)
          call calc_xc
          call calc_total_energy(.false.,disp_switch,iter)
          if ( mod(imix,2) == 0 ) then
             call perform_mixing(ML_1-ML_0+1,MSP_1-MSP_0+1,rho(ML_0,MSP_0) ,flag_conv,disp_switch)
             call normalize_density
             m=(ML_1-ML_0+1)*(MSP_1-MSP_0+1)
             call mpi_allgather(rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8,comm_spin,ierr)
             call calc_hartree(ML_0,ML_1,MSP,rho)
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
#ifdef _USPP_
        call getDij
#endif
       end if

       call watch(ct1,et1)
       if ( disp_switch ) write(*,*) "time(scf)",ct1-ct0,et1-et0

       call global_watch(flag_end)
       flag_exit = (flag_conv.or.flag_end)

       call watcht(disp_switch,"",0)
       call write_data(disp_switch,flag_exit)
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) exit

    end do ! iter

    if ( myrank == 0 ) then
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
    end if

    iter_final = iter

    deallocate( esp0 )

    if ( flag_end ) then
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       call end_mpi_parallel
       stop
    end if

    call gather_wf

  END SUBROUTINE calc_scf


END MODULE scf_module
