MODULE sweep_module

  use parallel_module, only: myrank, end_mpi_parallel
  use electron_module, only: Nfixed, Ndspin, Nspin, Nband, Nelectron
  use bz_module, only: weight_bz, Nbzsm
  use wf_module
  use cg_module, only: Ncg, conjugate_gradient
  use array_bound_module, only: ML_0,ML_1,MBZ_0,MBZ_1,MSP_0,MSP_1
  use gram_schmidt_module
  use io_module
  use total_energy_module, only: calc_with_rhoIN_total_energy
  use fermi_module
  use subspace_diag_module
  use esp_gather_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: calc_sweep

  real(8),allocatable :: esp0(:,:,:)

  real(8) :: Echk, Echk0
  real(8),parameter :: tol_Echk=1.d-12

CONTAINS

  SUBROUTINE calc_sweep( Diter, isw_gs, iter_final, disp_switch )
    implicit none
    integer,intent(IN)  :: Diter,isw_gs
    integer,intent(OUT) :: iter_final
    logical,intent(IN)  :: disp_switch
    integer :: iter,s,k,n,m,ierr
    real(8) :: ct0,et0,ct1,et1
    logical :: flag_exit, flag_end, flag_conv

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    Echk0     = 0.0d0

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    do iter=1,Diter

       if ( disp_switch ) then
          write(*,'(a40," sweep_iter=",i4)') repeat("-",40),iter
       end if

       call watch(ct0,et0)

       esp0=esp
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          call conjugate_gradient(ML_0,ML_1,Nband,k,s,Ncg,isw_gs &
                                 ,unk(ML_0,1,k,s),esp(1,k,s),res(1,k,s))
          call gram_schmidt(1,Nband,k,s)
          call subspace_diag(k,s)
       end do
       end do

       call esp_gather(Nband,Nbzsm,Nspin,esp)

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call calc_with_rhoIN_total_energy( .false., Echk )

       call watch(ct1,et1)
       if ( disp_switch ) call write_info_sweep(ct1-ct0,et1-et0)

! ---
       call conv_check( iter, flag_conv )

       call global_watch( .false., flag_end )

       flag_exit = (flag_end.or.flag_conv)
! ---
       call watcht(disp_switch,"",0)
       call write_data( disp_switch, flag_exit )
       call watcht(disp_switch,"io",1)

       if ( flag_exit ) exit

    end do ! iter

    iter_final = iter

    deallocate( esp0 )

    if ( flag_end ) then
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       call end_mpi_parallel
       stop
    end if

    call gather_wf

  END SUBROUTINE calc_sweep


  SUBROUTINE conv_check( iter, flag_conv )
    implicit none
    integer,intent(IN)  :: iter
    logical,intent(OUT) :: flag_conv
    if ( iter == 1 ) Echk0=0.0d0
    if ( abs(Echk-Echk0) < tol_Echk ) then
       flag_conv=.true.
    else
       flag_conv=.false.
       Echk0=Echk
    end if
  END SUBROUTINE conv_check


  SUBROUTINE write_info_sweep(ct,et)
    implicit none
    real(8),intent(IN) :: ct,et
    integer :: s,k,n
    write(*,'(a4,a6,a20,2a13,1x)') &
         "k","n","esp(n,k,s)","esp_err","occ(n,k,s)"
    do k=1,Nbzsm
    do n=max(1,nint(Nelectron/2)-20),min(nint(Nelectron/2)+80,Nband)
       write(*,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
            ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,Nspin)
    end do
    end do
    write(*,'(1x,"Echk,dif/tol =",g18.10,2x,g12.5," /",g12.5)') Echk, Echk-Echk0, tol_Echk
    write(*,*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
    write(*,*) "time(sweep)=",ct,et
  END SUBROUTINE write_info_sweep


END MODULE sweep_module
