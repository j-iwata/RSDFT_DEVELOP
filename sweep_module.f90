MODULE sweep_module

  use parallel_module
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
  PUBLIC :: calc_sweep, init_sweep

  real(8),allocatable :: esp0(:,:,:)

  integer :: iconv_check=1
  real(8) :: Echk, Echk0
  real(8) :: tol_Echk=1.d-12
  real(8) :: tol_esp=1.d-7
  real(8) :: max_esperr
  integer :: mb_ref

CONTAINS


  SUBROUTINE init_sweep( iconv_check_in, mb_ref_in, tol_in )
    implicit none
    integer,intent(IN) :: iconv_check_in, mb_ref_in
    real(8),intent(IN) :: tol_in
    iconv_check = iconv_check_in
    mb_ref = mb_ref_in
    select case( iconv_check )
    case( 1 )
       tol_Echk = tol_in
    case( 2 )
       tol_esp = tol_in
    end select
  END SUBROUTINE init_sweep


  SUBROUTINE calc_sweep( Diter, isw_gs, ierr_out, disp_switch )
    implicit none
    integer,intent(IN)  :: Diter,isw_gs
    integer,intent(OUT) :: ierr_out
    logical,intent(IN)  :: disp_switch
    integer :: iter,s,k,n,m
    real(8) :: ct0,et0,ct1,et1
    logical :: flag_exit, flag_end, flag_conv

    flag_end  = .false.
    flag_exit = .false.
    flag_conv = .false.
    Echk      = 0.0d0
    ierr_out  = 0

    allocate( esp0(Nband,Nbzsm,Nspin) ) ; esp0=0.0d0

    do iter=1,Diter

       if ( disp_switch ) then
          write(*,'(a40," sweep_iter=",i4)') repeat("-",40),iter
       end if

       call watch(ct0,et0)

       Echk0=Echk
       esp0 =esp
       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          call watcht(disp_switch,"",0)
          call conjugate_gradient(ML_0,ML_1,Nband,k,s,Ncg,isw_gs &
                                 ,unk(ML_0,1,k,s),esp(1,k,s),res(1,k,s))
          call watcht(disp_switch,"cg",1)
          call gram_schmidt(1,Nband,k,s)
          call watcht(disp_switch,"gs",1)
          call subspace_diag(k,s)
          call watcht(disp_switch,"diag",1)
       end do
       end do

       call esp_gather(Nband,Nbzsm,Nspin,esp)

       call calc_fermi(iter,Nfixed,Nband,Nbzsm,Nspin,Nelectron,Ndspin &
                      ,esp,weight_bz,occ,disp_switch)

       call calc_with_rhoIN_total_energy( .false., Echk )

       call watch(ct1,et1)

       call conv_check( iter, flag_conv )
       call global_watch( .false., flag_end )
       flag_exit = (flag_end.or.flag_conv)

       if ( disp_switch ) call write_info_sweep(ct1-ct0,et1-et0)
! ---
       call watcht(disp_switch,"",0)
       call write_data( disp_switch, flag_exit )
       call watcht(disp_switch,"io",1)
! ---
       if ( flag_exit ) exit

    end do ! iter

    deallocate( esp0 )

    ierr_out = iter

    if ( flag_end ) then
       ierr_out = -1
       if ( myrank == 0 ) write(*,*) "flag_end=",flag_end
       return
    end if

    if ( iter > Diter ) then
       ierr_out = -2
       if ( myrank == 0 ) write(*,*) "sweep not converged"
    end if

    call gather_wf

  END SUBROUTINE calc_sweep


  SUBROUTINE conv_check( iter, flag_conv )
    implicit none
    integer,intent(IN)  :: iter
    logical,intent(OUT) :: flag_conv
    select case( iconv_check )
    case( 1 )
         call conv_check_1( iter, flag_conv )
    case( 2 )
         call conv_check_2( flag_conv )
    end select
  END SUBROUTINE conv_check

  SUBROUTINE conv_check_1( iter, flag_conv )
    implicit none
    integer,intent(IN)  :: iter
    logical,intent(OUT) :: flag_conv
    if ( iter == 1 ) Echk0=0.0d0
    flag_conv=.false.
    if ( abs(Echk-Echk0) < tol_Echk ) flag_conv=.true.
  END SUBROUTINE conv_check_1

  SUBROUTINE conv_check_2( flag_conv )
    implicit none
    logical,intent(OUT) :: flag_conv
    integer :: ierr
    real(8) :: err0
    max_esperr = maxval( abs(  esp(1:mb_ref,MBZ_0,MSP_0:MSP_1) &
                             -esp0(1:mb_ref,MBZ_0,MSP_0:MSP_1) ) )
    call mpi_allreduce(max_esperr,err0,1,MPI_REAL8,MPI_MAX,comm_spin,ierr)
    call mpi_allreduce(err0,max_esperr,1,MPI_REAL8,MPI_MAX,comm_bzsm,ierr)
    flag_conv = .false.
    if ( max_esperr < tol_esp ) flag_conv = .true.
  END SUBROUTINE conv_check_2


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
    if ( iconv_check == 1 ) then
       write(*,'(1x,"Echk,dif/tol =",g18.10,2x,g12.5," /",g12.5)') &
            Echk, Echk-Echk0, tol_Echk
    else
       write(*,'(1x,"max_esperr/tol, mb_ref =",g12.5," /",g12.5,i7)') &
            max_esperr,tol_esp,mb_ref
    end if
    write(*,*) "sum(occ)=",(sum(occ(:,:,s)),s=1,Nspin)
    write(*,*) "time(sweep)=",ct,et
  END SUBROUTINE write_info_sweep


END MODULE sweep_module
