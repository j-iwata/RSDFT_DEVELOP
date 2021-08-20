module total_energy_module

  use rgrid_module, only: dV
  use hamiltonian_module
  use hartree_variables, only: Vh, E_hartree
  use xc_module, only: Vxc,E_exchange,E_correlation,Exc,E_exchange_exx,DCxc
  use eion_module, only: Eewald
  use wf_module, only: unk,esp,occ
  use localpot_module, only: Vloc
  use ps_local_module, only: Vion, const_ps_local
  use density_module, only: rho
  use parallel_module
  use fermi_module, only: efermi,Eentropy
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1 &
                               ,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
  use fock_module
  use var_sys_parameter, only: pp_kind
  use vdw_grimme_module
  use hamiltonian_ncol_module

  implicit none

  private
  public :: calc_total_energy
  public :: calc_with_rhoIN_total_energy
  public :: Ekin,Eion,Enlc,cnst

  integer :: scf_iter_

  real(8) :: Ekin,Eloc,Enlc,Eeig,Eion,Fene,Evdw
  real(8) :: cnst
  real(8) :: Etot_0=0.d0
  real(8) :: Ekin_0=0.d0
  real(8) :: Eloc_0=0.d0
  real(8) :: Enlc_0=0.d0
  real(8) :: Eeig_0=0.d0
  real(8) :: Eion_0=0.d0
  real(8) :: Ehat_0=0.d0
  real(8) :: Exc_0 =0.d0
  real(8) :: Ex_0  =0.d0
  real(8) :: Ec_0  =0.d0

  real(8) :: efermi_0  =0.d0
  real(8) :: Eentropy_0=0.d0
  real(8) :: Fene_0    =0.d0

  real(8) :: Ehwf
  real(8) :: Ehwf_0  = 0.d0
  real(8) :: Eloc_in = 0.d0
  real(8) :: Ehat_in = 0.d0
  real(8) :: Exc_in  = 0.d0
  real(8) :: Eion_in = 0.d0

  real(8) :: diff_etot = 0.d0
  real(8) :: Etot_wo_cnst = 0.d0
  real(8) :: Etot_1 = 0.0d0

contains


  subroutine calc_total_energy &
       ( flag_recalc_esp, Etot, Free_energy, unit_in, flag_ncol )
    implicit none
    logical,intent(in) :: flag_recalc_esp
    real(8),optional,intent(inout) :: Etot
    real(8),optional,intent(inout) :: Free_energy
    integer,optional,intent(in) :: unit_in
    logical,optional,intent(in) :: flag_ncol
    integer :: i,n,k,s,n1,n2,ierr,nb1,nb2,nn,unit
    real(8) :: s0(4),s1(4),uu,cnst,c1,c2
    real(8),allocatable :: esp0(:,:,:,:),esp1(:,:,:,:)
    real(8),allocatable :: esp0_Q(:,:,:),esp1_Q(:,:,:)
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: work(:,:)
    real(8),allocatable :: work00(:,:),zw1(:,:,:),zw2(:,:,:)
#else
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: work(:,:)
    complex(8),allocatable :: work00(:,:),zw1(:,:,:),zw2(:,:,:)
#endif
    include 'mpif.h'
    complex(8) :: ztmp,ztmp1

    call write_border( 1, " calc_total_energy(start)" )

    Etot_0 = 0.0d0
    if ( present(Etot) ) Etot_0 = Etot
    Etot_1 = 0.0d0
    Ekin = 0.d0
    Eloc = 0.d0
    Enlc = 0.d0
    Eeig = 0.d0
    Eion = 0.d0
    Fene = 0.d0
    Evdw = 0.d0

    n1 = ML_0
    n2 = ML_1

    if ( flag_recalc_esp ) then

       allocate( esp0(MB,MBZ,MSP,5) ) ; esp0=0.0d0
       allocate( esp0_Q(MB,MBZ,MSP) ) ; esp0_Q=0.0d0
       allocate( work(n1:n2,MB_d)   ) ; work=zero
       allocate( work00(n1:n2,MB_d) ) ; work00=zero
       allocate( zw1(n1:n2,1,MSP) ) ; zw1=zero
       allocate( zw2(n1:n2,1,MSP) ) ; zw2=zero

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0,MB_1,MB_d

          nb1=n
          nb2=min(nb1+MB_d-1,MB_1)
          nn =nb2-nb1+1

!---------------------------------------------------- kinetic

          work=zero
!$OMP parallel
          call op_kinetic( unk(:,nb1:nb2,k,s), work(:,1:nn), k )
!$OMP end parallel
          do i=nb1,nb2
#ifdef _DRSDFT_
          esp0(i,k,s,1)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
          esp0(i,k,s,1)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
          end do

!---------------------------------------------------- local

          work=zero
!$OMP parallel
          call op_localpot( unk(:,nb1:nb2,k,s), work(:,1:nn), s )
!$OMP end parallel
          do i=nb1,nb2
#ifdef _DRSDFT_
          esp0(i,k,s,2)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
          esp0(i,k,s,2)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
          end do

!---------------------------------------------------- nonlocal

          if ( pp_kind == 'USPP' ) then

             work=zero
             work00=zero
!$OMP parallel
             call op_nonlocal( unk(:,nb1:nb2,k,s), work(:,1:nn), nb1,k,s )
!$OMP end parallel
             do i=nb1,nb2
#ifdef _DRSDFT_
                esp0(i,k,s,3)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
                esp0_Q(i,k,s)=sum( unk(:,i,k,s)*work00(:,i-nb1+1) )*dV
#else
                esp0(i,k,s,3)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
                esp0_Q(i,k,s)=sum( conjg(unk(:,i,k,s))*work00(:,i-nb1+1) )*dV
#endif
             end do

          else if ( pp_kind == 'NCPP' ) then

             work=zero
!$OMP parallel
             call op_nonlocal( unk(:,nb1:nb2,k,s), work(:,1:nn), nb1,k,s )
!$OMP end parallel
             do i=nb1,nb2
#ifdef _DRSDFT_
                esp0(i,k,s,3)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
                esp0(i,k,s,3)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
             end do

          end if

!---------------------------------------------------- fock

          work=zero
          call op_fock( unk(:,nb1:nb2,k,s), work(:,1:nn), nb1,k,s )
          
          do i=nb1,nb2
#ifdef _DRSDFT_
             esp0(i,k,s,4)=sum( unk(:,i,k,s)*work(:,i-nb1+1) )*dV
#else
             esp0(i,k,s,4)=sum( conjg(unk(:,i,k,s))*work(:,i-nb1+1) )*dV
#endif
          end do

       end do ! n
       end do ! k
       end do ! s

       deallocate( work   )
       deallocate( work00 )

       if ( present(flag_ncol) ) then
          if ( flag_ncol ) then
             do k=MBZ_0,MBZ_1
             do n=MB_0 ,MB_1
                zw1(:,1,:)=unk(:,n,k,:)
                zw2=zero
#ifndef _DRSDFT_
                call hamiltonian_ncol( k, n1,n2, zw1, zw2 )
                esp0(n,k,1,5) = sum( conjg(zw1)*zw2 )*dV
#endif
             end do ! n
             end do ! k
             esp0(:,:,MSP,5)=esp0(:,:,1,5)
          end if
       end if

       allocate( esp1(MB,MBZ,MSP,5) )
       allocate( esp1_Q(MB,MBZ,MSP) )

       n=MB*MBZ*MSP*5
       call mpi_allreduce(esp0,esp1,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       esp1=esp1/np_fkmb
       n=MB*MBZ*MSP
       call mpi_allreduce(esp0_Q,esp1_Q,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       esp1_Q=esp1_Q/np_fkmb

       Ekin = sum( occ(:,:,:)*esp1(:,:,:,1) )
       Eloc = sum( occ(:,:,:)*esp1(:,:,:,2) )
        
       if ( pp_kind == 'USPP' ) then

          Enlc = sum( occ(:,:,:)*esp1_Q(:,:,:) )

       else if ( pp_kind == 'NCPP' ) then

          Enlc = sum( occ(:,:,:)*esp1(:,:,:,3) )

       endif

       esp(:,:,:)=esp1(:,:,:,1) &
                 +esp1(:,:,:,2) &
                 +esp1(:,:,:,3) &
                 +esp1(:,:,:,4)

       if ( present(flag_ncol) ) then
          if ( flag_ncol ) then
             esp(:,:,:)=0.0d0
             do s=1,MSP
                esp(:,:,1)=esp(:,:,1)+esp1(:,:,s,1) &
                                     +esp1(:,:,s,2) &
                                     +esp1(:,:,s,3) &
                                     +esp1(:,:,s,4)
             end do
             esp(:,:,1)=esp(:,:,1)+esp1(:,:,1,5)
             esp(:,:,MSP)=esp(:,:,1)
          end if
       end if

       deallocate( zw2 )
       deallocate( zw1 )
       deallocate( esp1 )
       deallocate( esp0 )
       deallocate( esp1_Q )
       deallocate( esp0_Q )

    end if ! flag_recalc_esp

    Eeig = sum( occ(:,:,:)*esp(:,:,:) )
    if ( present(flag_ncol) ) then
       if ( flag_ncol ) Eeig = sum( occ(:,:,1)*esp(:,:,1) )
    end if

    cnst = sum( occ(:,:,:) )*const_ps_local

    if ( present(flag_ncol) ) then
       if ( flag_ncol ) cnst = sum( occ(:,:,1) )*const_ps_local
    end if

    select case( pp_kind )
    case( "USPP" )

    s0(:)=0.d0
    s1(:)=0.d0
    do s=MSP_0,MSP_1
       do i=n1,n2
          s0(1) = s0(1) + rho(i,s)*Vloc(i,s) ! rho inlculdes Qij(r) terms
          s0(2) = s0(2) + rho(i,s)*Vion(i)
       end do
    end do
    s1(:)=s0(:)*dV

    call mpi_allreduce(s1,s0,2,mpi_real8,mpi_sum,comm_grid,ierr)
    call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,comm_spin,ierr)

    case( "NCPP" )

!    s0(:)=0.d0
!    do s=MSP_0,MSP_1
!    do k=MBZ_0,MBZ_1
!    do n=MB_0,MB_1
!       s1(:)=0.d0
!       do i=n1,n2
!          uu=abs(unk(i,n,k,s))**2
!          s1(1) = s1(1) + uu*Vloc(i,s)
!          s1(2) = s1(2) + uu*Vion(i)
!          s1(3) = s1(3) + uu*Vh(i)          ! not used?
!          s1(4) = s1(4) + uu*Vxc(i,s)       ! not used?
!       end do
!       s0(:)=s0(:)+occ(n,k,s)*s1(:)
!    end do
!    end do
!    end do
!    s1(:)=s0(:)*dV
!    call mpi_allreduce(s1,s0,4,mpi_real8,mpi_sum,comm_grid,ierr)
!    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,comm_band,ierr)
!    call mpi_allreduce(s1,s0,4,mpi_real8,mpi_sum,comm_bzsm,ierr)
!    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,comm_spin,ierr)
!    s0(:)=s0(:)*dV/np_fkmb
!    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,MPI_COMM_WORLD,ierr)

       c1=0
       c2=0
       do s=MSP_0,MSP_1
          do i=n1,n2
             c1 = c1 + rho(i,s)*Vloc(i,s)
             c2 = c2 + rho(i,s)*Vion(i)
          end do
       end do
       s1(1)=c1*dV
       s1(2)=c2*dV
       call MPI_Allreduce( s1,s0,2,MPI_REAL8,MPI_SUM,comm_grid,ierr )
       call MPI_Allreduce( s0,s1,2,MPI_REAL8,MPI_SUM,comm_spin,ierr )

    end select

    Eloc = s1(1)
    Eion = s1(2)

    call get_E_vdw_grimme( Evdw )

    Etot_1 = Eeig - Eloc + E_hartree + Exc + Eion + Eewald &
           - 2*E_exchange_exx + Evdw + cnst - DCxc

    Ehwf = Eeig - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald &
         - 2*E_exchange_exx + Evdw + cnst -DCxc

    Fene = Etot_1 - Eentropy

    diff_etot = Etot_1 - Etot_0
!    diff_etot = Etot_1 - Ehwf

    if ( present(Etot) ) Etot = Etot_1
    if ( present(Free_energy) ) Free_energy = Fene

    unit=99 ; if ( present(unit_in) ) unit=unit_in
    call write_info_total_energy( Etot_1, (myrank==0), unit )

    call check_disp_length( i, 0 )
    if ( i > 1 ) call write_info_total_energy( Etot_1, (myrank==0), 6 )

!    if ( present(flag_rewind) ) then
!       call write_info_total_energy( disp_switch, flag_rewind )
!    else
!       call write_info_total_energy( disp_switch, .false. )
!    end if
!    if ( disp_switch ) then
!       write(*,'(1x,"Total Energy   =",f16.8,2x,"(Hartree)")') Etot_1
!       write(*,'(1x,"Harris Energy  =",f16.8,2x,"(Hartree)")') Ehwf
!       write(*,'(1x,"difference    =",g13.5)') Etot_1-Ehwf
!       write(*,'(1x,"Total (Harris) Energy =",f16.8,2x,"(",f16.8,")" &
!            ,2x,"(Hartree)")') Etot_1, Ehwf
!       write(*,'(1x,"difference =",g13.5)') Etot_1-Ehwf
!    end if

    Ekin_0 = Ekin
    Eloc_0 = Eloc
    Enlc_0 = Enlc
    Eion_0 = Eion
    Ehat_0 = E_hartree
    Exc_0  = Exc
    Ex_0   = E_exchange
    Ec_0   = E_correlation
    Eeig_0 = Eeig

    efermi_0   = efermi
    Eentropy_0 = Eentropy
    Fene_0     = Fene

    call write_border( 1, " calc_total_energy(end)" )

  END SUBROUTINE calc_total_energy


  subroutine calc_with_rhoIN_total_energy( Etot, flag_ncol )
    use xc_module, only: E_exchange_exx, DCxc
    use parallel_module, only: comm_grid, comm_spin
    implicit none
    real(8),optional,intent(out) :: Etot
    logical,optional,intent(in)  :: flag_ncol
    real(8) :: sb(2),rb(2),Eeig_tmp
    integer :: s,ierr
    include 'mpif.h'
    call write_border( 1, " calc_with_rhoIN_total_energy(start)" )
    sb(:)=0.0d0
    do s=lbound(Vloc,2),ubound(Vloc,2)
      sb(1) = sb(1) + sum(Vloc(:,s)*rho(:,s))
      sb(2) = sb(2) + sum(Vion(:)*rho(:,s))
    end do
    call MPI_Allreduce(sb,rb,2,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    call MPI_Allreduce(rb,sb,2,MPI_REAL8,MPI_SUM,comm_spin,ierr)
    Eloc_in = sb(1)*dV
    Eion_in = sb(2)*dV
    Ehat_in = E_hartree
    Exc_in  = Exc
    Eeig_tmp=sum( occ(:,:,:)*esp(:,:,:) )
    cnst=sum(occ)*const_ps_local
    if ( present(flag_ncol) ) then
      if ( flag_ncol ) then
        Eeig_tmp=sum( occ(:,:,1)*esp(:,:,1) )
        cnst=sum(occ(:,:,1))*const_ps_local
      end if
    end if
    call get_E_vdw_grimme( Evdw )
    Etot = Eeig_tmp - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald &
         - 2*E_exchange_exx + cnst + Evdw - DCxc
    Etot_wo_cnst = Etot - cnst
    call write_border( 1, " calc_with_rhoIN_total_energy(end)" )
  end subroutine calc_with_rhoIN_total_energy


  SUBROUTINE write_info_total_energy( Etot, flag_write, u )
    implicit none
    real(8),intent(IN) :: Etot
    logical,intent(IN) :: flag_write
    integer,intent(IN) :: u
    if ( flag_write ) then
       if ( u == 6 ) then
          call write_border( 0, " total energy" )
       else
          rewind u
       end if
       write(u,*) "Total Energy ",Etot,diff_etot
       write(u,*) "Harris Energy",Ehwf
       write(u,*) "Free Energy  ",Fene
       write(u,*) "Ion-Ion                    ",Eewald
       write(u,*) "Local Potential            ",Eloc
       write(u,*) "Ion Local Potential        ",Eion
       if ( Enlc /= 0.0d0 ) then
          write(u,*) "Ion Nonlocal Potential     ",Enlc
       else
          write(u,*) "Ion Nonlocal Potential     ","Not computed"
       end if
       if ( Ekin /= 0.0d0 ) then
          write(u,*) "Kinetic Energy             ",Ekin
       else
          write(u,*) "Kinetic Energy             ","Not computed"
       end if
       write(u,*) "Hartree Energy             ",E_hartree
       write(u,*) "Exchange-Correlation Energy",Exc
       write(u,*) "Exchange Energy            ",E_exchange
       write(u,*) "Correlation Energy         ",E_correlation
       write(u,*) "Sum of eigenvalues         ",Eeig
       write(u,*) "Fermi energy               ",efermi
       write(u,*) "Entropy                    ",Eentropy
       write(u,*) "Constant of local part     ",const_ps_local*sum(occ)
       write(u,*) "Total Energy (wo constant) ",Etot_wo_cnst
       write(u,*) "Dispersion energy          ",Evdw
       if ( u == 6 ) call write_border( 0, "" )
    end if
!    u(:) = (/ 6, 99 /)
!    do i=1,2
!       if ( u(i) == 6 .and. .not.disp_switch ) cycle
!       if ( u(i) /= 6 .and. myrank /= 0 ) cycle
!       if ( u(i) /= 6 .and. myrank == 0 .and. flag_rewind ) rewind u(i)
!       if ( Evdw /= 0.0d0 ) write(u(i),*) '(VDW) ',Evdw
!       if ( const_ps_local /= 0.0d0 ) write(u(i),*) '(cnst)',const_ps_local*sum(occ)
!       write(u(i),*) '(EII) ',Eewald
!       write(u(i),*) '(KIN) ',Ekin, Ekin-Ekin_0
!       write(u(i),*) '(LOC) ',Eloc, Eloc-Eloc_0
!       write(u(i),*) '(NLC) ',Enlc, Enlc-Enlc_0
!       write(u(i),*) '(ION) ',Eion, Eion-Eion_0
!       write(u(i),*) '(HTR) ',E_hartree, E_hartree-Ehat_0
!       write(u(i),*) '(EXC) ',Exc,  Exc-Exc_0
!       write(u(i),*) '(EXG) ',E_exchange, E_exchange-Ex_0
!       write(u(i),*) '(COR) ',E_correlation, E_correlation-Ec_0
!       write(u(i),*) '(EIG) ',Eeig, Eeig-Eeig_0
!!       write(u(i),*) '(HWF) ',Ehwf, Ehwf-Etot
!!       write(u(i),*) '(TOT) ',Etot, Etot_0-Etot
!       write(u(i),*) '(efermi)  ',efermi, efermi-efermi_0
!       !write(u(i),*) '(entropy) ',Eentropy,Eentropy-Eentropy_0
!       !write(u(i),*) '(FreeEne) ',Fene,Fene-Fene_0
!       call flush(u(i))
!    end do
  END SUBROUTINE write_info_total_energy


END MODULE total_energy_module
