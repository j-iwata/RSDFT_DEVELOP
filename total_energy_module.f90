MODULE total_energy_module

  use rgrid_module, only: dV
  use electron_module, only: occ
  use hamiltonian_module
  use hartree_module, only: Vh, E_hartree
  use xc_module, only: Vxc,E_exchange,E_correlation,Exc,E_exchange_exx
  use ewald_module, only: Eewald
  use wf_module, only: unk,esp
  use localpot_module, only: Vloc
  use ps_local_module, only: Vion
  use density_module, only: rho
  use parallel_module
  use fermi_module, only: Efermi,Eentropy
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1 &
                               ,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
  use fock_module

  implicit none

  PRIVATE
  PUBLIC :: Etot,calc_total_energy,calc_with_rhoIN_total_energy

  real(8) :: Etot,Ekin,Eloc,Enlc,Eeig,Eion,Fene
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

CONTAINS


  SUBROUTINE calc_total_energy(flag_recalc_esp,disp_switch)
    logical,intent(IN) :: flag_recalc_esp,disp_switch
    integer :: i,n,k,s,n1,n2,ierr
    real(8) :: s0(4),s1(4),uu
    real(8),allocatable :: esp0(:,:,:,:),esp1(:,:,:,:)
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.d0
    real(8),allocatable :: work(:,:)
#else
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: work(:,:)
#endif
    include 'mpif.h'

    Etot = 0.d0
    Ekin = 0.d0
    Eloc = 0.d0
    Enlc = 0.d0
    Eeig = 0.d0
    Eion = 0.d0
    Fene = 0.d0

    n1 = ML_0
    n2 = ML_1

    if ( flag_recalc_esp ) then

       allocate( esp0(MB,MBZ,MSP,4) ) ; esp0=0.d0

       allocate( work(n1:n2,1) )

       do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0,MB_1
          work=zero
          call op_kinetic(k,unk(n1,n,k,s),work,n1,n2,n,n)
#ifdef _DRSDFT_
          esp0(n,k,s,1)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
          esp0(n,k,s,1)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif
          work=zero
          call op_localpot(s,n2-n1+1,1,unk(n1,n,k,s),work)
#ifdef _DRSDFT_
          esp0(n,k,s,2)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
          esp0(n,k,s,2)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif
          work=zero
          call op_nonlocal(k,unk(n1,n,k,s),work,n1,n2,n,n)
#ifdef _DRSDFT_
          esp0(n,k,s,3)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
          esp0(n,k,s,3)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif
          work=zero
          call op_fock(k,s,n1,n2,n,n,unk(n1,n,k,s),work)
#ifdef _DRSDFT_
          esp0(n,k,s,4)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
          esp0(n,k,s,4)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif
       end do
       end do
       end do

       deallocate( work )

       allocate( esp1(MB,MBZ,MSP,4) )

       n=MB*MBZ*MSP*4
       call mpi_allreduce(esp0,esp1,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

       Ekin = sum( occ(:,:,:)*esp1(:,:,:,1) )
       Eloc = sum( occ(:,:,:)*esp1(:,:,:,2) )
       Enlc = sum( occ(:,:,:)*esp1(:,:,:,3) )

       esp(:,:,:)=esp1(:,:,:,1)+esp1(:,:,:,2)+esp1(:,:,:,3)+esp1(:,:,:,4)

       deallocate( esp1 )
       deallocate( esp0 )

    end if

    Eeig = sum( occ(:,:,:)*esp(:,:,:) )

    s0(:)=0.d0
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0,MB_1
       s1(:)=0.d0
       do i=n1,n2
          uu=abs( unk(i,n,k,s) )**2
          s1(1) = s1(1) + uu*Vloc(i,s)
          s1(2) = s1(2) + uu*Vion(i)
          s1(3) = s1(3) + uu*Vh(i)
          s1(4) = s1(4) + uu*Vxc(i,s)
       end do
       s0(:)=s0(:)+occ(n,k,s)*s1(:)
    end do
    end do
    end do
    s0(:)=s0(:)*dV
    call mpi_allreduce(s0,s1,4,mpi_real8,mpi_sum,MPI_COMM_WORLD,ierr)

    Eloc = s1(1)
    Eion = s1(2)

    Etot = Eeig - Eloc + E_hartree + Exc + Eion + Eewald - 2*E_exchange_exx

    Ehwf = Eeig - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald - 2*E_exchange_exx

    Fene = Etot - Eentropy

    if ( disp_switch ) then
       write(*,*) '(EII) ',Eewald
       write(*,*) '(KIN) ',Ekin, Ekin-Ekin_0
       write(*,*) '(LOC) ',Eloc, Eloc-Eloc_0
       write(*,*) '(NLC) ',Enlc, Enlc-Enlc_0
       write(*,*) '(ION) ',Eion, Eion-Eion_0
       write(*,*) '(HTR) ',E_hartree, E_hartree-Ehat_0
       write(*,*) '(EXC) ',Exc,  Exc-Exc_0
       write(*,*) '(EXG) ',E_exchange, E_exchange-Ex_0
       write(*,*) '(COR) ',E_correlation, E_correlation-Ec_0
       write(*,*) '(EIG) ',Eeig, Eeig-Eeig_0
       write(*,*) '(HWF) ',Ehwf, Ehwf-Etot
       write(*,*) '(TOT) ',Etot, Etot_0-Etot
       write(*,*) '(efermi)  ',efermi, efermi-efermi_0
       write(*,*) '(entropy) ',Eentropy,Eentropy-Eentropy_0
       write(*,*) '(FreeEne) ',Fene,Fene-Fene_0
       write(41,*) abs(Etot_0-Etot)
    end if

    Etot_0 = Etot
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

  END SUBROUTINE calc_total_energy


  SUBROUTINE calc_with_rhoIN_total_energy(disp_switch)
    logical,intent(IN) :: disp_switch
    real(8) :: sb(2),rb(2),Eeig_tmp
    integer :: s,ierr
    sb(:)=0.d0
    do s=MSP_0,MSP_1
       sb(1) = sb(1) + sum(Vloc(:,s)*rho(:,s))
       sb(2) = sb(2) + sum(Vion(:)*rho(:,s))
    end do
    call mpi_allreduce(sb,rb,2,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    call mpi_allreduce(rb,sb,2,MPI_REAL8,MPI_SUM,comm_spin,ierr)
    Eloc_in = sb(1)*dV
    Eion_in = sb(2)*dV
    Ehat_in = E_hartree
    Exc_in  = Exc
    Eeig_tmp=sum( occ(:,:,:)*esp(:,:,:) )
    Ehwf = Eeig_tmp - Eloc_in + Ehat_in + Exc_in + Eion_in + Eewald
    if ( disp_switch ) then
       write(*,*) '(HWF) ',Ehwf, Ehwf_0-Ehwf
    end if
    Ehwf_0 = Ehwf
  END SUBROUTINE calc_with_rhoIN_total_energy


END MODULE total_energy_module
