MODULE localpot2_te_module

  use wf_module
  use array_bound_module
  use rgrid_module
  use localpot_module
  use kinetic_module
  use nonlocal_module
  use ewald_module
  use electron_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_te, etot_lpot2

  real(8) :: etot_lpot2
  real(8) :: Etot_0=0.d0
  real(8) :: Ekin_0=0.d0
  real(8) :: Eloc_0=0.d0
  real(8) :: Enlc_0=0.d0
  real(8) :: Eeig_0=0.d0
  real(8) :: Eion_0=0.d0
  real(8) :: Eh_0=0.d0
  real(8) :: Exc_0 =0.d0

CONTAINS

  SUBROUTINE localpot2_te(Eion,Eh,Exc,disp_switch)
    implicit none
    real(8),intent(IN) :: Eion,Eh,Exc
    logical,intent(IN) :: disp_switch
    integer :: n1,n2,n,k,s,ierr
    real(8) :: Etot,Ekin,Eloc,Enlc,Eeig
    real(8),allocatable :: esp0(:,:,:,:),esp1(:,:,:,:)
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: work(:,:)
#else
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: work(:,:)
#endif
    include 'mpif.h'

    n1 = ML_0
    n2 = ML_1

    allocate( esp0(MB,MBZ,MSP,3) ) ; esp0=0.0d0
    allocate( esp1(MB,MBZ,MSP,3) ) ; esp1=0.0d0
    allocate( work(n1:n2,1)      ) ; work=(0.0d0,0.0d0)

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
       call op_nonlocal(k,s,unk(n1,n,k,s),work,n1,n2,n,n)
#ifdef _DRSDFT_
       esp0(n,k,s,3)=sum( unk(:,n,k,s)*work(:,1) )*dV
#else
       esp0(n,k,s,3)=sum( conjg(unk(:,n,k,s))*work(:,1) )*dV
#endif

    end do
    end do
    end do

    deallocate( work )

    n=MB*MBZ*MSP*3
    call MPI_ALLREDUCE(esp0,esp1,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

    esp(:,:,:)=esp1(:,:,:,1)+esp1(:,:,:,2)+esp1(:,:,:,3)

    Eeig = sum( occ(:,:,:)*esp(:,:,:) )
    Ekin = sum( occ(:,:,:)*esp1(:,:,:,1) )
    Eloc = sum( occ(:,:,:)*esp1(:,:,:,2) )
    Enlc = sum( occ(:,:,:)*esp1(:,:,:,3) )

    deallocate( esp1, esp0 )

    Etot_lpot2 = Eeig - Eloc + Eion + Eh + Exc + Eewald
    Etot=Etot_lpot2

    if ( disp_switch ) then
       write(*,*) "--- localpot2_total-energy ---"
       write(*,*) '(EII) ',Eewald
       write(*,*) '(KIN) ',Ekin, Ekin-Ekin_0
       write(*,*) '(LOC) ',Eloc, Eloc-Eloc_0
       write(*,*) '(NLC) ',Enlc, Enlc-Enlc_0
       write(*,*) '(ION) ',Eion, Eion-Eion_0
       write(*,*) '(HTR) ',Eh, Eh-Eh_0
       write(*,*) '(EXC) ',Exc,  Exc-Exc_0
       write(*,*) '(EIG) ',Eeig, Eeig-Eeig_0
       write(*,*) '(TOT) ',Etot, Etot_0-Etot
    end if

    Etot_0 = Etot
    Ekin_0 = Ekin
    Eloc_0 = Eloc
    Enlc_0 = Enlc
    Eion_0 = Eion
    Eh_0   = Eh
    Exc_0  = Exc
    Eeig_0 = Eeig

  END SUBROUTINE localpot2_te

END MODULE localpot2_te_module
