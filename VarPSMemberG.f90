MODULE VarPSMemberG
  implicit none

  integer,allocatable :: nl(:)
  integer,allocatable :: nr(:,:)
  integer,allocatable :: nl3v(:,:)
  integer,allocatable :: l2v(:,:,:)
  
  real(8),allocatable :: qrL(:,:,:,:)
  real(8),allocatable :: dqrL(:,:,:,:,:)

  complex(8),allocatable :: QG(:,:,:)

  real(8),allocatable :: QRij(:,:)

  integer,allocatable :: no(:,:)

  real(8),allocatable :: ddi(:,:,:,:)

  real(8),allocatable :: qqr(:,:,:,:)
  real(8),allocatable :: qqc(:,:,:,:)

  real(8),allocatable :: rabr2(:)

  integer :: k1max,k2max,k3max
  
  integer,allocatable :: k1_to_k2(:,:)
  integer,allocatable :: k1_to_k3(:,:)
  integer,allocatable :: k1_to_iorb(:,:,:)
  integer,allocatable :: N_k1(:)

  integer,allocatable :: k2_to_iorb(:,:,:)
  integer,allocatable :: N_k2(:)
  integer,allocatable :: icheck_k2(:)

  integer,allocatable :: Q_NRps(:,:)
  real(8),allocatable :: Q_Rps(:,:)




CONTAINS

!------------------------------------------
  SUBROUTINE allocateKtoK( k1max,k2max,nki,Rrefmax,Lrefmax )
    implicit none
    integer,intent(IN) :: k1max,k2max,nki,Rrefmax,Lrefmax
    
    if ( allocated(k1_to_k2) ) then
      call deallocateKtoK
    end if

    allocate( k1_to_k2(k1max,nki)     ) ; k1_to_k2(:,:)=0
    allocate( k1_to_k3(k1max,nki)     ) ; k1_to_k3(:,:)=0
    allocate( k1_to_iorb(2,k1max,nki) ) ; k1_to_iorb(:,:,:)=0
    allocate( N_k1(nki)               ) ; N_k1(:)=0

    allocate( k2_to_iorb(2,k2max,nki)   ) ; k2_to_iorb(:,:,:)=0
    allocate( N_k2(nki)               ) ; N_k2(:)=0

    allocate( icheck_k2(nki)          ) ; icheck_k2(:)=0

    allocate( qqc(Rrefmax,Rrefmax,Lrefmax,nki) ) ; qqc(:,:,:,:)=0.d0

    return
  END SUBROUTINE allocateKtoK

!------------------------------------------
  SUBROUTINE deallocateKtoK
    implicit none
    
    deallocate( k1_to_k2 )
    deallocate( k1_to_k3 )
    deallocate( k1_to_iorb )
    deallocate( N_k1 )

    deallocate( k2_to_iorb )
    deallocate( N_k2 )

    deallocate( icheck_k2 )

    deallocate( qqc )

    return
  END SUBROUTINE deallocateKtoK

!------------------------------------------
  SUBROUTINE checkKtoK

  END SUBROUTINE checkKtoK

!------------------------------------------

  SUBROUTINE allocateQRps( k2max,nki )
    implicit none
    integer,intent(IN) :: k2max,nki

    if ( allocated(Q_NRps) ) then
      call deallocateQRps
    end if

    allocate( Q_NRps(k2max,nki) ) ; Q_NRps(:,:)=0
    allocate( Q_Rps(k2max,nki)  ) ; Q_Rps(:,:)=0.d0

    return
  END SUBROUTINE allocateQRps

!------------------------------------------
  SUBROUTINE deallocateQRps
    implicit none

    deallocate( Q_NRps )
    deallocate( Q_Rps  )

    return
  END SUBROUTINE deallocateQRps

END MODULE VarPSMemberG
