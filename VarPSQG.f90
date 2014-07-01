MODULE VarPSQG
  
  implicit none

  integer,allocatable :: k1_to_k2(:,:)
  integer,allocatable :: k1_to_k3(:,:)
  integer,allocatable :: k1_to_iorb(:,:,:)
  integer,allocatable :: N_k1(:)

  real(8),allocatable :: qqc(:,:,:,:)

CONTAINS

!------------------------------------------
  SUBROUTINE allocateKtoK( k1max,nki,Rrefmax,Lrefmax )
    implicit none
    
    if ( allocated(k1_to_k2) ) then
        deallocate( k1_to_k2 )
        deallocate( k1_to_k3 )
        deallocate( k1_to_iorb )
        deallocate( N_k1 )
        deallocate( qcc )
    end if

    allocate( k1_to_k2(k1max,nki)     ) ; k1_to_k2(:,:)=0
    allocate( k1_to_k3(k1max,nki)     ) ; k1_to_k3(:,:)=0
    allocate( k1_to_iorb(2,k1max,nki) ) ; k1_to_iorb(:,:,:)=0
    allocate( N_k1(nki)               ) ; N_k1(:)=0
    allocate( qcc(Rrefmax,Rrefmax,Lrefmax,nki) ) ; qcc(:,:,:,:)=0.d0

    return
  END SUBROUTINE allocateKtoK

!------------------------------------------
  SUBROUTINE deallocateKtoK
    implicit none
    
    deallocate( k1_to_k2 )
    deallocate( k1_to_k3 )
    deallocate( k1_to_iorb )
    deallocate( N_k1 )
    deallocate( qcc )

    return
  END SUBROUTINE deallocateKtoK

!------------------------------------------
  SUBROUTINE checkKtoK

  END SUBROUTINE checkKtoK

!------------------------------------------


END MODULE VarPSQG
