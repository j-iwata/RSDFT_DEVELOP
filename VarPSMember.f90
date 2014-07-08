MODULE VarPSMember
  implicit none

  integer :: Nelement_PP
  integer :: Nelement_local
  integer,allocatable :: ippform(:)
  character(30),allocatable :: file_ps(:)
  real(8),allocatable :: rad(:,:),rab(:,:)  ! rad1(:,:)
  real(8),allocatable :: vql(:,:),viod(:,:,:)   ! dvql(:,:),dviod(:,:,:)
  real(8),allocatable :: cdc(:,:),cdd(:,:)
  real(8),allocatable :: anorm(:,:)

  real(8),allocatable :: Rps(:,:)   ! Rps0(:,:),Rps1(:,:)
  integer,allocatable :: lo(:,:),inorm(:,:),NRps(:,:)   ! NRps0(:,:),NRps1(:,:)

  integer,allocatable :: Mr(:),norb(:)
  real(8),allocatable :: parloc(:,:)
  real(8),allocatable :: Zps(:)
  real(8),allocatable :: cdd_coef(:,:,:)

  integer :: Rrefmax,Lrefmax

  integer,allocatable :: NRps0(:,:)
  real(8),allocatable :: Rps0(:,:)

  integer :: max_psgrd=0,max_psorb=0,max_ngauss=0

  real(8),allocatable :: psi(:,:,:,:)
  real(8),allocatable :: phi(:,:,:,:)
  real(8),allocatable :: beta(:,:,:,:)


CONTAINS
!------------------------------------------
  SUBROUTINE allocateRps
    implicit none
    integer :: m

    if (allocated( NRps0 )) then
      call deallocateRps
    end if
    m=maxval( norb )
    allocate( NRps0(m,Nelement_local) ) ; NRps0=0
    allocate( Rps0(m,Nelement_local)  ) ; Rps0=0.d0
    NRps0(:,:)=NRps(:,:)
    Rps0(:,:)=Rps(:,:)

    return
  END SUBROUTINE allocateRps

!------------------------------------------
  SUBROUTINE deallocateRps
    implicit none

    deallocate( NRps0 )
    deallocate( Rps0  )

    return
  END SUBROUTINE deallocateRps

!------------------------------------------
    SUBROUTINE allocatePSI
        implicit none
        
        if allocated( psi ) then
            call deallocatePSI
        end if


    END SUBROUTINE allocatePSI

!------------------------------------------
    SUBROUTINE deallocatePSI
        implicit none
        
        deallocate( psi  )
        deallocate( phi  )
        deallocate( beta )
        
        return
    END SUBROUTINE deallocatePSI

END MODULE VarPSMember
