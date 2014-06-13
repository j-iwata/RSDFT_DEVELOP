MODULE VarPSMember
  implicit none

  integer :: Nelement_PP
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

  integer :: max_psgrd=0,max_psorb=0,max_ngauss=0


END MODULE VarPSMember
