MODULE VarParaPSnonLocG
  implicit none

  integer :: Mqr
  integer :: MMJJ_Q,MMJJ_t_Q
  integer :: MAXMJJ

  real(8),allocatable :: QRij(:,:)
  integer,allocatable :: JJ_MAP_Q(:,:,:)
  integer,allocatable :: MJJ_MAP_Q(:)
  integer,allocatable :: MMJ_Q(:)
  integer,allocatable :: nl_rank_map_Q(:)

  integer,allocatable :: qr_nsend(:)
  integer,allocatable :: sendmap_Q(:,:),recvmap_Q(:,:)

CONTAINS

END MODULE VarParaPSnonLocG
