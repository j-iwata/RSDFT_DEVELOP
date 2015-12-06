MODULE localpot2_variables

  implicit none

  PRIVATE
  PUBLIC :: fecut_loc, Ngrid_dense, Igrid_dense, dV_dense &
       ,vion_nl, vh_nl, vxc_nl, vloc_dense, vloc_dense_old, rho_nl &
       ,Ndens_loc, nitp_0, nitp_1, Clag1, Clag2, Clag3

  real(8),allocatable :: vion_nl(:,:,:)
  real(8),allocatable :: vh_nl(:,:,:)
  real(8),allocatable :: vxc_nl(:,:,:)
  real(8),allocatable :: rho_nl(:,:,:)
  real(8),allocatable :: vloc_dense(:,:,:)
  real(8),allocatable :: vloc_dense_old(:,:,:)

  integer :: Ngrid_dense(0:3),Ndens_loc,nitp_0,nitp_1
  integer :: Igrid_dense(2,3)
  real(8),allocatable :: Clag1(:,:),Clag2(:,:),Clag3(:,:)
  real(8) :: dV_dense

  real(8) :: fecut_loc

END MODULE localpot2_variables
