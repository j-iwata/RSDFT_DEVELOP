MODULE localpot2_variables

  implicit none

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

CONTAINS

  SUBROUTINE init_localpot2(n1,n2,n3)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    integer :: m1,m2,m3,m1_0,m1_1,m2_0,m2_1,m3_0,m3_1
    Ngrid_dense(1) = n1*Ndens_loc
    Ngrid_dense(2) = n2*Ndens_loc
    Ngrid_dense(3) = n3*Ndens_loc
    Ngrid_dense(0) = Ngrid_dense(1)*Ngrid_dense(2)*Ngrid_dense(3)
    m1_0=Igrid_dense(1,1)
    m1_1=Igrid_dense(2,1)
    m2_0=Igrid_dense(1,2)
    m2_1=Igrid_dense(2,2)
    m3_0=Igrid_dense(1,3)
    m3_1=Igrid_dense(2,3)
    allocate( vion_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ; vion_nl=0.0d0
    allocate(   vh_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;   vh_nl=0.0d0
    allocate(  vxc_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;  vxc_nl=0.0d0
    allocate(  rho_nl(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ;  rho_nl=0.0d0
    allocate( vloc_dense(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) ) ; vloc_dense=0.0d0
    allocate( vloc_dense_old(m1_0:m1_1,m2_0:m2_1,m3_0:m3_1) )
    vloc_dense_old=0.0d0
  END SUBROUTINE init_localpot2

END MODULE localpot2_variables
