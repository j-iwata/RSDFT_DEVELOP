MODULE localpot2_variables

  implicit none

  real(8),allocatable :: vion_nl(:,:,:)
  real(8),allocatable :: vh_nl(:,:,:)
  real(8),allocatable :: vxc_nl(:,:,:)
  real(8),allocatable :: rho_nl(:,:,:)

  integer :: Ngrid_dense(3),Ndens_loc,nitp_0,nitp_1
  real(8),allocatable :: Clag1(:,:),Clag2(:,:),Clag3(:,:)

CONTAINS

  SUBROUTINE init_localpot2(n1,n2,n3)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    integer :: m1,m2,m3
    Ngrid_dense(1) = n1*Ndens_loc
    Ngrid_dense(2) = n2*Ndens_loc
    Ngrid_dense(3) = n3*Ndens_loc
    m1=Ngrid_dense(1)
    m2=Ngrid_dense(2)
    m3=Ngrid_dense(3)
    allocate( vion_nl(0:m1-1,0:m2-1,0:m3-1) ) ; vion_nl=0.0d0
    allocate( vh_nl(0:m1-1,0:m2-1,0:m3-1)   ) ; vh_nl=0.0d0
    allocate( vxc_nl(0:m1-1,0:m2-1,0:m3-1)  ) ; vxc_nl=0.0d0
    allocate( rho_nl(0:m1-1,0:m2-1,0:m3-1)  ) ; rho_nl=0.0d0
  END SUBROUTINE init_localpot2

END MODULE localpot2_variables
