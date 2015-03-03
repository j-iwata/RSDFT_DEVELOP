MODULE xc_variables

  use grid_module, only: grid

  implicit none

  PRIVATE
  PUBLIC :: xc, init_type_xc

  type xc
     real(8) :: Exc
     real(8) :: Ex
     real(8) :: Ec
     real(8),allocatable :: Vxc(:,:)
     integer :: mm,m0,m1
     integer :: nn,n0,n1
  end type xc

CONTAINS

  SUBROUTINE init_type_xc( MSP_0, MSP_1, r, f )
    implicit none
    integer,intent(IN) :: MSP_0,MSP_1
    type(grid),intent(IN)  :: r
    type(xc),intent(OUT) :: f
    f%m0 = r%SubGrid(1,0)
    f%m1 = r%SubGrid(2,0)
    f%n0 = MSP_0
    f%n1 = MSP_1
    if ( .not.allocated(f%Vxc) ) then
       allocate( f%Vxc(f%m0:f%m1,f%n0:f%n1) )
    end if
    f%Vxc=0.0d0
  END SUBROUTINE init_type_xc

END MODULE xc_variables
