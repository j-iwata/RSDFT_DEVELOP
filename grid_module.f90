MODULE grid_module

  use rgrid_variables
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: grid, init_grid
  PUBLIC :: get_map_3d_to_1d

  type grid
     integer :: NumGrid(0:3)
     integer :: SubGrid(2,0:3)
     real(8) :: SizeGrid(3)
     real(8) :: dV
  end type grid

CONTAINS

  SUBROUTINE init_grid( g )
    implicit none
    type(grid) :: g
    g%NumGrid(0:3)  = Ngrid(0:3)
    g%SizeGrid(1:3) = Hgrid(1:3)
    g%SubGrid(1:2,0:3) = Igrid(1:2,0:3)
    g%dV = dV
  END SUBROUTINE init_grid

  SUBROUTINE get_map_3d_to_1d( LLL )
     implicit none
     integer,allocatable,intent(OUT) :: LLL(:,:,:)
     integer :: i,i1,i2,i3,ierr
     allocate( LLL(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) ) ; LLL=0
     i=Igrid(1,0)-1
     do i3=Igrid(1,3),Igrid(2,3)
     do i2=Igrid(1,2),Igrid(2,2)
     do i1=Igrid(1,1),Igrid(2,1)
        i=i+1
        LLL(i1,i2,i3)=i
     end do
     end do
     end do
     call MPI_ALLREDUCE( MPI_IN_PLACE, LLL, size(LLL), MPI_INTEGER, MPI_SUM, comm_grid, ierr )
  END SUBROUTINE get_map_3d_to_1d


END MODULE grid_module
