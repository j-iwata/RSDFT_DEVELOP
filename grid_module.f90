MODULE grid_module

  use rgrid_variables

  implicit none

  PRIVATE
  PUBLIC :: grid, init_grid

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

END MODULE grid_module
