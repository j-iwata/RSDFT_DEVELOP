MODULE ewald_variables

  implicit none

  PRIVATE
  PUBLIC :: eta, mg, mr, LG, LR, ipair, mpair, rrcut, ecut

  real(8) :: eta
  integer :: mg,mr
  integer,allocatable :: LG(:,:),LR(:,:)
  integer :: mpair
  integer,allocatable :: ipair(:,:)
  real(8) :: rrcut, ecut

END MODULE ewald_variables
