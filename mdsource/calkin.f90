!-----------------------------------------------------------------------
!     Kinetic energy of nuclei
!-----------------------------------------------------------------------
subroutine calkin(kine)
  use atom_module, only: Natom
  use cpmd_variables, only: pmass,AMU,Velocity,iatom
  implicit none
  real(8),intent(inout) :: kine
  integer :: i
  real(8) :: pm
  kine=0.0d0
  do i=1,Natom
     pm=pmass(iatom(i))*AMU
     kine=kine+(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
  enddo
  kine=kine*0.5d0
  return
end subroutine calkin
