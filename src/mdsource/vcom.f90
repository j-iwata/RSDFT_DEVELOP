!-----------------------------------------------------------------------
!     Make center of mass motion off
!-----------------------------------------------------------------------

SUBROUTINE vcom( Velocity )

  use atom_module, only: Natom,ki_atom,zn_atom
  use cpmd_variables, only: pmass,AMU
  use calc_kine_temp_module, only: calc_kine

  implicit none
  real(8),intent(INOUT) :: Velocity(3,Natom)
  real(8) rescale,pm,kine0,kine,MG,VG(3)
  integer i

  call calc_kine( Velocity, kine0 )

  MG=0.0d0
  VG=0.0d0

  do i=1,Natom
     pm = pmass( zn_atom(ki_atom(i)) )*AMU
     MG = MG + pm
     VG(1) = VG(1) + Velocity(1,i)*pm
     VG(2) = VG(2) + Velocity(2,i)*pm
     VG(3) = VG(3) + Velocity(3,i)*pm
  end do

  VG(1:3)=VG(1:3)/MG
  do i=1,Natom
     Velocity(1,i) = Velocity(1,i) - VG(1)
     Velocity(2,i) = Velocity(2,i) - VG(2)
     Velocity(3,i) = Velocity(3,i) - VG(3)
  end do

  call calc_kine( Velocity, kine )

  rescale = sqrt( kine0/kine )

  Velocity(:,:) = Velocity(:,:)*rescale

END SUBROUTINE vcom
