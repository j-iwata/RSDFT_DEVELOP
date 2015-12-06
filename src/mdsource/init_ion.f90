!---------------------------------
! Convert to Cartesian coordinate
!---------------------------------
SUBROUTINE init_ion
  use aa_module
  use atom_module
  use cpmd_variables, only: Rion
  implicit none
  integer :: a,i
  do i=1,3
     do a=1,Natom
        if ( aa_atom(i,a)>=1.d0 ) then
           aa_atom(i,a)=aa_atom(i,a)-1.d0
        else if ( aa_atom(i,a)<0.d0 ) then
           aa_atom(i,a)=aa_atom(i,a)+1.d0
        end if
     end do
  end do
  do a=1,Natom
     Rion(1,a)=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
     Rion(2,a)=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
     Rion(3,a)=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
  end do
  return
END SUBROUTINE init_ion
