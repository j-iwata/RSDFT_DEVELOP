MODULE eion_mol_module

  PRIVATE
  PUBLIC :: calc_eion_mol

CONTAINS

  SUBROUTINE calc_eion_mol(Eion)
    use pseudopot_module, only: Zps
    use atom_module, only: Natom,ki_atom,aa_atom
    real(8),intent(OUT) :: Eion
    real(8) :: r
    integer :: i,j,ik,jk
    Eion=0.d0
    do i=1,Natom
       ik=ki_atom(i)
       do j=1,i-1
          jk=ki_atom(j)
          r=sqrt( ( aa_atom(1,i) - aa_atom(1,j) )**2 &
                 +( aa_atom(2,i) - aa_atom(2,j) )**2 &
                 +( aa_atom(3,i) - aa_atom(3,j) )**2 )
          Eion = Eion + Zps(ik)*Zps(jk)/r
       end do
    end do
    Eion = 0.5d0*Eion
  END SUBROUTINE calc_eion_mol

END MODULE eion_mol_module
