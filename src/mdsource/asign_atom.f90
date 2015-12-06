!-----------------------------------------------------------------------
!     Asign atomic numer
!-----------------------------------------------------------------------
subroutine asign_atom
  use pseudopot_module, only: file_ps
  use atom_module, only: ki_atom,Natom
  use cpmd_variables, only: iatom,batm
  implicit none 
  integer :: i,k
  character(2) :: bion
  do i=1,Natom
     bion=file_ps(ki_atom(i))
     do k=1,36
        if ( bion == batm(k) ) then 
           iatom(i)=k
           exit 
        endif
     enddo
     if ( k > 36 ) then
        write(*,*) "atomic number can not be asigned",i
        write(*,*) "file_ps, bion= ",file_ps(ki_atom(i)),bion
        stop
     end if
  enddo
  return
end subroutine asign_atom
