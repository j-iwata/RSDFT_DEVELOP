MODULE virtualH_module

  use atom_module, only: zn_atom
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: virtualH

  real(8) :: Zvh=0.0d0

CONTAINS


  SUBROUTINE read_virtualH
    implicit none
    call IOTools_readReal8Keyword( "VIRTUALH", Zvh )
  END SUBROUTINE read_virtualH


  SUBROUTINE virtualH( Zion, Vlocal )
    implicit none
    real(8),intent(INOUT) :: Zion(:), Vlocal(:,:)
    integer :: i
    call read_virtualH
    do i=1,size(Zion)
       if ( zn_atom(i) == 1 ) then
          if ( Zvh /= 0.0d0 ) then
             Vlocal(:,i) = (Zvh/Zion(i))*Vlocal(:,i)
             Zion(i) = Zvh
          end if
       end if
    end do
  END SUBROUTINE virtualH


END MODULE virtualH_module
