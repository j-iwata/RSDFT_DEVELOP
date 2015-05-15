MODULE BasicFunctions

  use BasicTypeFactory

  implicit none

  PRIVATE
  PUBLIC :: allocateArray, allocateGSArray

CONTAINS

  SUBROUTINE allocateArray( gs )
    implicit none
    type( GSArray ),optional :: gs
    if ( present(gs) ) then
       allocate( gs%val(gs%g_prange%head:gs%g_prange%tail, gs%s_srange%head:gs%s_srange%tail) )
    endif
  END SUBROUTINE allocateArray

  SUBROUTINE allocateGSArray( gs )
    implicit none
    integer :: m1,m2,n1,n2
    type( GSArray ),optional :: gs
    m1=gs%g_range%head
    m2=gs%g_range%tail
    n1=gs%s_range%head
    n2=gs%s_range%tail
    allocate( gs%val(m1:m2,n1:n2) ) ; gs%val=0.0d0
  END SUBROUTINE allocateGSArray

END MODULE BasicFunctions
