MODULE BasicFunctions
   USE BasicTypeFactory
   implicit none

   PRIVATE
   PUBLIC :: allocateArray

CONTAINS

   SUBROUTINE allocateArray( gs )
     implicit none
     type( GSArray ),optional :: gs
     if ( present(gs) ) then
             allocate( gs%val(gs%g_prange%head:gs%g_prange%tail, gs%s_srange%head:gs%s_srange%tail) )
      endif
   END SUBROUTINE allocateArray
END MODULE BasicFunctions
