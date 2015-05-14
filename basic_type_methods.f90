MODULE BasicTypeMethods
  USE BasicTypeFactory
  implicit none

  PRIVATE
  PUBLIC :: allocateGS
  PUBLIC :: allocateGBKS
  PUBLIC :: allocaterGBKS
  PUBLIC :: allocatecGBKS

#ifdef REAL_VER
  double precision,parameter :: zero = 0.d0
#elif defined COMPLEX_VER
  complex(kind(0d0)),parameter :: zero = (0.d0,0.d0)
#endif
  complex(kind(0d0)),parameter :: z0 = (0.d0,0.d0)

CONTAINS
  SUBROUTINE allocateGS( gs )
    implicit none
    type( GSArray    ) ::  gs
    allocate( gs%val(gs%g_prange%head:gs%g_prange%tail, gs%s_srange%head:gs%s_srange%tail) )
    gs%val = 0.d0
  END SUBROUTINE allocateGS

  SUBROUTINE allocateGBKS( gbks )
    implicit none
    type( GBKSArray  ) ::  gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = zero
  END SUBROUTINE allocateGBKS

  SUBROUTINE allocaterGBKS( gbks )
    implicit none
    type( rGBKSArray ) :: gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = 0.d0
  END SUBROUTINE allocaterGBKS

  SUBROUTINE allocatecGBKS( gbks )
    implicit none
    type( cGBKSArray ) :: gbks
    allocate( gbks%val(gbks%g_prange%head:gbks%g_prange%tail, gbks%b_prange%head:gbks%b_prange%tail, gbks%k_prange%head:gbks%k_prange%tail, gbks%s_prange%head:gbks%s_prange%tail) )
    gbks%val = z0
  END SUBROUTINE allocatecGBKS
END MODULE BasicTypeMethods
