MODULE BasicTypeFactory

  implicit none

  PRIVATE
  PUBLIC :: get_range_size

  type,PUBLIC :: ArrayRange1D
    sequence
    integer :: head
    integer :: tail
    integer :: size
    integer :: head_global
    integer :: tail_global
    integer :: size_global
  end type ArrayRange1D

  type,PUBLIC :: ArrayRange3D
    sequence
    type ( ArrayRange1D ) :: x
    type ( ArrayRange1D ) :: y
    type ( ArrayRange1D ) :: z
    integer :: size
  end type ArrayRange3D

  type,PUBLIC :: ArrayRange1D_v2
    sequence
    integer :: head
    integer :: tail
    integer :: size
 end type ArrayRange1D_v2

  type,PUBLIC :: ArrayRange3D_v2
    sequence
    type ( ArrayRange1D ) :: r(3)
    integer :: size
  end type ArrayRange3D_v2

  type,PUBLIC :: pArrayRange1D
     sequence
     type( ArrayRange1D ) :: globl
     type( ArrayRange1D ) :: local
     type( ArrayRange1D ) :: alloc
  end type pArrayRange1D

  type,PUBLIC :: pArrayRange3D
     sequence
     type( pArrayRange1D ) :: r(3)
  end type pArrayRange3D

  type,PUBLIC :: Array1D
    sequence
#ifdef REAL_VER
    double precision,allocatable :: array(:)
#elif defined COMPLEX_VER
    complex(kind(0d0)),allocatable :: array(:)
#endif
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type Array1D

  type,PUBLIC :: Array3D
    sequence
#ifdef REAL_VER
    double precision,allocatable :: array(:,:,:)
#elif defined COMPLEX_VER
    complex(kind(0d0)),allocatable :: array(:,:,:)
#endif
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type Array3D

  type,PUBLIC :: rArray1D
    sequence
    double precision,allocatable :: array(:)
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type rArray1D

  type,PUBLIC :: rArray3D
    sequence
    double precision,allocatable :: array(:,:,:)
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type rArray3D

  type,PUBLIC :: cArray1D
    sequence
    complex(kind(0d0)),allocatable :: array(:)
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type cArray1D

  type,PUBLIC :: cArray3D
    sequence
    complex(kind(0d0)),allocatable :: array(:,:,:)
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type cArray3D

  type,PUBLIC :: GSArray
     sequence
     type( ArrayRange1D ) :: g_srange, g_prange
     type( ArrayRange1D ) :: s_srange, s_prange
     type( ArrayRange1D ) :: g_range
     type( ArrayRange1D ) :: s_range
     real(8),allocatable :: val(:,:)
  end type GSArray

  type,PUBLIC :: GSArray_v2
     sequence
     type( pArrayRange1D ) :: g_range
     type( pArrayRange1D ) :: s_range
     real(8),allocatable :: val(:,:)
  end type GSArray_v2

  type,PUBLIC :: GBKSArray
     sequence
     type( ArrayRange1D ) :: g_srange, g_prange
     type( ArrayRange1D ) :: b_srange, b_prange
     type( ArrayRange1D ) :: k_srange, k_prange
     type( ArrayRange1D ) :: s_srange, s_prange
#ifdef REAL_VER
     real(8),allocatable :: val(:,:,:,:)
#else 
     complex(kind(0d0)),allocatable :: val(:,:,:,:)
#endif
  end type GBKSArray

  type,PUBLIC :: rGBKSArray
     sequence
     type( ArrayRange1D ) :: g_srange, g_prange
     type( ArrayRange1D ) :: b_srange, b_prange
     type( ArrayRange1D ) :: k_srange, k_prange
     type( ArrayRange1D ) :: s_srange, s_prange
     real(8),allocatable :: val(:,:,:,:)
  end type rGBKSArray

  type,PUBLIC :: cGBKSArray
     sequence
     type( ArrayRange1D ) :: g_srange, g_prange
     type( ArrayRange1D ) :: b_srange, b_prange
     type( ArrayRange1D ) :: k_srange, k_prange
     type( ArrayRange1D ) :: s_srange, s_prange
     complex(kind(0d0)),allocatable :: val(:,:,:,:)
  end type cGBKSArray

CONTAINS

  SUBROUTINE get_range_size( a )
    implicit none
    type( ArrayRange1D ) :: a
    a%size = a%tail - a%size + 1
  END SUBROUTINE get_range_size

END MODULE BasicTypeFactory
