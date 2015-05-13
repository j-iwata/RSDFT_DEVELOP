MODULE BasicTypeFactory
  implicit none
  type ArrayRange1D
    sequence
    integer :: head
    integer :: tail
    integer :: size
  end type

  type ArrayRange3D
    sequence
    type ( ArrayRange1D ) :: x
    type ( ArrayRange1D ) :: y
    type ( ArrayRange1D ) :: z
    integer :: size
  end type
  
  type Array1D
    sequence
#ifdef REAL_VER
    double precision,allocatable :: array(:)
#elif defined COMPLEX_VER
    complex(kind(0d0)),allocatable :: array(:)
#endif
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type

  type Array3D
    sequence
#ifdef REAL_VER
    double precision,allocatable :: array(:,:,:)
#elif defined COMPLEX_VER
    complex(kind(0d0)),allocatable :: array(:,:,:)
#endif
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type

  type rArray1D
    sequence
    double precision,allocatable :: array(:)
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type

  type rArray3D
    sequence
    double precision,allocatable :: array(:,:,:)
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type

  type cArray1D
    sequence
    complex(kind(0d0)),allocatable :: array(:)
    type ( ArrayRange1D ) :: s_range
    type ( ArrayRange1D ) :: p_range
  end type

  type cArray3D
    sequence
    complex(kind(0d0)),allocatable :: array(:,:,:)
    type ( ArrayRange3D ) :: s_range
    type ( ArrayRange3D ) :: p_range
  end type

  type GSArray
     sequence
     type( ArrayRange1D ) :: g_srange, g_prange
     type( ArrayRange1D ) :: s_srange, s_prange
     real(8),allocatable :: val(:,:)
  end type

END MODULE BasicTypeFactory
