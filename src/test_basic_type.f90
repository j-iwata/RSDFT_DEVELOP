PROGRAM Test_BasicType
  USE TestAssert
  USE BasicTypeFactory
  implicit none
  type( ArrayRange1D ) :: a,b
  type( ArrayRange3C ) :: c
  type( ArrayRange3D ) :: d
  call assert_init
  a%head = 1
  a%tail = 4
  call getSize1D( a )
  call assert( a%size==4, 'get ArrayRange1D size' )
  c%r(1)%head = 1
  c%r(1)%tail = 3
  c%r(2)%head = 1
  c%r(2)%tail = 3
  c%r(3)%head = 1
  c%r(3)%tail = 3
  call getSize3C( c )
  call assert( c%r(1)%size==3, 'get ArrayRange3C size element' )
  call assert( c%size==27, 'get ArrayRange3C size' )
  d%x%head = 1
  d%x%tail = 3
  d%y%head = 1
  d%y%tail = 3
  d%z%head = 1
  d%z%tail = 3
  call getSize3D( d )
  call assert( d%x%size==3, 'get ArrayRange3D size element' )
  call assert( d%size==27, 'get ArrayRange3D size' )

  call assert_finalize('BasicType test')
END PROGRAM Test_BasicType
