PROGRAM Test_BasicType
  use BasicTypeFactory
  implicit none
  type( Array1D ) :: a,b
  complex(kind(0d0)),parameter :: z0 = ( 0.d0, 0.d0 )
  integer :: i
  type( ArrayRange2D ) :: xy
  a%s_range%head = 1
  a%s_range%tail = 4
  call getSize1D( a%s_range )
  allocate( a%val(a%s_range%size ) ) ; a%val = z0
  do i = a%s_range%head, a%s_range%tail
    a%val(i) = i
  enddo
  do i = a%s_range%head, a%s_range%tail
    write(*,*) a%val(i)
  enddo
  b%s_range = a%s_range
  write(*,*) b%s_range%head, b%s_range%tail
!  do i = b%s_range%head, b%s_range%tail
!    write(*,*) b%val(i)
!  enddo
  xy%r(1)%head = 1
  xy%r(1)%tail = 4
  xy%r(2)%head = 1
  xy%r(2)%tail = 4

END PROGRAM Test_BasicType
