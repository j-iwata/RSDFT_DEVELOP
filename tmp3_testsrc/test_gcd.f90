PROGRAM test

  implicit none

  integer :: m0,n0,m,n,mtmp,loop
  integer,parameter :: maxloop=100

  write(*,*) "m,n="
  read(*,*) m0,n0

  if ( m0 >= n0 ) then
     m=m0
     n=n0
  else
     m=n0
     n=m0
  end if

  do loop=1,maxloop

     if ( n == 0 ) exit

     mtmp = n
     n = mod(m,n)
     m = mtmp

  end do

  write(*,*) m, "   loop",loop
  write(*,*) m0/m,n0/m

END PROGRAM test
