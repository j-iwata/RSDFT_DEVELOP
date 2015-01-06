!
! from NUMERICAL RECIPES
!
SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  implicit none
  integer :: n
  real(8) :: x(n),y(n),yp1,ypn,y2(n)
  integer,parameter :: nmax=50000
  integer :: i,k
  real(8) :: p,qn,sig,un,u(nmax)
  integer,save :: unit=300

  if ( yp1 > 1.d30 ) then
     y2(1) = 0.d0
      u(1) = 0.d0
  else
     y2(1) = -0.5d0
      u(1) = ( 3.d0/(x(2)-x(1)) )*( (y(2)-y(1))/(x(2)-x(1)) - yp1 )
  end if

  do i=2,n-1
     sig = ( x(i)-x(i-1) )/( x(i+1)-x(i-1) )
     p = sig*y2(i-1)+2.d0
     y2(i) = ( sig-1.d0 )/p
     u(i) = ( 6.d0*( (y(i+1)-y(i))/(x(i+1)-x(i)) &
                    -(y(i)-y(i-1))/(x(i)-x(i-1)) )/(x(i+1)-x(i-1))-sig*u(i-1) )/p
  end do

  if ( ypn > 1.d30 ) then
     qn = 0.d0
     un = 0.d0
  else
     qn = 0.5d0
     un = ( 3.d0/(x(n)-x(n-1)) )*( ypn - (y(n)-y(n-1))/(x(n)-x(n-1)) )
  end if
  y2(n) = ( un - qn*u(n-1) )/( qn*y2(n-1) + 1.d0 )
  do k=n-1,1,-1
     y2(k) = y2(k)*y2(k+1) + u(k)
  end do

END SUBROUTINE spline

SUBROUTINE splint(xa,ya,y2a,n,x,y)
  implicit none
  integer :: n
  real(8) :: x,y,xa(n),y2a(n),ya(n)
  integer :: k,khi,klo
  real(8) :: a,b,h
  klo=1
  khi=n
1 if ( khi-klo > 1 ) then
     k=(khi+klo)/2
     if ( xa(k) > x ) then
        khi=k
     else
        klo=k
     end if
     goto 1
  end if
  h=xa(khi)-xa(klo)
  if ( h == 0.d0 ) stop "splint"
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+( (a**3-a)*y2a(klo)+(b**3-b)*y2a(khi) )*(h**2)/6.d0
END SUBROUTINE splint
