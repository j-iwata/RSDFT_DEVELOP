!------------------------------------------------------------!
! Clebsh-Gordan coefficient (for any angular momentum)       !
!------------------------------------------------------------!
! This program is a copy from                                !
!   S. Okabe                                                 !
!   Simyulation butsurigaku 5, Ryousiron - Undou to Houhou - !
!   pp348-350, Kindai-kagaku sya (1992)                      !
!------------------------------------------------------------!
FUNCTION clgr( xj1, xm1, xj2, xm2, xj3, xm3 )

  implicit none

  real(8) :: clgr
  real(8),intent(IN) :: xj1,xm1,xj2,xm2,xj3,xm3
  integer,parameter :: maxj=199, maxj1=maxj+1, maxj2=maxj+2, maxc=10302
  real(8) :: cg(maxj1)
  real(8),allocatable :: cga(:)
  real(8),save :: c(maxc)
  integer,save :: nn(maxj2)
  integer :: ifirst=0

  integer :: mj1,n,m,jmn,i,m3,m2,m1,j3,j2,j1,j123,n1,mk
  integer :: k,i1,l,i2,i3,n0,kn,ne,jm1,jmp
  real(8) :: x,xjj,yj1,z2,z1,xx,xk,ym2,ym1,yj2,yjj2,yjj1,cj

  if ( ifirst == 0 ) then

     mj1 = maxj + 1

     nn(1)=0

     do n=1,mj1
        nn(n+1) = nn(n) + (n+1)/2
     end do

     c(1) = 1.0d0

     i=1

     do n=1,mj1
        n0=(n/2)*2-n
        kn=n/2+1
        do k=1,kn
           i=i+1
           i2=nn(n)+k
           i3=i2-1
           c(i)=1.0d0
           if ( k /= 1 ) c(i)=c(i2)+c(i3)
           if ( n0 == 0 .and. k == kn ) c(i)=c(i3)+c(i3)
        end do
     end do

     ifirst = 1

  end if

  clgr = 0.0d0

  if ( nint(xj1+xj2+xj3) >= maxj ) then
     write(*,*) "The values of the angluar momenta are too large !"
     stop "stop@clgr(1)"
  end if

  allocate( cga(nint(2*xj1+1)) ) ; cga=0.0d0

  jm1 = nint( xm1+xj1+1 )

  if ( xj2 <= xj1 ) then

     ne=0
     yj1 = xj1
     yj2 = xj2

  else

     ne=1
     yj1 = xj2
     yj2 = xj1

  end if

  yjj1 = yj1*( yj1 + 1.0d0 )
  yjj2 = yj2*( yj2 + 1.0d0 )
  xjj = xj3*( xj3 + 1.0d0 ) - yjj1 - yjj2

  j1 = nint( yj1 + yj1 )
  j2 = nint( yj2 + yj2 )
  j3 = nint( xj3 + xj3 )
  m1 = nint( xm1 + xm1 )
  m2 = nint( xm2 + xm2 )
  m3 = nint( xm3 + xm3 )

  if ( j3 < abs(m3) .or. j3 > j1+j2 .or. j3 < abs(j1-j2) .or. &
       mod(j1+j2+j3,2) == 1 .or. mod(m1+m2+m3,2) == 1 .or. &
       mod(j3+m3,2) == 1 .or. mod(nint(xj1+xj1)+m1,2) == 1 ) then

     kn = nint( xj1+xj1+1 )
     do i=1,kn
        cga(i)=0.0d0
     end do
     deallocate( cga )
     return

  end if

  j123 = ( j1 + j2 + j3 )/2
  jmp = ( j3 + m3 )/2
  jmn = ( j3 - m3 )/2

  if ( j3 == 0 ) then
     x = 1.0d0/sqrt(yj1+yj1+1.0d0)
     kn = j1 + 1
     do m=1,kn
        k=(kn-m)/2
        n=k+k-kn+m
        cga(m)=(-1.0d0)**n * x
     end do
     clgr = cga(jm1)
     deallocate( cga )
     return
  end if

  n=j123+1
  k=j2+1
  i=nn(n+1)+k+1
  l=n-k-k

  if ( l < 0 ) i=i+l
  cj = (xj3+xj3+1.0d0)/(yj2+yj2+1.0d0)/c(i)

  kn=j1+1
  do k=1,kn
     cg(k)=0.0d0
  end do
  if ( m3 > j1-j2 ) goto 6

  n=jmn
  k=j123-j1
  i1=nn(n+1)+k+1
  l=n-k-k
  if ( l < 0 ) i1=i1+l

  n=j123-jmn
  k=jmp
  i2=nn(n+1)+k+1
  l=n-k-k
  if ( l < 0 ) i2=i2+l

  m=n+1
  cg(m)=sqrt( cj*c(i1)*c(i2) )

  kn=j2
  kn=min(n,kn)
  if ( m3 == 0 ) kn=j2/2
  if ( kn == 0 ) goto 7

  ym1 = xm3 + yj2
  ym2 = -yj2

  cg(m-1)=cg(m)*(xjj-2.0d0*ym1*ym2)/sqrt((yjj1-ym1*ym1+ym1)*(yjj2-ym2*ym2-ym2))

  if ( kn == 1 ) goto 7

  do k=2,kn
     xk=k-1
     ym2=xk-yj2
     ym1=xm3-ym2
     xx=ym1*ym2
     z1=yjj1-ym1*ym1
     z2=yjj2-ym2*ym2
     mk=m-k
     cg(mk) = ( cg(mk+1)*(xjj-xx-xx)-cg(mk+2)*sqrt((z1-ym1)*(z2+ym2)) ) &
          /sqrt((z1+ym1)*(z2-ym2))
  end do

  goto 7

6 continue

  n1 = j123-jmp
  k=jmn
  i1=nn(n1+1)+k+1
  l=n1-k-k

  if ( l < 0 ) i1=i1+l

  n=jmp
  k=j123-j2-jmn
  i2=nn(n+1)+k+1
  l=n-k-k

  if ( l < 0 ) i2=i2+l

  m=k+1

  kn = j123-j3
  kn = (kn/2)*2-kn

  cg(m) = sqrt( cj*c(i1)*c(i2) )*(-1)**kn

  kn = j2
  kn = min( n1, kn )

  if ( kn == 0 ) goto 14

  ym1 = xm3-yj2
  ym2 = yj2

  cg(m+1)=cg(m)*(xjj-2.0d0*ym1*ym2) &
       /sqrt((yjj1-ym1*ym1-ym1)*(yjj2-ym2*ym2+ym2))

  if ( kn == 1 ) goto 14

  do k=2,kn
     xk=k-1
     ym2=yj2-xk
     ym1=xm3-ym2
     xx=ym1*ym2
     z1=yjj1-ym1*ym1
     z2=yjj2-ym2*ym2
     mk=m+k
     cg(mk)=( cg(mk-1)*(xjj-xx-xx)-cg(mk-2)*sqrt((z1+ym1)*(z2-ym2)) ) &
          /sqrt((z1-ym1)*(z2+ym2))
  end do

  goto 14
7 continue

  if ( m3 /= 0 ) goto 14

  i1 = (j1-j2)/2+1
  i2 = (j1+1)/2
  k=j123/2

  x=(-1.0d0)**(k+k-j123)
  k=j1+2
  do m=i1,i2
     cg(m)=x*cg(k-m)
  end do

14 continue

  if ( ne == 0 ) goto 16

  l=j1+1
  kn=j2+1
  n=j123-j3
  cj=(-1.0d0)**((n/2)*2-n)
  i=j123-jmn+2
  i1=i/2
  kn=min(kn,i1)

  do m=1,kn
     m1=i-m
     if ( m1 < 1 .or. m1 > l ) then
        cg(m)=0.0d0
     else
        xx=cg(m1)*cj
        cg(m1)=cg(m)*cj
        cg(m)=xx
     end if
  end do

16 continue

  kn = nint( xj1+xj1+1 )

  do m=1,kn
     cga(m) = cg(m)
  end do

  if ( nint(xj1) >= nint(abs(xm1)) .and. nint(xj2) >= nint(abs(xm2)) .and. &
       nint(xm1+xm2-xm3) == 0 ) clgr=cga(jm1)

  deallocate( cga )

  return
END FUNCTION clgr
