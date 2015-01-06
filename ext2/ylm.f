!--------1---------2---------3---------4---------5---------6---------7--
! Spherical Harmonic Function
!
      FUNCTION Ylm(x,y,z,L,M_in)
      implicit none
      real(8) :: Ylm
      real(8),intent(IN) :: x,y,z
      integer,intent(IN) :: L,M_in
      real(8) :: r,r2,Clm,ct,phi,c,pi4
      integer :: k1,k2,k,M
      real(8),parameter :: ep=1.d-10,zero=0.d0
      real(8),parameter :: half=0.5d0,one=1.d0,two=2.d0,four=4.d0
      real(8),parameter :: pi=3.141592653589793238d0

      INTERFACE
         FUNCTION plgndr(L,M,x)
         real(8) :: plgndr,x
         integer :: L,M
         intent(IN) :: L,M
         END FUNCTION plgndr
      END INTERFACE

      Ylm = zero
      M   = abs(M_in)
      pi4 = four*pi

      if( L<0 .or. M>L )then
         write(*,*) "L,M_in=",L,M_in
         stop "Bad arguments (stop at Ylm)"
      end if

      if( L==0 )then
         Ylm=one/sqrt(pi4)
         return
      end if

      if( x==0.d0 .and. y==0.d0 .and. z==0.d0 )then
         return
      end if

      c=one
      do k=L-M+1,L+M
         c=c*dble(k)
      end do
      Clm=sqrt((two*L+one)/(c*pi4))

      r  = sqrt(x*x+y*y+z*z)
      ct = z/r
      if( abs(x)<ep )then
         phi=half*pi*sign(one,y)
      else
         phi=atan(y/x)
         if ( x<0 ) phi=phi+pi
      end if

      if ( M_in>0 ) then

         Ylm = sqrt(two)*Clm*plgndr(L,M,ct)*cos(phi*M)

      else if ( M_in==0 ) then

         Ylm = Clm*plgndr(L,M,ct)

      else if ( mod(M,2)==0 ) then

         Ylm =-sqrt(two)*Clm*plgndr(L,M,ct)*sin(phi*M)

      else

         Ylm = sqrt(two)*Clm*plgndr(L,M,ct)*sin(phi*M)

      end if

      return
      END FUNCTION Ylm

!--------1---------2---------3---------4---------5---------6---------7--
! Associated Legendre Polynomial Plm(x)
!     (from Numerical Recipes)
!
      FUNCTION plgndr(L,M,x)
      implicit none
      real(8) :: plgndr,x ,pmm,somx2,fact,pmmp1,pll,tmp
      real(8),parameter :: ep=1.d-13,one=1.d0,two=2.d0
      integer :: L,M      ,i,ll
      intent(IN) :: L,M

      tmp=abs(abs(x)-one) ; if ( tmp<ep ) x=sign(one,x)

      if( m<0 .or. m>l .or. abs(x)>one )then
         stop "Bad arguments (stop at plgndr)"
      end if
      pmm=one
      if( m>0 )then
         somx2=sqrt((one-x)*(one+x))
         fact=one
         do i=1,M
            pmm=-pmm*fact*somx2
            fact=fact+two
         end do
      end if
      if( L==M )then
         plgndr=pmm
      else
         pmmp1=x*(two*m+one)*pmm
         if( L==M+1 )then
            plgndr=pmmp1
         else
            do ll=M+2,L
               pll=(x*(two*ll-one)*pmmp1-(ll+M-one)*pmm)/(ll-M)
               pmm=pmmp1
               pmmp1=pll
            end do
            plgndr=pll
         end if
      end if
      return
      END FUNCTION plgndr
