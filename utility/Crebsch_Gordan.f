!--------1---------2---------3---------4---------5---------6---------7--

      implicit none
      integer :: l1_i,l2_i,m1_i,m2_i,L_i,M_i,j,i
      real(8) :: l1,l2,m1,m2,L,M,l1_max
      real(8) :: rm,pi4,bm1,bm2,bM
      real(8) :: c,c1,c2,c3,c4,c5
      complex(8) :: am1,am2,aM,zc1,zc

      INTERFACE
         FUNCTION CG(j1,j2,m1,m2,jj,mm)
         real(8) :: CG
         integer,intent(IN) :: j1,j2,m1,m2,jj,mm
         END FUNCTION CG
      END INTERFACE

      pi4    = atan(1.d0)*16.d0
      l2_i   = 1
      l1_max = 3

      i=0
      do j=1,3
         select case(j)
         case(1) ; m2_i=1
         case(2) ; m2_i=-1
         case(3) ; m2_i=0
         end select
         do l1_i=0,l1_max
         do L_i=abs(l1_i-l2_i),l1_i+l2_i
         do M_i=-L_i,L_i
         do m1_i=-l1_i,l1_i

            l1 = l1_i
            l2 = l2_i
            L  = L_i
            m1 = m1_i
            m2 = m2_i
            M  = M_i

            select case(m1_i)
            case(1:)
               am1 = 1.d0/sqrt(2.d0)
               bm1 = 1.d0
            case(0)
               am1 = 1.d0
               bm1 = 0.d0
            case(:-1)
               am1 = dcmplx(0.d0,-1.d0/sqrt(2.d0))
               bm1 =-1.d0
            end select

            select case(m2_i)
            case(1:)
               am2 = 1.d0/sqrt(2.d0)
               bm2 = 1.d0
            case(0)
               am2 = 1.d0
               bm2 = 0.d0
            case(:-1)
               am2 = dcmplx(0.d0,-1.d0/sqrt(2.d0))
               bm2 =-1.d0
            end select

            select case(M_i)
            case(1:)
               aM = 1.d0/sqrt(2.d0)
               bM = 1.d0
            case(0)
               aM = 1.d0
               bM = 0.d0
            case(:-1)
               aM = dcmplx(0.d0,-1.d0/sqrt(2.d0))
               bM =-1.d0
            end select

            rm=sqrt( (2.d0*l1+1.d0)*(2.d0*l2+1.d0)/pi4/(2.d0*L+1.d0) )
     &           *CG(l1_i,l2_i,0,0,L_i,0)

            zc1=am1*am2*aM*rm
            c2=(bM+bm1*bm2)*CG(l1_i,l2_i,m1_i,m2_i,L_i,M_i)
            c3=(-1)**M*(1.d0+bM*bm1*bm2)
     &           *CG(l1_i,l2_i,m1_i,m2_i,L_i,-M_i)
            c4=(-1)**m1_i*(bm2+bM*bm1)
     &           *CG(l1_i,l2_i,-m1_i,m2_i,L_i,M_i)
            c5=(-1)**m2_i*(bm1+bM*bm2)
     &           *CG(l1_i,l2_i,m1_i,-m2_i,L_i,M_i)
c            c2=(-1)**M_i*( 1 + (-1)**(m1_i+m2_i-M_i+l1_i+l2_i-L_i)
c     &           *bm1*bm2*bM)*CG(l1_i,l2_i,m1_i,m2_i,L_i,-M_i)
c            c3=( bM + (-1)**(m1_i+m2_i+M_i+l1_i+l2_i-L_i)
c     &           *bm1*bm2)*CG(l1_i,l2_i,m1_i,m2_i,L_i,M_i)
c            c4=(-1)**m2_i*(bM*bm2+(-1)**(m1_i-m2_i+M_i+l1_i+l2_i-L_i)
c     &           *bm1)*CG(l1_i,l2_i,m1_i,-m2_i,L_i,M_i)
c            c5=(-1)**m1_i*( (-1)**(l1_i+l2_i-L_i)*bm1*bM
c     &                     +(-1)**(-m1_i+m2_i+M_i)*bm2)
c     &           *CG(l1_i,l2_i,m1_i,-m2_i,L_i,-M_i)

            zc=zc1*(c2+c3+c4+c5)
            if ( abs(zc)==0.d0 ) cycle
            i=i+1
            write(*,'(1x,i3,2x,4i3,1x,2i3,1x,2f20.16)')
     &           i,l1_i,m1_i,L_i,M_i,l2_i,m2_i,real(zc),aimag(zc)

         end do
         end do
         end do
         end do
      end do

      stop
      end


!--------1---------2---------3---------4---------5---------6---------7--

      FUNCTION CG(j1_in,j2_in,m1_in,m2_in,jj_in,mm_in)
      implicit none
      real(8) :: CG
      integer,intent(IN) :: j1_in,j2_in,m1_in,m2_in,jj_in,mm_in
      real(8),parameter :: one=1.d0,two=2.d0
      real(8) :: j1,j2,m1,m2,jj,mm

      CG=0.d0

      if ( j1_in<0 .or. j2_in<0 .or. jj_in<0 ) return
      if ( abs(m1_in)>j1_in .or.
     &     abs(m2_in)>j2_in .or.
     &     abs(mm_in)>jj_in     ) return
      if ( m1_in+m2_in /= mm_in ) return
      if ( .not.(abs(j1_in-j2_in)<=jj_in.and.jj_in<=j1_in+j2_in)
     &     ) return

!
! The following is a temporary restriction
!
      if ( j2_in/=1 ) then
         write(*,*) "WARNING! j2/=1 is not available yet"
         return
      end if

!
! --- CG for j2=1 ---
!
      j1 = j1_in
      j2 = j2_in
      jj = jj_in
      m1 = m1_in
      m2 = m2_in
      mm = mm_in

      if ( jj_in==j1_in+1 ) then

         select case(m2_in)
         case( 1 )
            CG=sqrt( (j1+mm)*(j1+mm+one)/(two*j1+one)/(two*j1+two) )
         case( 0 )
            CG=sqrt( (j1-mm+one)*(j1+mm+one)/(two*j1+one)/(j1+one) )
         case(-1 )
            CG=sqrt( (j1-mm)*(j1-mm+one)/(two*j1+one)/(two*j1+two) )
         end select

      else if ( jj_in==j1_in ) then

         select case(m2_in)
         case( 1 )
            CG=-sqrt( (j1+mm)*(j1-mm+one)/(two*j1)/(j1+one) )
         case( 0 )
            CG=mm/sqrt( j1*(j1+one) )
         case(-1 )
            CG=sqrt( (j1-mm)*(j1+mm+one)/(two*j1)/(j1+one) )
         end select

      else if ( jj_in==j1_in-1 ) then

         select case(m2_in)
         case( 1 )
            CG=sqrt( (j1-mm)*(j1-mm+one)/(two*j1)/(two*j1+one) )
         case( 0 )
            CG=-sqrt( (j1-mm)*(j1+mm)/j1/(two*j1+one) )
         case(-1 )
            CG=sqrt( (j1+mm+one)*(j1+mm)/(two*j1)/(two*j1+one) )
         end select

      end if

      return
      END FUNCTION CG
