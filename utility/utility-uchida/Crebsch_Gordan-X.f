!--------1---------2---------3---------4---------5---------6---------7--

      implicit none
      integer :: l1_i,l2_i,m1_i,m2_i,L_i,M_i,j,i
      real(8) :: l1,l2,m1,m2,L,M,l1_max,l2_max
      real(8) :: rm,pi4,bm1,bm2,bM
      real(8) :: c,c1,c2,c3,c4,c5
      complex(8) :: am1,am2,aM,zc1,zc
 
      Integer :: i_0

      INTERFACE
         FUNCTION CG2(j1,j2,m1,m2,jj,mm)
         real(8) :: CG2
         integer,intent(IN) :: j1,j2,m1,m2,jj,mm
         END FUNCTION CG2
      END INTERFACE

      pi4    = atan(1.d0)*16.d0
 

      l1_max = 2
      l2_max = 2

      i=0
      i_0=0

       do l2_i=0,l2_max
       do m2_i=-l2_i,l2_i 

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
     &           *CG2(l1_i,l2_i,0,0,L_i,0)

            zc1=am1*am2*aM*rm
            c2=(bM+bm1*bm2)*CG2(l1_i,l2_i,m1_i,m2_i,L_i,M_i)
            c3=(-1)**M*(1.d0+bM*bm1*bm2)
     &           *CG2(l1_i,l2_i,m1_i,m2_i,L_i,-M_i)
            c4=(-1)**m1_i*(bm2+bM*bm1)
     &           *CG2(l1_i,l2_i,-m1_i,m2_i,L_i,M_i)
            c5=(-1)**m2_i*(bm1+bM*bm2)
     &           *CG2(l1_i,l2_i,m1_i,-m2_i,L_i,M_i)
c            c2=(-1)**M_i*( 1 + (-1)**(m1_i+m2_i-M_i+l1_i+l2_i-L_i)
c     &           *bm1*bm2*bM)*CG2(l1_i,l2_i,m1_i,m2_i,L_i,-M_i)
c            c3=( bM + (-1)**(m1_i+m2_i+M_i+l1_i+l2_i-L_i)
c     &           *bm1*bm2)*CG2(l1_i,l2_i,m1_i,m2_i,L_i,M_i)
c            c4=(-1)**m2_i*(bM*bm2+(-1)**(m1_i-m2_i+M_i+l1_i+l2_i-L_i)
c     &           *bm1)*CG2(l1_i,l2_i,m1_i,-m2_i,L_i,M_i)
c            c5=(-1)**m1_i*( (-1)**(l1_i+l2_i-L_i)*bm1*bM
c     &                     +(-1)**(-m1_i+m2_i+M_i)*bm2)
c     &           *CG2(l1_i,l2_i,m1_i,-m2_i,L_i,-M_i)

            zc=zc1*(c2+c3+c4+c5)
            i_0=i_0+1

            if ( abs(zc)==0.d0 ) cycle
            i=i+1

!            write(*,'(1x,i3,2x,4i3,1x,2i3,1x,2f20.16)')
!     &           i,l1_i,m1_i,L_i,M_i,l2_i,m2_i,real(zc),aimag(zc)

            write(*,'(A,6(i3,A),f20.15,A)')
     &      'C_ijLM(',L_i,',',M_i,',',l1_i,',',m1_i,',',l2_i,','
     &     ,m2_i,'  )= ',real(zc),'d0'




         end do
         end do
         end do
         end do
      end do
      end do

      Write(*,*) 'i,i_0= ',i,i_0
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

!      if ( j1_in<0 .or. j2_in<0 .or. jj_in<0 ) return
!      if ( abs(m1_in)>j1_in .or.
!     &     abs(m2_in)>j2_in .or.
!     &     abs(mm_in)>jj_in     ) return
!      if ( m1_in+m2_in /= mm_in ) return
!     if ( .not.(abs(j1_in-j2_in)<=jj_in.and.jj_in<=j1_in+j2_in)
!     &     ) return

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

       FUNCTION CG2(j1,j2,m1,m2,j3,m3)
      !-----------------------------------
       Implicit NONE

       Integer :: j1,m1,j2,m2,j3,m3
       Integer :: s,z,N
       Real(8) :: CG2
       Real(8) :: v
       Real(8) :: A1
       Real(8) :: A2
       Real(8) :: A3
       Real(8) :: A4
       Real(8) :: A5
       Real(8) :: A6
       Real(8) :: A7
       Real(8) :: A8
       Real(8) :: A9
       Real(8) :: A10
       Real(8) :: A11

       Real(8) :: F1
       Real(8) :: F2
       Real(8) :: F3

       Real(8) :: C1
       Real(8) :: C2
       Real(8) :: C3
       Real(8) :: C4
       Real(8) :: C5
       Real(8) :: C6

 
       Integer :: I1
       Integer :: I2
       Integer :: I3
       Integer :: I4
       Integer :: I5

       Integer :: counter

       counter=0


!--
       CG2=0.d0

      if ( j1<0 .or. j2<0 .or. j3<0 ) return
      if ( abs(m1)>j1 .or.
     &     abs(m2)>j2 .or.
     &     abs(m3)>j3  ) return
      if ( m1+m2 /= m3 ) return
      if ( .not.(abs(j1-j2)<=j3.and.j3<=j1+j2)
     &     ) return
 
!-------------------------------------------------- 
 
       s=j1+j2+j3 

       A1=2*j3+1
       call factorial(A2,s-2*j3)
       call factorial(A3,s-2*j2)
       call factorial(A4,s-2*j1)
       call factorial(A5,s+1)

       call factorial(A6,j1+m1)
       call factorial(A7,j1-m1)
       call factorial(A8,j2+m2)
       call factorial(A9,j2-m2)
       call factorial(A10,j3+m3)
       call factorial(A11,j3-m3)      
       
       F1=sqrt(A1*A2*A3*A4/A5*A6*A7*A8*A9*A10*A11)    

!--
       F2=0.d0
       
       Do z=0,10000
!--


       I1=j1+j2-j3-z     
       I2=j1-m1-z     
       I3=j2+m2-z   
       I4=j3-j2+m1+z     
       I5=j3-j1-m2+z
 
       if (        I1<0 
     &         .or.I2<0 
     &         .or.I3<0 
     &         .or.I4<0 
     &         .or.I5<0 ) cycle

 
      call factorial(C1,I1)
      call factorial(C2,I2)
      call factorial(C3,I3)
      call factorial(C4,I4)
      call factorial(C5,I5)
      call factorial(C6,z)

      F3=((-1)**z)/(C6*C1*C2*C3*C4*C5) 

      F2=F2+F3

      End Do !z


!------------------------------------------------- 
 800   continue


!--------------------------------------------------
 600   continue

       CG2=F1*F2
  
!-------------------------------------------------- 
 900   continue
 
       return 

       End FUNCTION CG2
!--------------------------------
       Subroutine factorial(v,N)
      !------------------------- 
       Implicit NONE
     
       real(8) :: v   
       integer :: N,i

       v=1.d0
 
       if (N<0) then
         write(*,*) 'N<0',N
         stop
       end if

       Do i=1,N
        v=V*dble(i)
       End Do

       Return
       End Subroutine factorial
!--------------------------------
