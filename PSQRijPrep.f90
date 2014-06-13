MODULE PSQRijPrep
  implicit none
  
  
  complex(8),allocatable :: qaL(:,:) !qaL(k3max,Lrefmax)

  complex(8),parameter :: z0=(0.d0,0.d0),z1=(1.d0,0.d0),zi=(0.d0,1.d0)

CONTAINS

  SUBROUTINE prep_QRij_p12
    implicit none

  END SUBROUTINE prep_QRij_p12

!---------------------------------------------------------------------------------------

  SUBROUTINE get_qaL ( r,x,y,z )
    implicit none

    real(8) :: r,x,y,z
    real(8),parameter :: sq3=sqrt(3.d0)
    real(8),parameter :: sq5=sqrt(5.d0)

!----------------------------
    qaL(:,:)   = z0

    qaL(1,1)   = z1
    qaL(3,1)   = z1
    qaL(6,1)   = z1
    qaL(10,1)  = z1
    qaL(15,1)  = z1
    qaL(21,1)  = z1
    qaL(28,1)  = z1
    qaL(36,1)  = z1
    qaL(45,1)  = z1

!--      
    if ( r <= 1.0D-10 ) return
!--      
    x = x/r 
    y = y/r 
    z = z/r 

    qaL(1,1)=1.d0
    qaL(2,1)=sq3*zi*y
    qaL(3,1)=1.d0
    qaL(3,2)=(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/2.d0
    qaL(4,1)=-sq3*zi*z
    qaL(5,1)=0.d0
    qaL(5,2)=3.d0*y*z
    qaL(6,1)=1.d0
    qaL(6,2)=-(3.d0*z**2-1.d0)
    qaL(7,1)=sq3*zi*x
    qaL(8,1)=0.d0
    qaL(8,2)=-3.d0*x*y
    qaL(9,1)=0.d0
    qaL(9,2)=3.d0*x*z
    qaL(10,1)=1.d0
    qaL(10,2)=-(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0)/2.d0
    qaL(11,1)=-sqrt(15.d0)*x*y
    qaL(12,1)=3.d0*sq5*zi*x/5.d0
    qaL(12,2)=(   3.d0*sq5*zi*x * (5.d0*x**2-15.d0*y**2+5.d0*z**2-1.d0) )/20.d0
    qaL(13,1)=0.d0
    qaL(13,2)=3.d0*sq5*zi*x*y*z
    qaL(14,1)=3.d0*sq5*zi*y/5.d0
    qaL(14,2)=-(3.d0*sq5*zi*y*(15.d0*x**2 -5.d0*y**2-5.d0*z**2+1.d0))/(20.d0)
    qaL(15,1)=1.d0
    qaL(15,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(15,3)=-(3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4-35.d0*z**4+30.d0*z**2-3.d0))/56.d0
    qaL(16,1)=3.d0*sqrt(5.d0/3.0d0)*y*z
    qaL(17,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(17,2)=-(3.d0*sq5*zi*z *(5.d0*x**2-5.d0*y**2+5.d0*z**2-3.d0))/10.d0
    qaL(18,1)=(3.d0*sq5*zi*y)/5.d0
    qaL(18,2)=-(3.d0*zi*y*(5.d0*z**2-1.d0))/sq5
    qaL(19,1)=0.d0
    qaL(19,2)=3.d0*sq5*zi*x*y*z
    qaL(20,1)=0.d0
    qaL(20,2)=(15.d0*x*z)/7.d0
    qaL(20,3)=(15.d0*x*z *(7.d0*x**2-21.d0*y**2+7.d0*z**2-3.d0))/28.d0
    qaL(21,1)=1.d0
    qaL(21,2)=(5.d0*(3.d0*x**2-3.d0*y**2-3.d0*z**2+1.d0))/14.d0
    qaL(21,3)=-(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2+35.d0*z**4-30.d0*z**2+3.d0))/14.d0
    qaL(22,1)=-(sq5*(3.d0*z**2-1.d0))/2.d0
    qaL(23,1)=-(sqrt(15.d0)*zi*y)/5.d0
    qaL(23,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*y*(5.d0*z**2-1.d0)
    qaL(24,1)=-(2.d0*sqrt(15.d0)*zi*z)/5.d0
    qaL(24,2)=9.d0*zi*z*(5.d0*z**2-3.d0)/(2.d0*sqrt(15.d0))
    qaL(25,1)=-(sqrt(15.d0)*zi*x)/5.d0
    qaL(25,2)=-1.5d0*sqrt(3.d0/5.d0)*zi*x*(5.d0*z**2-1.d0)
    qaL(26,1)=0.d0
    qaL(26,2)=10.d0*sq3*x*y/7.d0
    qaL(26,3)=15.d0*sq3*x*y*(7.d0*z**2-1.d0)/14.d0
    qaL(27,1)=0.d0
    qaL(27,2)=(15.d0*y*z)/(7.d0*sq3)
    qaL(27,3)=-15.d0*sq3*y*z*(7.d0*z**2-3.d0)/14.d0
    qaL(28,1)=1.d0
    qaL(28,2)=-(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(28,3)=(9.d0*(35.d0*z**4-30.d0*z**2+3.d0))/28.d0
    qaL(29,1)=(3.d0*sqrt(5.d0/3.d0)*x*z)
    qaL(30,1)=0.d0
    qaL(30,2)=3.d0*sq5*zi*x*y*z
    qaL(31,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(31,2)=-(3.d0*zi*x*(5.d0*z**2-1.d0))/sq5
    qaL(32,1)=-(3.d0*sq5*zi*z)/5.d0
    qaL(32,2)=0.3d0*sq5*zi*z*(5.d0*x**2-5.d0*y**2-5.d0*z**2+3.d0)
    qaL(33,1)=0.d0
    qaL(33,2)=(15.d0*y*z)/7.d0
    qaL(33,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(34,1)=0.d0
    qaL(34,2)=-(15.d0*x*y)/7.d0
    qaL(34,3)=(15.d0*x*y*(7.d0*z**2-1.d0))/7.d0
    qaL(35,1)=0.d0
    qaL(35,2)=(15.d0*x*z)/(7.d0*sq3)
    qaL(35,3)=-15.d0*sq3*x*z*(7.d0*z**2-3.d0)/14.d0
    qaL(36,1)=1.d0
    qaL(36,2)=-5.d0*(3.d0*x**2-3.d0*y**2+3.d0*z**2-1.d0)/14.d0
    qaL(36,3)=(3.d0*(35.d0*x**2*z**2-5.d0*x**2-35.d0*y**2*z**2 +5.d0*y**2-35.d0*z**4+30.d0*z**2-3.d0))/14.d0
    qaL(37,1)=-(3.d0*sq5*(x**2-y**2))/(2.d0*sq3)
    qaL(38,1)=-(3.d0*sqrt(10.d0)*zi*y)/(5.d0*sqrt(2.d0))
    qaL(38,2)=-0.15d0*sq5*zi*y *(15.d0*x**2-5.d0*y**2+5.d0*z**2-1.d0)
    qaL(39,1)=0.d0
    qaL(39,2)=1.5d0*sq5*zi*z*(x**2-y**2)
    qaL(40,1)=(3.d0*sq5*zi*x)/5.d0
    qaL(40,2)=-0.15d0*sq5*zi*x  *(5.d0*x**2-15.d0*y**2-5.d0*z**2+1.d0)
    qaL(41,1)=0.d0
    qaL(41,2)=0.d0
    qaL(41,3)=7.5d0*x*y*(x**2-y**2)
    qaL(42,1)=0.d0
    qaL(42,2)=-(15.d0*y*z)/7.d0
    qaL(42,3)=-15.d0*y*z *(21.d0*x**2-7.d0*y**2+7.d0*z**2-3.d0)/28.d0
    qaL(43,1)=0.d0
    qaL(43,2)=(15.d0*(x**2-y**2))/(7.d0*sq3)
    qaL(43,3)=15.d0*sq3 *(7.d0*x**2*z**2-x**2-7.d0*y**2*z**2+y**2)/28.d0
    qaL(44,1)=0.d0
    qaL(44,2)=(15.d0*x*z)/7.d0
    qaL(44,3)=-15.d0*x*z *(7.d0*x**2-21.d0*y**2-7.d0*z**2+3.d0)/28.d0
    qaL(45,1)=1.d0
    qaL(45,2)=(5.d0*(3.d0*z**2-1.d0))/7.d0
    qaL(45,3)=3.d0*(35.d0*x**4-210.d0*x**2*y**2 +35.d0*y**4+35.d0*z**4-30.d0*z**2+3.d0)/56.d0

!--    THIS will be done when qaL is used.
!      qaL(:,:)=qaL(:,:)/( (4.d0*Pi)*(-zi)**L )

    return

  END SUBROUTINE get_qaL

END MODULE PSQRijPrep
