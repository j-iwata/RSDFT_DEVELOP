MODULE modified_bessel_module

  implicit none

  PRIVATE
  PUBLIC :: bessi,bessk

CONTAINS

  FUNCTION bessi(n,x)
    integer :: n
    real(8) :: bessi,x
    integer,parameter :: IACC=40
    real(8),parameter :: BIGNO=1.0d10, BIGNI=1.0d-10
    integer :: j,m
    real(8) :: bi,bim,bip,tox
    if ( n == 0 ) then
       bessi = bessi0(x)
       return
    else if ( n == 1 ) then
       bessi = bessi1(x)
       return
    else
       if ( x == 0.d0 ) then
          bessi=0.d0
       else
          tox=2.0d0/abs(x)
          bip=0.0d0
          bi =1.0d0
          bessi=0.0d0
          m=2*((n+int(sqrt(dble(IACC*n)))))
          do j=m,1,-1
             bim=bip+dble(j)*tox*bi
             bip=bi
             bi=bim
             if ( abs(bi) > BIGNO ) then
                bessi = bessi*BIGNI
                bi=bi*BIGNI
                bip=bip*BIGNI
             endif
             if ( j == n ) bessi=bip
          enddo
          bessi=bessi*bessi0(x)/bi
          if ( x < 0.0d0 .and. mod(n,2) == 1 ) bessi=-bessi
       endif
    end if
    return
  END FUNCTION bessi

  FUNCTION bessk(n,x)
    integer :: n
    real(8) :: bessk,x
    integer :: j
    real(8) :: bk,bkm,bkp,tox
    if ( n == 0 ) then
       bessk = bessk0(x)
    else if ( n == 1 ) then
       bessk = bessk1(x)
    else
       tox=2.d0
       bkm=bessk0(x)
       bk=bessk1(x)
       do j=1,n-1
          bkp=bkm+j*tox*bk
          bkm=bk
          bk=bkp
       end do
       bessk=bk
       return
    end if
  END FUNCTION bessk

  FUNCTION bessi0(x)
    real(8) :: bessi0,x
    real(8) :: ax,y
    real(8),parameter :: p1=1.0d0, p2=3.5156229d0, p3=3.0899424d0
    real(8),parameter :: p4=1.2067492d0, p5=0.2659732d0
    real(8),parameter :: p6=0.36768d-1, p7=0.45813d-2
    real(8),parameter :: q1=0.39894228d0, q2=0.1328592d-1
    real(8),parameter :: q3=0.225319d-2 , q4=-0.157565d-2
    real(8),parameter :: q5=0.916281d-2 , q6=-0.2057706d-1
    real(8),parameter :: q7=0.2635537d-1, q8=-0.1647633d-1, q9=0.392377d-2
    if ( abs(x) < 3.75d0 ) then
       y=(x/3.75d0)**2
       bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    else
       ax=abs(x)
       y=3.75d0/ax
       bessi0 = ( exp(ax)/sqrt(ax) ) &
            *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
    end if
    return
  END FUNCTION bessi0

  FUNCTION bessk0(x)
    real(8) :: bessk0,x
    real(8),parameter :: p1=-0.57721566d0, p2=0.42278420d0, p3=0.23069756d0
    real(8),parameter :: p4= 0.3488590d-1, p5=0.262698d-2 , p6=0.10750d-3
    real(8),parameter :: p7= 0.74d-5
    real(8),parameter :: q1= 1.25331414d0, q2=-0.7832358d-1, q3=0.2189568d-1
    real(8),parameter :: q4=-0.1062446d-1, q5=0.587872d-2
    real(8),parameter :: q6=-0.251540d-2 , q7=0.53208d-3
    real(8) :: y
    if ( x <= 2.d0 ) then
       y=x*x/4.d0
       bessk0 = ( -log(x/2.d0)*bessi0(x) ) &
            + (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       y=(2.d0/x)
       bessk0 = ( exp(-x)/sqrt(x) ) &
            *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
    end if
    return
  END FUNCTION bessk0

  FUNCTION bessi1(x)
    real(8) :: bessi1,x
    real(8) :: ax,y
    real(8),parameter :: p1=0.5d0, p2=0.87890594d0, p3=0.51498869d0
    real(8),parameter :: p4=0.15084934d0, p5=0.2658733d-1
    real(8),parameter :: p6=0.301532d-2, p7=0.32411d-3
    real(8),parameter :: q1= 0.39894228d0, q2=-0.3988024d-1
    real(8),parameter :: q3=-0.362018d-2 , q4= 0.163801d-2
    real(8),parameter :: q5=-0.1031555d-1, q6= 0.2282967d-1
    real(8),parameter :: q7=-0.2895312d-1, q8= 0.1787654d-1, q9=-0.420059d-2
    if ( abs(x) < 3.75d0 ) then
       y=(x/3.75d0)**2
       bessi1 = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       ax=abs(x)
       y=3.75d0/ax
       bessi1 = ( exp(ax)/sqrt(ax) ) &
            *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
    end if
    return
  END FUNCTION bessi1

  FUNCTION bessk1(x)
    real(8) :: bessk1,x
    real(8),parameter :: p1=1.0d0, p2=0.15443144d0, p3=-0.67278579d0
    real(8),parameter :: p4=-0.18156897d0, p5=-0.1919402d-1
    real(8),parameter :: p6=-0.110404d-2, p7=-0.4686d-4
    real(8),parameter :: q1=1.25331414d0, q2=0.23498619d0
    real(8),parameter :: q3=-0.3655620d-1, q4=0.1504268d-1
    real(8),parameter :: q5=-0.780353d-2, q6=0.325614d-2, q7=-0.68245d-3
    real(8) :: y
    if ( x <= 2.d0 ) then
       y=x*x/4.d0
       bessk1 = ( log(x/2.d0)*bessi1(x) ) &
            + (1.d0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       y=2.d0/x
       bessk1 = ( exp(-x)/sqrt(x) ) &
            *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
    end if
    return
  END FUNCTION bessk1

END MODULE modified_bessel_module
