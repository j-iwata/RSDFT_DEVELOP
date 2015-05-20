MODULE bberf_module

  PRIVATE
  PUBLIC :: bberf

CONTAINS

#ifndef _COPYRIGHTFREE_
!----------------------------------------------------
! The copyright of the following function is unknown.
!----------------------------------------------------

  function bberf(arg)
!   the main computation evaluates near-minimax approximations
!   from "rational chebyshev approximations for the error function"
!   by w. j. cody, math. comp., 1969, pp. 631-638.  this
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  the accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!   approximate values for some important machines are:
!
!                          xmin       xinf        xneg     xsmall
!
!  cdc 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
!  cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
!  ieee (ibm/xt,
!    sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
!  ieee (ibm/xt,
!    sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
!  ibm 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
!  univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
!  vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
!  vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16
!
!
!                          xbig       xhuge       xmax
!
!  cdc 7600      (s.p.)  25.922      8.39e+6     1.80x+293
!  cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
!  ieee (ibm/xt,
!    sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
!  ieee (ibm/xt,
!    sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
!  ibm 195       (d.p.)  13.306      1.90d+8     7.23e+75
!  univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
!  vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
!  vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307
!
!  author: w. j. cody
!          mathematics and computer science division
!          argonne national laboratory
!          argonne, il 60439
!  modified by J.Yamauchi
!     1997-09-02
!
!------------------------------------------------------------------
    real*8 bberf, &
         arg,del,four,half,one,result,sixten,sqrpi, &
         thresh,x,xbig,xden,xnum,xsmall,y,ysq,zero, &
         a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,c6,c7,c8,c9, &
         d1,d2,d3,d4,d5,d6,d7,d8,p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5
!------------------------------------------------------------------
!  mathematical constants
!------------------------------------------------------------------
    parameter (four=4.d0,one=1.0d0,half=0.5d0,zero=0.0d0, &
         sqrpi=5.6418958354775628695d-1,thresh=0.46875d0, &
         sixten=16.0d0)
!------------------------------------------------------------------
!  machine-dependent constants
!------------------------------------------------------------------
! for Sun(double)
    parameter( xsmall=1.11d-16, xbig=26.543d0)
! for CRAY
!      parameter( xsmall=7.11e-15, xbig=75.326)
!------------------------------------------------------------------
!  coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
    parameter (a1=3.16112374387056560d00,a2=1.13864154151050156d02, &
         a3=3.77485237685302021d02,a4=3.20937758913846947d03, &
         a5=1.85777706184603153d-1)
    parameter (b1=2.36012909523441209d01,b2=2.44024637934444173d02, &
         b3=1.28261652607737228d03,b4=2.84423683343917062d03)
!------------------------------------------------------------------
!  coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
    parameter ( c1=5.64188496988670089d-1,c2=8.88314979438837594d0, &
         c3=6.61191906371416295d01,c4=2.98635138197400131d02, &
         c5=8.81952221241769090d02,c6=1.71204761263407058d03, &
         c7=2.05107837782607147d03,c8=1.23033935479799725d03, &
         c9=2.15311535474403846d-8)
    parameter (d1=1.57449261107098347d01,d2=1.17693950891312499d02, &
         d3=5.37181101862009858d02,d4=1.62138957456669019d03, &
         d5=3.29079923573345963d03,d6=4.36261909014324716d03, &
         d7=3.43936767414372164d03,d8=1.23033935480374942d03)
!------------------------------------------------------------------
!  coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
    parameter (p1=3.05326634961232344d-1,p2=3.60344899949804439d-1, &
         p3=1.25781726111229246d-1,p4=1.60837851487422766d-2, &
         p5=6.58749161529837803d-4,p6=1.63153871373020978d-2)
    parameter (q1=2.56852019228982242d00,q2=1.87295284992346047d00, &
         q3=5.27905102951428412d-1,q4=6.05183413124413191d-2, &
         q5=2.33520497626869185d-3)
!------------------------------------------------------------------
    x = arg
    y = abs(x)
    if (y .le. xsmall) then
!------------------------------------------------------------------
!  evaluate  erf  for  |x| <= xsmall
!------------------------------------------------------------------
       result = x * a4 / b4
!------------------------------------------------------------------
!  evaluate  erf  for  xsmall < |x| <= 0.46875
!------------------------------------------------------------------
    else if (y .le. thresh) then
       ysq = y * y
       xnum = a5*ysq
       xden = ysq
       xnum = (((xnum + a1) * ysq + a2) * ysq + a3) * ysq
       xden = (((xden + b1) * ysq + b2) * ysq + b3) * ysq
       result = x * (xnum + a4) / (xden + b4)
!------------------------------------------------------------------
!  evaluate  erfc  for 0.46875 < |x| <= 4.0
!------------------------------------------------------------------
    else if (y .le. four) then
       xnum = c9*y
       xden = y
       xnum = (xnum + c1) * y
       xden = (xden + d1) * y
       xnum = (xnum + c2) * y
       xden = (xden + d2) * y
       xnum = (xnum + c3) * y
       xden = (xden + d3) * y
       xnum = (xnum + c4) * y
       xden = (xden + d4) * y
       xnum = (xnum + c5) * y
       xden = (xden + d5) * y
       xnum = (xnum + c6) * y
       xden = (xden + d6) * y
       xnum = (xnum + c7) * y
       xden = (xden + d7) * y
       result = (xnum + c8) / (xden + d8)
       ysq = aint(y*sixten)/sixten
       del = (y-ysq)*(y+ysq)
       result = exp(-ysq*ysq) * exp(-del) * result
       result = (half - result) + half
       result = sign(result,x)
!------------------------------------------------------------------
!     evaluate  erfc  for xbig >= |x| > 4.0
!------------------------------------------------------------------
    else if (y.le.xbig) then
       ysq = one / (y * y)
       xnum = p6*ysq
       xden = ysq
       xnum = (xnum + p1) * ysq
       xden = (xden + q1) * ysq
       xnum = (xnum + p2) * ysq
       xden = (xden + q2) * ysq
       xnum = ((xnum + p3) * ysq + p4) * ysq
       xden = ((xden + q3) * ysq + q4) * ysq
       result = ysq *(xnum + p5) / (xden + q5)
       result = (sqrpi -  result) / y
       ysq = aint(y*sixten)/sixten
       del = (y-ysq)*(y+ysq)
       result = exp(-ysq*ysq) * exp(-del) * result
       result = (half - result) + half
       result = sign(result,x)
!------------------------------------------------------------------
!     evaluate  erfc  for |x| > xbig
!------------------------------------------------------------------
    else
       result = (half - zero) + half
       result = sign(result,x)
    end if
    bberf=result
    return
!---------- last card of calerf ----------
  end function bberf

#elif _COPYRIGHTFREE_
!------------------------------------------------------
! The following functions are not copies of something.
! They are written by J. Iwata.
!------------------------------------------------------

  FUNCTION bberf( x )
    implicit none
    real(8) :: bberf_ji
    real(8),intent(IN) :: x
    bberf_ji=1.0d0-ccerf(x)
  END FUNCTION bberf

  FUNCTION ccerf( x )

    implicit none
    real(8) :: ccerf
    real(8),intent(IN) :: x
    real(8),parameter :: pi = 3.141592653589793d0
    real(8) :: y,alpha,tny,eps,delta
    real(8) :: ai,bi,f0,f1,C0,C1,D0,D1,a0
    integer :: nmax,i

    nmax  = 100000
    eps   = epsilon(1.0d0)
    tny   = tiny(1.0d0)
    alpha = 0.5d0

    y = x*x

    f0 = tny
    C0 = f0
    D0 = 0.0d0
    a0 = 1.0d0

    do i=1,nmax
       bi = y + 2.0d0*i - 1.0d0 - alpha
       ai = -(i-1)*(i-1-alpha) + a0
       ai = -(i-1)*(i-1.5d0) + a0
       D1 = bi + ai*D0
       if ( D1 == 0.0d0 ) D1=tny
       C1 = bi + ai/C0
       if ( C1 == 0.0d0 ) C1=tny
       D1 = 1.0d0/D1
       delta = C1*D1
       f1 = f0*delta
       if ( abs(delta-1.0d0) < eps ) exit
       f0 = f1
       C0 = C1
       D0 = D1
       a0 = 0.0d0
    end do

    ccerf = f1*exp(-y)*sqrt(y/pi)

  END FUNCTION ccerf
#endif

END MODULE bberf_module
