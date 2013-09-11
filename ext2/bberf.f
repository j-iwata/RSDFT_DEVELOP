      function bberf(arg)
c   the main computation evaluates near-minimax approximations
c   from "rational chebyshev approximations for the error function"
c   by w. j. cody, math. comp., 1969, pp. 631-638.  this
c   transportable program uses rational functions that theoretically
c   approximate  erf(x)  and  erfc(x)  to at least 18 significant
c   decimal digits.  the accuracy achieved depends on the arithmetic
c   system, the compiler, the intrinsic functions, and proper
c   selection of the machine-dependent constants.
c   approximate values for some important machines are:
c
c                          xmin       xinf        xneg     xsmall
c
c  cdc 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
c  cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
c  ieee (ibm/xt,
c    sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
c  ieee (ibm/xt,
c    sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
c  ibm 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
c  univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
c  vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
c  vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16
c
c
c                          xbig       xhuge       xmax
c
c  cdc 7600      (s.p.)  25.922      8.39e+6     1.80x+293
c  cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
c  ieee (ibm/xt,
c    sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
c  ieee (ibm/xt,
c    sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
c  ibm 195       (d.p.)  13.306      1.90d+8     7.23e+75
c  univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
c  vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
c  vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307
c
c  author: w. j. cody
c          mathematics and computer science division
c          argonne national laboratory
c          argonne, il 60439
c  modified by J.Yamauchi
c     1997-09-02
c
c------------------------------------------------------------------
      real*8 bberf,
     1     arg,del,four,half,one,result,sixten,sqrpi,
     2     thresh,x,xbig,xden,xnum,xsmall,y,ysq,zero,
     3     a1,a2,a3,a4,a5,b1,b2,b3,b4,c1,c2,c3,c4,c5,c6,c7,c8,c9,
     4     d1,d2,d3,d4,d5,d6,d7,d8,p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5
c------------------------------------------------------------------
c  mathematical constants
c------------------------------------------------------------------
      parameter (four=4.d0,one=1.0d0,half=0.5d0,zero=0.0d0,
     1     sqrpi=5.6418958354775628695d-1,thresh=0.46875d0,
     2     sixten=16.0d0)
c------------------------------------------------------------------
c  machine-dependent constants
c------------------------------------------------------------------
c for Sun(double)
      parameter( xsmall=1.11d-16, xbig=26.543d0)
c for CRAY
c      parameter( xsmall=7.11e-15, xbig=75.326)
c------------------------------------------------------------------
c  coefficients for approximation to  erf  in first interval
c------------------------------------------------------------------
      parameter (a1=3.16112374387056560d00,a2=1.13864154151050156d02,
     1       a3=3.77485237685302021d02,a4=3.20937758913846947d03,
     2       a5=1.85777706184603153d-1)
      parameter (b1=2.36012909523441209d01,b2=2.44024637934444173d02,
     1       b3=1.28261652607737228d03,b4=2.84423683343917062d03)
c------------------------------------------------------------------
c  coefficients for approximation to  erfc  in second interval
c------------------------------------------------------------------
      parameter ( c1=5.64188496988670089d-1,c2=8.88314979438837594d0,
     1     c3=6.61191906371416295d01,c4=2.98635138197400131d02,
     2     c5=8.81952221241769090d02,c6=1.71204761263407058d03,
     3     c7=2.05107837782607147d03,c8=1.23033935479799725d03,
     4     c9=2.15311535474403846d-8)
      parameter (d1=1.57449261107098347d01,d2=1.17693950891312499d02,
     1     d3=5.37181101862009858d02,d4=1.62138957456669019d03,
     2     d5=3.29079923573345963d03,d6=4.36261909014324716d03,
     3     d7=3.43936767414372164d03,d8=1.23033935480374942d03)
c------------------------------------------------------------------
c  coefficients for approximation to  erfc  in third interval
c------------------------------------------------------------------
      parameter (p1=3.05326634961232344d-1,p2=3.60344899949804439d-1,
     1     p3=1.25781726111229246d-1,p4=1.60837851487422766d-2,
     2     p5=6.58749161529837803d-4,p6=1.63153871373020978d-2)
      parameter (q1=2.56852019228982242d00,q2=1.87295284992346047d00,
     1     q3=5.27905102951428412d-1,q4=6.05183413124413191d-2,
     2     q5=2.33520497626869185d-3)
c------------------------------------------------------------------
      x = arg
      y = abs(x)
      if (y .le. xsmall) then
c------------------------------------------------------------------
c  evaluate  erf  for  |x| <= xsmall
c------------------------------------------------------------------
        result = x * a4 / b4
c------------------------------------------------------------------
c  evaluate  erf  for  xsmall < |x| <= 0.46875
c------------------------------------------------------------------
      else if (y .le. thresh) then
        ysq = y * y
        xnum = a5*ysq
        xden = ysq
        xnum = (((xnum + a1) * ysq + a2) * ysq + a3) * ysq
        xden = (((xden + b1) * ysq + b2) * ysq + b3) * ysq
        result = x * (xnum + a4) / (xden + b4)
c------------------------------------------------------------------
c  evaluate  erfc  for 0.46875 < |x| <= 4.0
c------------------------------------------------------------------
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
c------------------------------------------------------------------
c     evaluate  erfc  for xbig >= |x| > 4.0
c------------------------------------------------------------------
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
c------------------------------------------------------------------
c     evaluate  erfc  for |x| > xbig
c------------------------------------------------------------------
      else
        result = (half - zero) + half
        result = sign(result,x)
      end if
      bberf=result
      return
c---------- last card of calerf ----------
      end
