c     $Id: simc.F,v 1.2 1997/06/25 05:07:08 skimu Exp $
c
c     simc: local potential generator
C           MINPACK version
c
      subroutine simc(rad,vin,rc,zv,parloc,mesh)
      implicit none
      integer*4 mesh
      real*8 rad(mesh),vin(mesh),parloc(4),zv,rc
c
      integer*4 maxnum
      parameter (maxnum=2000)
      common /parusc/ rads,vins,zvs,wgt
      real*8 rads(maxnum),vins(maxnum),zvs,wgt(maxnum)
c
      external uscfit
c
      real*8 lambda
      parameter (lambda = 3.5d0)
      real*8 pi
c
      integer*4 k,num
      real*8 nxtsmp
      real*8 x(3),fvec(maxnum),fjac(maxnum,3)
      integer*4 info
      real*8 x1ini,x2ini,x3ini
      parameter (x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0)
      real*8 tol
      parameter(tol=1.0d-5)
      real*8 rmax,vmax,vrzmin,smpstp
      integer*4 nummin
      parameter(rmax=10.0d0, vmax=100.d0, vrzmin=3.0d-6)
      parameter(smpstp=0.2d0, nummin=6)
      integer*4 ipvt(3)
      integer*4 lwa
      parameter(lwa=5*3+maxnum)
      real*8 wa(lwa)
c
      pi = 4.0d0*atan(1.0d0)
c
      num=0
      nxtsmp=0.d0
      do 11 k=1,mesh
        if(rad(k).gt.nxtsmp) then
          nxtsmp = nxtsmp + smpstp
          if(abs(vin(k)).le.vmax) then
c**   1.0d0*rc ---- 1.5d0*rc in original
            num=num+1
            if (num .gt. maxnum) then
              write(6,*) 'simc: Too many sample points.'
              stop
            end if
            rads(num)=rad(k)
            vins(num)=vin(k)
            wgt(num)=1.d0-dexp(-(1.2d0*rad(k)/rc)**lambda)
          end if
        if ((abs(vin(k)*rad(k)+zv).lt.vrzmin .or. rad(k).gt.rmax)
     &         .and. num.gt.nummin) then
          goto 100
        end if
        end if
   11 continue
  100 continue
      zvs = zv
      x(1) = x1ini
      x(2) = x2ini
      x(3) = x3ini
      call lmder1(uscfit,num,3,x,fvec,fjac,maxnum,tol,info,ipvt,wa,
     &                  lwa)
      write(6,*) 'lmder1:',info
      if (info.eq.0 .or. info.eq.4 .or. info.eq.5 .or. info.eq.6 .or.
     &     info.eq.7) then
        write(6,*) 'simc: Not converged.'
        write(6,*) 'x(1) = ',x(1)
        write(6,*) 'x(2) = ',x(2)
        write(6,*) 'x(3) = ',x(3)
        write(6,*) 'k, rads(k),vins(k),wgt(k)'
        do 12 k=1,num
          write(6,*)k,rads(k),vins(k),wgt(k)
   12   continue
        stop
      end if
      if (x(2).lt.0.0d0 .or. x(3).lt.0.0d0) then
        write(6,*)'simc: illegally converged.'
        stop
      end if
      parloc(1) = x(1)
      parloc(2) = x(2)
      parloc(3) = 1.0d0 - x(1)
      parloc(4) = x(3)
      return
      end
c
c     fitting function for simc
c
      subroutine uscfit(m,n,x,fvec,fjac,ldfjac,iflag)
      implicit none
      integer*4 m,n,ldfjac,iflag
      real*8 x(n),fvec(m),fjac(ldfjac,n)
      real*8 bberf
      integer*4 maxnum
      parameter (maxnum=2000)
      common /parusc/ rad,vin,zv,wgt
      real*8 rad(maxnum),vin(maxnum),zv,wgt(maxnum)

      real*8 pi

      integer*4 i

      pi = 4.0d0*atan(1.0d0)

      if (x(2).lt.0.0d0) x(2)=0.0d0
      if (x(3).lt.0.0d0) x(3)=0.0d0
      if (iflag .eq. 1) then
        do 13 i=1,m
          fvec(i) = (- zv/rad(i)*
     &         (x(1)*bberf(sqrt(x(2))*rad(i))
     &          + (1.0d0-x(1))*bberf(sqrt(x(3))*rad(i)))
     &         - vin(i))*sqrt(wgt(i))
   13   continue
      else if (iflag .eq. 2) then
        do 14 i=1,m
          fjac(i,1) = -zv/rad(i)
     &         *(bberf(sqrt(x(2)*rad(i)))-bberf(sqrt(x(3)*rad(i))))
     &         *sqrt(wgt(i))
          fjac(i,2) = -zv/sqrt(pi*x(2))*x(1)
     &         *exp(-x(2)*rad(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zv/sqrt(pi*x(3))*(1.0d0-x(1))
     &         *exp(-x(3)*rad(i)**2)*sqrt(wgt(i))
   14   continue
      else
        write(6,*) 'Error in vlfit: iflag must be 1 or 2.'
        stop
      end if
      return
      end
