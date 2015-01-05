MODULE simc_module

  implicit none

  PRIVATE
  PUBLIC :: simc

  INTERFACE
     FUNCTION bberf(x)
       real(8) :: bberf,x
     END FUNCTION bberf
  END INTERFACE

CONTAINS


  SUBROUTINE simc(rad,vin,rc,zv,parloc,mesh)
!     $Id: simc.F,v 1.2 1997/06/25 05:07:08 skimu Exp $
!
!     simc: local potential generator
!           MINPACK version
!
    implicit none
    integer :: mesh
    real(8) :: rad(mesh),vin(mesh),parloc(4),zv,rc

    integer,parameter :: maxnum = 2000

    common /parusc/ rads,vins,zvs,wgt
    real*8 rads(maxnum),vins(maxnum),zvs,wgt(maxnum)

    real(8) :: lambda = 3.5d0
    real(8) :: pi

    integer :: k,num
    real(8) :: nxtsmp
    real(8) :: x(3),fvec(maxnum),fjac(maxnum,3)
    integer :: info
    real(8) :: x1ini,x2ini,x3ini
    parameter (x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0)
    real(8),parameter :: tol=1.0d-5
    real(8) :: rmax,vmax,vrzmin,smpstp
    integer :: nummin
    parameter(rmax=10.0d0, vmax=100.d0, vrzmin=3.0d-6)
    parameter(smpstp=0.2d0, nummin=6)
    integer :: ipvt(3)
    integer :: lwa
    parameter(lwa=5*3+maxnum)
    real(8) :: wa(lwa)

    pi = 4.0d0*atan(1.0d0)

    num=0
    nxtsmp=0.d0
    do k=1,mesh
       if ( rad(k)>nxtsmp ) then
          nxtsmp = nxtsmp + smpstp
          if ( abs(vin(k)) <= vmax ) then
!**   1.0d0*rc ---- 1.5d0*rc in original
             num=num+1
             if ( num > maxnum ) then
                write(6,*) 'simc: Too many sample points.'
                stop
             end if
             rads(num)=rad(k)
             vins(num)=vin(k)
             wgt(num)=1.d0-dexp(-(1.2d0*rad(k)/rc)**lambda)
          end if
          if ( (abs(vin(k)*rad(k)+zv)<vrzmin .or. rad(k)>rmax) &
               .and. num>nummin) exit
       end if
    end do
    zvs = zv
    x(1) = x1ini
    x(2) = x2ini
    x(3) = x3ini
    call lmder1(uscfit,num,3,x,fvec,fjac,maxnum,tol,info,ipvt,wa,lwa)
!    write(6,*) 'lmder1:',info
    if ( info==0 .or. info==4 .or. info==5 .or. info==6 .or. info==7 ) then
       write(6,*) 'simc: Not converged.'
       write(6,*) 'x(1) = ',x(1)
       write(6,*) 'x(2) = ',x(2)
       write(6,*) 'x(3) = ',x(3)
       write(6,'(A5,3A20)') 'k','rads(k)','vins(k)','wgt(k)'
       do k=1,num
          write(6,'(I5,3g20.12)')k,rads(k),vins(k),wgt(k)
       end do
       stop
    end if
    if ( x(2)<0.0d0 .or. x(3)<0.0d0 ) then
       write(6,*)'simc: illegally converged.'
       stop
    end if
    parloc(1) = x(1)
    parloc(2) = x(2)
    parloc(3) = 1.0d0 - x(1)
    parloc(4) = x(3)
    return
  END SUBROUTINE simc
!
!     fitting function for simc
!
  SUBROUTINE uscfit(m,n,x,fvec,fjac,ldfjac,iflag)
    implicit none
    integer :: m,n,ldfjac,iflag
    real(8) :: x(n),fvec(m),fjac(ldfjac,n)
    integer,parameter :: maxnum=2000
    common /parusc/ rad,vin,zv,wgt
    real(8) :: rad(maxnum),vin(maxnum),zv,wgt(maxnum)
    real(8) :: pi
    integer :: i

    pi = 4.0d0*atan(1.0d0)

    if ( x(2) < 0.0d0 ) x(2)=0.0d0
    if ( x(3) < 0.0d0 ) x(3)=0.0d0
    if ( iflag == 1 ) then
       do i=1,m
          fvec(i) = (- zv/rad(i)*(x(1)*bberf(sqrt(x(2))*rad(i)) &
               + (1.0d0-x(1))*bberf(sqrt(x(3))*rad(i))) - vin(i))*sqrt(wgt(i))
       end do
    else if ( iflag == 2 ) then
       do i=1,m
          fjac(i,1) = -zv/rad(i) &
               *(bberf(sqrt(x(2)*rad(i)))-bberf(sqrt(x(3)*rad(i)))) &
               *sqrt(wgt(i))
          fjac(i,2) = -zv/sqrt(pi*x(2))*x(1) &
               *exp(-x(2)*rad(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zv/sqrt(pi*x(3))*(1.0d0-x(1)) &
               *exp(-x(3)*rad(i)**2)*sqrt(wgt(i))
       end do
    else
       write(6,*) 'Error in vlfit: iflag must be 1 or 2.'
       stop
    end if
    return
  END SUBROUTINE uscfit


END MODULE simc_module
