! $Id: fft3d.f90,v 1.9 2001/01/22 07:01:04 yosimoto Exp $
module fft_3d
  integer,parameter,private:: nradix=4
  integer,dimension(nradix),parameter,private:: radix=(/3, 4, 2, 5/)
  integer,parameter,private:: nfmax = 30
  private:: fact3425,genw,fftv,bftv,fft_tp
  private:: fftrd2,bftrd2,fftrd3,bftrd3,fftrd4,bftrd4,fftrd5,bftrd5
  !
  type fft3_struct
     integer:: nx,ny,nz,nxyz,lx,ly,lz,lxyz
     integer:: nfx,nfy,nfz
     integer,dimension(nfmax):: facx,facy,facz
     real(kind=8),pointer,dimension(:,:):: wx,wy,wz
  end type fft3_struct
contains
!
  subroutine fact3425(n,nf,fac)
    implicit none
    integer,intent(in):: n
    integer,intent(out):: nf,fac(nfmax)
    !
    integer:: p,q,r,x,j
    !
    p = n
    fac(:) = 0
    !
    j = 1
    nf = 0
    x = radix(j)
    do while (p .ne. 1)
       q = p/x
       r = p-x*q
       if (r .eq. 0) then
          nf = nf + 1
          fac(nf) = x
          p = q
       else if (j .lt. nradix) then
          j=j+1
          x = radix(j)
       else
          write(6,*) n, ' can not be factorized with try radixes'
          stop
       end if
    end do
    return
  end subroutine fact3425

  subroutine genw(n,nf,fac,w)
    implicit none
    integer,intent(in):: n,nf
    integer:: fac(nfmax)
    real(kind=8),intent(out):: w(2,n)
    !
    real(kind=8),parameter:: tpi=2.0d0*3.14159265358979323846d0
    integer:: iw,darg,ir,rad,narg,k,p
    real(kind=8):: argh,argi,arg
    !
    iw = 1
    darg = n
    do ir = 1,nf
       argh = tpi/darg
       rad = fac(ir)
       narg = darg/rad
       do p = 0,narg-1
          do k = 1,rad-1
             argi = argh*k
             arg = argi*p
             w(1,iw) = cos(arg)
             w(2,iw) = -sin(arg)
             iw = iw + 1
          end do
       end do
       darg = darg/rad
    end do
    return
  end subroutine genw

  subroutine fftrd2(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,2),xi(nvk,np,2)
    real(kind=8),intent(out):: yr(nvk,2,np),yi(nvk,2,np)
    real(kind=8),intent(in):: w(2,np)
    !
    real(kind=8):: wr1,wi1,tr,ti
    integer:: p,ivk
    !
    do p = 1,np
       wr1 = w(1,p)
       wi1 = w(2,p)
       do ivk = 1,nvk
          yr(ivk,1,p) = xr(ivk,p,1) + xr(ivk,p,2)
          yi(ivk,1,p) = xi(ivk,p,1) + xi(ivk,p,2)
          !
          tr = xr(ivk,p,1) - xr(ivk,p,2)
          ti = xi(ivk,p,1) - xi(ivk,p,2)
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
       end do
    end do
    return
  end subroutine fftrd2
  !
  subroutine bftrd2(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,2),xi(nvk,np,2)
    real(kind=8),intent(out):: yr(nvk,2,np),yi(nvk,2,np)
    real(kind=8),intent(in):: w(2,np)
    !
    real(kind=8):: wr1,wi1,tr,ti
    integer:: p,ivk
    !
    do p = 1,np
       wr1 = w(1,p)
       wi1 = w(2,p)
       do ivk = 1,nvk
          yr(ivk,1,p) = xr(ivk,p,1) + xr(ivk,p,2)
          yi(ivk,1,p) = xi(ivk,p,1) + xi(ivk,p,2)
          !
          tr = xr(ivk,p,1) - xr(ivk,p,2)
          ti = xi(ivk,p,1) - xi(ivk,p,2)
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
       end do
    end do
    return
  end subroutine bftrd2
  !
  subroutine fftrd3(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,3),xi(nvk,np,3)
    real(kind=8),intent(out):: yr(nvk,3,np),yi(nvk,3,np)
    !  real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,2*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,3)
          ti1 = xi(ivk,p,2) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = sp3*(xi(ivk,p,2) - xi(ivk,p,3))
          ti3 = -sp3*(xr(ivk,p,2) - xr(ivk,p,3))
          yr(ivk,1,p) = xr(ivk,p,1) + tr1
          yi(ivk,1,p) = xi(ivk,p,1) + ti1
          tr = tr2 + tr3
          ti = ti2 + ti3
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr2 - tr3
          ti = ti2 - ti3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
       end do
    end do
    return
  end subroutine fftrd3
  !
  subroutine bftrd3(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,3),xi(nvk,np,3)
    real(kind=8),intent(out):: yr(nvk,3,np),yi(nvk,3,np)
    !real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,2*np)
    !
    real(kind=8),parameter:: sp3 = 0.866025403784438647d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,tr1,ti1,tr2,ti2,tr3,ti3,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,3)
          ti1 = xi(ivk,p,2) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - 0.5d0*tr1
          ti2 = xi(ivk,p,1) - 0.5d0*ti1
          tr3 = -sp3*(xi(ivk,p,2) - xi(ivk,p,3))
          ti3 = sp3*(xr(ivk,p,2) - xr(ivk,p,3))
          yr(ivk,1,p) = xr(ivk,p,1) + tr1
          yi(ivk,1,p) = xi(ivk,p,1) + ti1
          tr = tr2 + tr3
          ti = ti2 + ti3
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr2 - tr3
          ti = ti2 - ti3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
       end do
    end do
    return
  end subroutine bftrd3
  !
  subroutine fftrd4(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,4),xi(nvk,np,4)
    real(kind=8),intent(out):: yr(nvk,4,np),yi(nvk,4,np)
    !    real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,3*np)
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,3)
          ti1 = xi(ivk,p,1) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - xr(ivk,p,3)
          ti2 = xi(ivk,p,1) - xi(ivk,p,3)
          tr3 = xr(ivk,p,2) + xr(ivk,p,4)
          ti3 = xi(ivk,p,2) + xi(ivk,p,4)
          tr4 = xi(ivk,p,2) - xi(ivk,p,4)
          ti4 = xr(ivk,p,4) - xr(ivk,p,2)
          yr(ivk,1,p) = tr1 + tr3
          yi(ivk,1,p) = ti1 + ti3
          tr = tr2 + tr4
          ti = ti2 + ti4
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr1 - tr3
          ti = ti1 - ti3
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = tr2 - tr4
          ti = ti2 - ti4
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
       end do
    end do
    return
  end subroutine fftrd4
  !
  subroutine bftrd4(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,4),xi(nvk,np,4)
    real(kind=8),intent(out):: yr(nvk,4,np),yi(nvk,4,np)
    !    real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,3*np)
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,1) + xr(ivk,p,3)
          ti1 = xi(ivk,p,1) + xi(ivk,p,3)
          tr2 = xr(ivk,p,1) - xr(ivk,p,3)
          ti2 = xi(ivk,p,1) - xi(ivk,p,3)
          tr3 = xr(ivk,p,2) + xr(ivk,p,4)
          ti3 = xi(ivk,p,2) + xi(ivk,p,4)
          tr4 = xi(ivk,p,4) - xi(ivk,p,2)
          ti4 = xr(ivk,p,2) - xr(ivk,p,4)
          yr(ivk,1,p) = tr1 + tr3
          yi(ivk,1,p) = ti1 + ti3
          tr = tr2 + tr4
          ti = ti2 + ti4
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr1 - tr3
          ti = ti1 - ti3
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = tr2 - tr4
          ti = ti2 - ti4
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
       end do
    end do
    return
  end subroutine bftrd4
  !
  subroutine fftrd5(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,5),xi(nvk,np,5)
    real(kind=8),intent(out):: yr(nvk,5,np),yi(nvk,5,np)
    !real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,4*np)
    !
    real(kind=8),parameter:: sp25 = 0.951056516295153572d0
    real(kind=8),parameter:: sq54 = 0.559016994374947424d0
    real(kind=8),parameter:: ss = 0.618033988749894848d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5
    real(kind=8):: tr6,ti6,tr7,ti7,tr8,ti8,tr9,ti9,tr10,ti10,tr11,ti11
    real(kind=8):: tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,5)
          ti1 = xi(ivk,p,2) + xi(ivk,p,5)
          tr2 = xr(ivk,p,3) + xr(ivk,p,4)
          ti2 = xi(ivk,p,3) + xi(ivk,p,4)
          tr3 = sp25*(xr(ivk,p,2) - xr(ivk,p,5))
          ti3 = sp25*(xi(ivk,p,2) - xi(ivk,p,5))
          tr4 = sp25*(xr(ivk,p,3) - xr(ivk,p,4))
          ti4 = sp25*(xi(ivk,p,3) - xi(ivk,p,4))
          tr5 = tr1 + tr2
          ti5 = ti1 + ti2
          tr6 = sq54*(tr1 - tr2)
          ti6 = sq54*(ti1 - ti2)
          tr7 = xr(ivk,p,1) - 0.25d0*tr5
          ti7 = xi(ivk,p,1) - 0.25d0*ti5
          tr8 = tr7 + tr6
          ti8 = ti7 + ti6
          tr9 = tr7 - tr6
          ti9 = ti7 - ti6
          tr10 = ti3 + ss*ti4
          ti10 = - tr3 - ss*tr4
          tr11 = ss*ti3 - ti4
          ti11 = - ss*tr3 + tr4
          yr(ivk,1,p) = xr(ivk,p,1) + tr5
          yi(ivk,1,p) = xi(ivk,p,1) + ti5
          tr = tr8 + tr10
          ti = ti8 + ti10
          yr(ivk,2,p) = wr1*tr - wi1*ti
          yi(ivk,2,p) = wr1*ti + wi1*tr
          tr = tr9 + tr11
          ti = ti9 + ti11
          yr(ivk,3,p) = wr2*tr - wi2*ti
          yi(ivk,3,p) = wr2*ti + wi2*tr
          tr = tr9 - tr11
          ti = ti9 - ti11
          yr(ivk,4,p) = wr3*tr - wi3*ti
          yi(ivk,4,p) = wr3*ti + wi3*tr
          tr = tr8 - tr10
          ti = ti8 - ti10
          yr(ivk,5,p) = wr4*tr - wi4*ti
          yi(ivk,5,p) = wr4*ti + wi4*tr
       end do
    end do
    return
  end subroutine fftrd5
  !
  !
  subroutine bftrd5(nvk,np,w,xr,xi,yr,yi)
    implicit none
    integer,intent(in):: nvk,np
    real(kind=8),intent(in):: xr(nvk,np,5),xi(nvk,np,5)
    real(kind=8),intent(out):: yr(nvk,5,np),yi(nvk,5,np)
    !real(kind=8),intent(in):: w(2,np)
    real(kind=8),intent(in):: w(2,4*np)
    !
    real(kind=8),parameter:: sp25 = 0.951056516295153572d0
    real(kind=8),parameter:: sq54 = 0.559016994374947424d0
    real(kind=8),parameter:: ss = 0.618033988749894848d0
    !
    real(kind=8):: wr1,wi1,wr2,wi2,wr3,wi3,wr4,wi4
    real(kind=8):: tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5
    real(kind=8):: tr6,ti6,tr7,ti7,tr8,ti8,tr9,ti9,tr10,ti10,tr11,ti11
    real(kind=8):: tr,ti
    integer:: p,ivk,iw
    !
    iw = 1
    do p = 1,np
       wr1 = w(1,iw)
       wi1 = w(2,iw)
       iw = iw+1
       wr2 = w(1,iw)
       wi2 = w(2,iw)
       iw = iw+1
       wr3 = w(1,iw)
       wi3 = w(2,iw)
       iw = iw+1
       wr4 = w(1,iw)
       wi4 = w(2,iw)
       iw = iw+1
       do ivk = 1,nvk
          tr1 = xr(ivk,p,2) + xr(ivk,p,5)
          ti1 = xi(ivk,p,2) + xi(ivk,p,5)
          tr2 = xr(ivk,p,3) + xr(ivk,p,4)
          ti2 = xi(ivk,p,3) + xi(ivk,p,4)
          tr3 = sp25*(xr(ivk,p,2) - xr(ivk,p,5))
          ti3 = sp25*(xi(ivk,p,2) - xi(ivk,p,5))
          tr4 = sp25*(xr(ivk,p,3) - xr(ivk,p,4))
          ti4 = sp25*(xi(ivk,p,3) - xi(ivk,p,4))
          tr5 = tr1 + tr2
          ti5 = ti1 + ti2
          tr6 = sq54*(tr1 - tr2)
          ti6 = sq54*(ti1 - ti2)
          tr7 = xr(ivk,p,1) - 0.25d0*tr5
          ti7 = xi(ivk,p,1) - 0.25d0*ti5
          tr8 = tr7 + tr6
          ti8 = ti7 + ti6
          tr9 = tr7 - tr6
          ti9 = ti7 - ti6
          tr10 = - ti3 - ss*ti4
          ti10 = tr3 + ss*tr4
          tr11 = - ss*ti3 + ti4
          ti11 = ss*tr3 - tr4
          yr(ivk,1,p) = xr(ivk,p,1) + tr5
          yi(ivk,1,p) = xi(ivk,p,1) + ti5
          tr = tr8 + tr10
          ti = ti8 + ti10
          yr(ivk,2,p) = wr1*tr + wi1*ti
          yi(ivk,2,p) = wr1*ti - wi1*tr
          tr = tr9 + tr11
          ti = ti9 + ti11
          yr(ivk,3,p) = wr2*tr + wi2*ti
          yi(ivk,3,p) = wr2*ti - wi2*tr
          tr = tr9 - tr11
          ti = ti9 - ti11
          yr(ivk,4,p) = wr3*tr + wi3*ti
          yi(ivk,4,p) = wr3*ti - wi3*tr
          tr = tr8 - tr10
          ti = ti8 - ti10
          yr(ivk,5,p) = wr4*tr + wi4*ti
          yi(ivk,5,p) = wr4*ti - wi4*tr
       end do
    end do
    return
  end subroutine bftrd5
  subroutine fftv(nv,n,nf,fac,w,lxyz,dr,di,er,ei)
    implicit none
    integer,intent(in):: nv,n,lxyz
    integer,intent(in):: nf,fac(nfmax)
    real(kind=8),intent(in):: w(2,n)
    real(kind=8),dimension(lxyz):: dr,di,er,ei
    !
    integer:: iw,nvk,np,uif,if
    !
    iw = 1
    nvk = nv
    np = n/fac(1)
    uif = mod(nf,2)
    do if = 1,nf-uif,2
       if (fac(if) .eq. 2) then
          call fftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 3) then
          call fftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 4) then
          call fftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 5) then
          call fftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       iw=iw+np*(fac(if)-1)
       nvk=nvk*fac(if)
       np=np/fac(if+1)
       if (fac(if+1) .eq. 2) then
          call fftrd2(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 3) then
          call fftrd3(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 4) then
          call fftrd4(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 5) then
          call fftrd5(nvk,np,w(1,iw),er,ei,dr,di)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       iw=iw+np*(fac(if+1)-1)
       nvk=nvk*fac(if+1)
       if (if+2 .le. nf) then
          np=np/fac(if+2)
       end if
    end do
    if (uif .eq. 1) then
       if (fac(nf) .eq. 2) then
          call fftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 3) then
          call fftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 4) then
          call fftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 5) then
          call fftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       dr = er
       di = ei
    end if

    return
  end subroutine fftv
  !
  subroutine bftv(nv,n,nf,fac,w,lxyz,dr,di,er,ei)
    implicit none
    integer,intent(in):: nv,n,lxyz
    integer,intent(in):: nf,fac(nfmax)
    real(kind=8),intent(in):: w(2,n)
    real(kind=8),dimension(lxyz):: dr,di,er,ei
    !
    integer:: iw,nvk,np,uif,if
    !
    iw = 1
    nvk = nv
    np = n/fac(1)
    uif = mod(nf,2)
    do if = 1,nf-uif,2
       if (fac(if) .eq. 2) then
          call bftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 3) then
          call bftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 4) then
          call bftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(if) .eq. 5) then
          call bftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       iw=iw+np*(fac(if)-1)
       nvk=nvk*fac(if)
       np=np/fac(if+1)
       if (fac(if+1) .eq. 2) then
          call bftrd2(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 3) then
          call bftrd3(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 4) then
          call bftrd4(nvk,np,w(1,iw),er,ei,dr,di)
       else if (fac(if+1) .eq. 5) then
          call bftrd5(nvk,np,w(1,iw),er,ei,dr,di)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       iw=iw+np*(fac(if+1)-1)
       nvk=nvk*fac(if+1)
       if (if+2 .le. nf) then
          np=np/fac(if+2)
       end if
    end do
    if (uif .eq. 1) then
       if (fac(nf) .eq. 2) then
          call bftrd2(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 3) then
          call bftrd3(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 4) then
          call bftrd4(nvk,np,w(1,iw),dr,di,er,ei)
       else if (fac(nf) .eq. 5) then
          call bftrd5(nvk,np,w(1,iw),dr,di,er,ei)
       else
          write(6,*)'unsupported factor'
          stop
       end if
       dr = er
       di = ei
    end if

    return
  end subroutine bftv
  !
  subroutine fft_pad_zero(nx,ny,nz,lx,ly,lz,dat)
    implicit none
    integer,intent(in):: nx,ny,nz,lx,ly,lz
    real(kind=8):: dat(lx,ly,lz,2)
    !
    integer:: j,k
    do k=1,2
       do j=nz+1,lz
          dat(:,:,j,k) = 0.0d0
       end do
       do j=nx+1,lx
          dat(j,:,:,k) = 0.0d0
       end do
       do j=ny+1,ly
          dat(:,j,:,k) = 0.0d0
       end do
    end do
    return
  end subroutine fft_pad_zero
  !
  subroutine fft_tp(lxy,lz,src,dst)
    implicit none
    integer,intent(in):: lxy,lz
    real(kind=8),intent(in):: src(lxy,lz,2)
    real(kind=8),intent(out):: dst(lz,lxy,2)
    !
    integer:: i,j
! directives for Origin2xxx
!DIR$ UNROLL(8)
    do i=1,lz
!DIR$ UNROLL(8)
       do j=1,lxy
          dst(i,j,1) = src(j,i,1)
       end do
    end do
!DIR$ UNROLL(8)
    do i=1,lz
!DIR$ UNROLL(8)
       do j=1,lxy
          dst(i,j,2) = src(j,i,2)
       end do
    end do
    return
  end subroutine fft_tp
  !
  subroutine fft3_init(mx,my,mz,lmx,lmy,lmz,fs)
    implicit none
    integer,intent(in):: mx,my,mz,lmx,lmy,lmz
    type (fft3_struct):: fs
    !
    fs%nx = mx
    fs%ny = my
    fs%nz = mz
    fs%nxyz = mx*my*mz
    fs%lx = lmx
    fs%ly = lmy
    fs%lz = lmz
    fs%lxyz = lmx*lmy*lmz
    call fact3425(mx,fs%nfx,fs%facx)
    call fact3425(my,fs%nfy,fs%facy)
    call fact3425(mz,fs%nfz,fs%facz)
    allocate(fs%wx(2,mx),fs%wy(2,my),fs%wz(2,mz))
    call genw(mx,fs%nfx,fs%facx,fs%wx)
    call genw(my,fs%nfy,fs%facy,fs%wy)
    call genw(mz,fs%nfz,fs%facz,fs%wz)
    return
  end subroutine fft3_init

  subroutine fft3_done(fs)
    implicit none
    type (fft3_struct):: fs
    !
    deallocate(fs%wx,fs%wy,fs%wz)
    return
  end subroutine fft3_done
#ifdef FFTVEC
  ! FFT routine for Vector machine
  subroutine fft3_fw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    !jy 2003-08-27
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,dat)
    !
    call fftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))
    !
    call fft_tp(fs%lx*fs%ly, fs%lz, dat, tmp)
    !
    call fftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
         fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2))
    !
    call fft_tp(fs%lz*fs%lx, fs%ly, tmp, dat)
    !
    call fftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))
    !
    call fft_tp(fs%ly*fs%lz, fs%lx, dat, tmp)
    !
    dat = tmp/fs%nxyz
    return
  end subroutine fft3_fw

  subroutine fft3_bw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    !jy 2003-08-27
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,dat)
    !
    call bftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))
    !
    call fft_tp(fs%lx*fs%ly, fs%lz, dat, tmp)
    !
    call bftv(fs%lz*fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, &
         fs%lxyz, tmp(1,1), tmp(1,2), dat(1,1), dat(1,2))
    !
    call fft_tp(fs%lz*fs%lx, fs%ly, tmp, dat)
    !
    call bftv(fs%ly*fs%lz, fs%nx, fs%nfx, fs%facx, fs%wx, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))
    !
    call fft_tp(fs%ly*fs%lz, fs%lx, dat, tmp)
    dat = tmp
    return
  end subroutine fft3_bw
#else
  ! FFT routine for Scalar machine: architecture using CACHE
  subroutine fft3_fw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: iv,lxy
    real(kind=8):: scl
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    !jy 2003-08-27
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,dat)
    !
    do iv=0, fs%ly*fs%nz - 1
       call fftv(1, fs%nx, fs%nfx, fs%facx, fs%wx, fs%lx, &
            dat(iv*fs%lx+1,1), dat(iv*fs%lx+1,2), &
            tmp(iv*fs%lx+1,1), tmp(iv*fs%lx+1,2))
    end do
    !
    lxy = fs%lx*fs%ly
    do iv=0, fs%nz - 1
       call fftv(fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, lxy, &
            dat(iv*lxy+1,1), dat(iv*lxy+1,2), &
            tmp(iv*lxy+1,1), tmp(iv*lxy+1,2))
    end do
    !
    call fftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))

    scl = 1.0d0/fs%nxyz
    dat = scl*dat
    return
  end subroutine fft3_fw

  subroutine fft3_bw(fs,dat,tmp)
    implicit none
    type (fft3_struct):: fs
    real(kind=8):: dat(fs%lxyz,2)
    real(kind=8):: tmp(fs%lxyz,2)
    integer:: iv,lxy
    !
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,tmp)
    !jy 2003-08-27
    call fft_pad_zero(fs%nx,fs%ny,fs%nz,fs%lx,fs%ly,fs%lz,dat)
    !
    do iv=0, fs%ly*fs%nz - 1
       call bftv(1, fs%nx, fs%nfx, fs%facx, fs%wx, fs%lx, &
            dat(iv*fs%lx+1,1), dat(iv*fs%lx+1,2), &
            tmp(iv*fs%lx+1,1), tmp(iv*fs%lx+1,2))
    end do
    !
    lxy = fs%lx*fs%ly
    do iv=0, fs%nz - 1
       call bftv(fs%lx, fs%ny, fs%nfy, fs%facy, fs%wy, lxy, &
            dat(iv*lxy+1,1), dat(iv*lxy+1,2), &
            tmp(iv*lxy+1,1), tmp(iv*lxy+1,2))
    end do
    !
    call bftv(fs%lx*fs%ly, fs%nz, fs%nfz, fs%facz, fs%wz, &
         fs%lxyz, dat(1,1), dat(1,2), tmp(1,1), tmp(1,2))
    return
  end subroutine fft3_bw
#endif
end module fft_3d
