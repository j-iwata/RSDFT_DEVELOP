!--------1---------2---------3---------4---------5---------6---------7--
!======================================================= Pseudopotential

      SUBROUTINE init_ps_mol
      use global_variables
      implicit none

      integer :: i,j,k,L,ik,a,m,lm,iorb,g,morb,lm0
      integer :: NRc,m1,m2,MMr,lma,m0,mm,iloc(1)
      integer,allocatable :: Jtmp(:,:,:),loc(:)
      real(8) :: x,y,z,r,c1,c2,uVrl,r1
      real(8) :: qc,sum0,maxerr,c,lambda0,eta0
      real(8) :: p1,p2,p3,p4,sb0x,sb1x,sb0y,sb1y,sb2,sb3
      real(8) :: Rc,dy0,dy,y0,pi4,const,r2
      real(8),parameter :: dr=2.d-3
      real(8),parameter :: eps=1.d-20
      real(8),parameter :: ep0=1.d-8
      real(8),parameter :: ep1=1.d-10
      real(8),allocatable :: tmp1(:),tmp(:),vork(:),wtmp(:)
      real(8),allocatable :: dvork(:,:,:),vshort(:),vlong(:)
      real(8),allocatable :: dvshort(:),dvlong(:)
      real(8) :: fac,drtmp
      real(8) :: mem,memax
      logical :: iflag
      integer,parameter :: nmsk=201
      real(8) :: xm(0:nmsk),maskr(0:nmsk),dxm
      real(8),allocatable :: wm(:,:,:)
      integer :: iu
      real(8) :: ctime0,ctime1,etime0,etime1
#ifdef TEST
      integer :: nqtmp,nrtmp
      real(8) :: q,dq,sum1
      real(8),allocatable :: qtmp(:),rtmp(:),dtmp(:)
#endif

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf
         real(8),intent(IN) :: x
         END FUNCTION bberf
      END INTERFACE

      write(*,'(a60," init_ps")') repeat("-",60)

      if ( pselect<1 .or. 3<pselect ) then
         write(*,*) "pselect=4 is not supported"
         call stop_program
      end if

      mem   = 0.d0
      memax = 0.d0
      pi4   = 16.d0*atan(1.d0)
      iu    = 70

      lo(:,:)      = 0
      rad(:,:)     = 0.d0
      rab(:,:)     = 0.d0
      vql(:,:)     = 0.d0
      viod(:,:,:)  = 0.d0
      dvql(:,:)    = 0.d0
      dviod(:,:,:) = 0.d0
      inorm(:,:)   = 0
      anorm(:,:)   = 0.d0
      cdc(:,:)     = 0.d0
      cdd(:,:)     = 0.d0
      Mr(:)        = 0
      norb(:)      = 0
      Rps(:,:)     = 0.d0
      NRps(:,:)    = 0
      Lref(:)      = 0
      parloc(:,:)    = 0.d0
      Rcloc(:)     = 0.d0
      NRcloc(:)    = 0.d0

!
! --- Read atomic pseudopotential data ---
! ( TM, KY and PSV data format is available. )
! ( PP data is read from UNIT=34,35,... for the 1st,2nd,... element. )
!

      do ik=1,MKI
         g=33+ik
         select case(ippform(ik))
         case default
            open(g,file=file_ps(ik),status='old')
            call KY_format( g,nrmax,lpsmax,lo(:,ik) &
     &           ,rad(:,ik),rab(:,ik),vql(:,ik),viod(:,:,ik) &
     &           ,dvql(:,ik),dviod(:,:,ik) &
     &           ,inorm(:,ik),anorm(:,ik),cdc(:,ik),cdd(:,ik),Mr(ik) &
     &           ,norb(ik),Zps(ik),Rps(:,ik),NRps(:,ik),Lref(ik) )
         case(1)
            open(g,file=file_ps(ik),form='unformatted',status='old')
            call TM_format( g,nrmax,lpsmax,lo(:,ik) &
     &           ,rad(:,ik),rab(:,ik),vql(:,ik),viod(:,:,ik) &
     &           ,dvql(:,ik),dviod(:,:,ik),inorm(:,ik),anorm(:,ik) &
     &           ,cdc(:,ik),cdd(:,ik),Mr(ik),norb(ik),Zps(ik),Rps(:,ik) &
     &           ,NRps(:,ik),Lref(ik) )
         case(2)
            open(g,file=file_ps(ik),form='formatted',status='old')
            call PSV_format( g,nrmax,lpsmax,lo(:,ik) &
     &           ,rad(:,ik),rab(:,ik) &
     &           ,vql(:,ik),viod(:,:,ik) &
     &           ,dvql(:,ik),dviod(:,:,ik) &
     &           ,inorm(:,ik),anorm(:,ik),cdc(:,ik),cdd(:,ik),Mr(ik) &
     &           ,norb(ik),Zps(ik),Rps(:,ik),NRps(:,ik),Lref(ik) &
     &           ,parloc(1,ik))
         end select
         close(g)
      end do


!      do ik=1,MKI
!         rewind iu
!         write(*,*) "--- write vql,viod (read) ---",ik,iu
!         MMr=Mr(ik)
!         morb=norb(ik)
!         do i=1,MMr
!            write(iu,'(1x,g16.8,1x,5g16.8)') rad(i,ik),vql(i,ik),viod(i,1:morb,ik)
!         end do
!         iu=iu+1
!      end do
!      stop

!
! Cutoff
!
      qc=qf*qcfac
      if (DISP_SWITCH) then
         write(*,*) "qf,qc,qcfac =",qf,qc,qcfac
      end if
      if ( qc<=0.d0 ) qc=qf

! ----------------------------------------------------------------------
! ------------------------------ Nonlocal ------------------------------
! ---------------------------------------------------------------------- 

      select case(pselect)
      case default

         write(*,*) "*** Potential filtering (pselect=1) ***"

         NRps0(:,:) = NRps(:,:)
         NRps1(:,:) = NRps(:,:)
         Rps0(:,:)  = Rps(:,:)
         Rps1(:,:)  = Rps(:,:)
         rad1(:,:)  = rad(:,:)

      case(2)

         write(*,*) "*** Potential filtering (pselect=2) ***"

         MMr  = maxval(Mr)
         morb = maxval(norb)

!- allocate ---------------------------------------------------
         allocate( wm(MMr,morb,MKI) ) ; wm=0.d0
         mem=mem+bdreal*(MMr*morb*MKI) ; memax=max(mem,memax)
!--------------------------------------------------------------

         NRps0(:,:) = NRps(:,:)
         Rps0(:,:)  = Rps(:,:)

         do ik=1,MKI
            MMr=Mr(ik)
            do iorb=1,norb(ik)
               Rc=Rps(iorb,ik)*rcfac
               iloc=minloc( abs(rad(1:MMr,ik)-Rc) )
               NRc=iloc(1) ; if ( rad(NRc,ik)<Rc ) NRc=NRc+1
               if ( NRc>MMr ) then
                  write(*,*) "NRc,MMr=",NRc,MMr
                  stop
               end if
               NRps(iorb,ik) = NRc
               Rps(iorb,ik)  = rad(NRc,ik)
            end do
         end do

         write(*,*) "rcfac =",rcfac
         write(*,'(1x,a3,a5,3a6,2a18)') "ik","iorb","NRps0","NRps","Mr","Rps0","Rps"
         do ik=1,MKI
            do iorb=1,norb(ik)
               write(*,'(1x,i3,i4,1x,3i6,2f18.14)') &
                    ik,iorb,NRps0(iorb,ik),NRps(iorb,ik),Mr(ik),Rps0(iorb,ik),Rps(iorb,ik)
            end do
         end do

         do ik=1,MKI
         do iorb=1,norb(ik)
            NRc = NRps(iorb,ik)
            Rc  = Rps(iorb,ik)
            call makemaskf
            c=1.d0/Rc
            maxerr=0.d0
            do i=1,NRc
               x=rad(i,ik)*c
               if ( x<=dxm ) then
                  y0=1.d0 ; dy0=0.d0
               else
                  m0=int(x/dxm)
                  dy0=1.d10
                  do m=1,20
                     m1=max(m0-m,1)       ; mm=m1-(m0-m)
                     m2=min(m0+m+mm,nmsk) ; mm=(m0+m+mm)-m2
                     call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                     if ( abs(dy)<dy0 ) then
                        y0=y ; dy0=abs(dy)
!                        if ( dy0<ep0 ) exit
                     end if
                  end do
               end if
!               if ( y0>1.d0 ) y0=1.d0
!               if ( y0<0.d0 ) y0=0.d0
               wm(i,iorb,ik)=y0
               maxerr=max(maxerr,dy0)
            end do
            write(*,*) "err(maskf)=",maxerr
         end do
         end do

! --- Potential filtering (pselect=2) ---

         NRps1(:,:) = NRps(:,:)
         Rps1(:,:)  = Rps(:,:)

         do ik=1,MKI
         do iorb=1,norb(ik)
            NRps(iorb,ik)=Rps(iorb,ik)/dr
            write(*,*) "regular mesh (ik,iorb)",ik,iorb
            if ( NRps(iorb,ik)>nrmax ) stop
         end do
         end do
         MMr=max( maxval(Mr),maxval(NRps) )
         do ik=1,MKI
            do i=1,MMr
               rad1(i,ik)=(i-1)*dr
            end do
         end do

!--------------------- MEMO ---------------------------!
! NRps0, Rps0, rad  ---> original (log mesh)           !
! NRps1, Rps1, rad  ---> enlarged (rcfac,log mesh)     !
! NRps , Rps , rad1 ---> enlarged (rcfac,regular mesh) !
!--------------------- MEMO ---------------------------!


!- allocate ---------------------------------------------------------
         NRc=maxval(NRps0)
         allocate( vork(NRc) ) ; vork=0.d0 ; mem=mem+bdreal*NRc*2.d0
         allocate( tmp(NRc)  ) ; tmp=0.d0  ; mem=mem+bdreal*NRc*2.d0
         memax=max(mem,memax)
!--------------------------------------------------------------------

         const=2.d0/Pi

         do ik=1,MKI
         do iorb=1,norb(ik)

            L   = lo(iorb,ik)
            NRc = NRps0(iorb,ik)

            do j=1,NRc
               vork(j)=rad(j,ik)*viod(j,iorb,ik)*rab(j,ik)/wm(j,iorb,ik)
            end do

            do i=1,NRps(iorb,ik)
               r=rad1(i,ik)
               select case(L)
               case(0)
                  if ( r==0.d0 ) then
                     r1=rad(1,ik)
                     if ( r1==0.d0 ) then
                        tmp(1)=qc*qc*qc/3.d0
                     else
                        tmp(1)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
                     end if
                     do j=2,NRc
                        r1=rad(j,ik)
                        tmp(j)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
                     end do
                  else
                     do j=1,NRc
                        r1=rad(j,ik)
                        if ( r1==0.d0 ) then
                           tmp(j)=sin(qc*r)/(r*r*r)-qc*cos(qc*r)/(r*r)
                        else if ( r1==r ) then
                           tmp(j)=(2*qc*r-sin(2.d0*qc*r))/(4*r*r*r)
                        else
                           tmp(j)=( sin(qc*(r-r1))/(r-r1)-sin(qc*(r+r1))/(r+r1) )/(2.d0*r*r1)
                        end if
                     end do
                  end if
               case(1)
                  if ( r==0.d0 ) then
                     viod(i,iorb,ik)=0.d0
                     cycle
                  else
                     do j=1,NRc
                        r1=rad(j,ik)
                        if ( r1==0.d0 ) then
                           tmp(j)=0.d0
                        else if ( r1==r ) then
                           sb0x=sin(qc*r)/(qc*r)
                           sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                           tmp(j)=(2*qc*r-sin(2.d0*qc*r))/(4*r*r*r)-qc*qc*sb0x*sb1x/r
                        else
                           sb0x=sin(qc*r)/(qc*r)
                           sb0y=sin(qc*r1)/(qc*r1)
                           sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                           sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                           tmp(j)=( r1*sb0y*sb1x-r*sb0x*sb1y )*qc*qc/(r*r-r1*r1)
                        end if
                     end do
                  end if
               case(2)
                  if ( r==0.d0 ) then
                     viod(i,iorb,ik)=0.d0
                     cycle
                  else
                     do j=1,NRc
                        r1=rad(j,ik)
                        if ( r1==0.d0 ) then
                           tmp(j)=0.d0
                        else if ( r1==r ) then
                           sb1x=sin(qc*r)/(qc*qc*r*r)-cos(qc*r)/(qc*r)
                           tmp(j)=(2.d0*qc*r-sin(2*qc*r))/(4.d0*r*r*r)-3.d0*qc*sb1x*sb1x/(r*r)
                        else
                           sb0x=sin(qc*r)/(qc*r)
                           sb0y=sin(qc*r1)/(qc*r1)
                           sb1x=sb0x/(qc*r)-cos(qc*r)/(qc*r)
                           sb1y=sb0y/(qc*r1)-cos(qc*r1)/(qc*r1)
                           tmp(j)=( r*sb0y*sb1x-r1*sb0x*sb1y )*qc*qc/(r*r-r1*r1)-3.d0*qc/(r*r1)*sb1x*sb1y
                        end if
                     end do
                  end if
               case default
                  write(*,*) "PP for L>2 is not implemented."
                  stop
               end select
               tmp(1:NRc)=tmp(1:NRc)*vork(1:NRc)
               call simp(tmp(1:NRc),sum0,2)
               viod(i,iorb,ik)=sum0*const
            end do

            call makemaskf
            c=1.d0/Rps(iorb,ik)
            maxerr=0.d0
            do i=1,NRps(iorb,ik)
               x=rad1(i,ik)*c
               if ( x<=dxm ) then
                  y0=1.d0 ; dy0=0.d0
               else
                  m0=int(x/dxm)
                  dy0=1.d10
                  do m=1,20
                     m1=max(m0-m,1)       ; mm=m1-(m0-m)
                     m2=min(m0+m+mm,nmsk) ; mm=(m0+m+mm)-m2
                     call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                     if ( abs(dy)<dy0 ) then
                        y0=y ; dy0=abs(dy)
!c                        if ( abs(dy)<ep0 ) exit
                     end if
                  end do
               end if
!c               if ( y0>1.d0 ) y0=1.d0
!c               if ( y0<0.d0 ) y0=0.d0
               maxerr=max(maxerr,dy0)
               viod(i,iorb,ik)=y0*viod(i,iorb,ik)
            end do
            write(*,*) "err(maskf2)=",maxerr

         end do ! iorb
         end do ! ik

!- deallocate -------------------------------
         mem=mem-bdreal*size(vork)
         mem=mem-bdreal*size(tmp)
         mem=mem-bdreal*size(wm)
         deallocate( vork,tmp,wm )
!--------------------------------------------

!         do ik=1,MKI
!            rewind iu
!            write(*,*) "masked,regular mesh (iu,ik)",iu,ik
!            MMr=Mr(ik)
!            morb=norb(ik)
!            do i=1,MMr
!               write(iu,'(1x,g16.8,1x,5g16.7)') rad1(i,ik),viod(i,1:morb,ik)
!            end do
!            iu=iu+1
!         end do

      end select ! pselect

!
! --- derivative of viod() ---
!

      dviod(:,:,:)=0.d0

      do ik=1,MKI
      do iorb=1,norb(ik)
         L=lo(iorb,ik)
         NRc=NRps(iorb,ik)
         maxerr=0.d0
         do i=1,NRc
            dy0=1.d10
            do m=1,20
               m1=max(i-m,1)
               m2=min(i+m,NRc)
               call dpolint( rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,rad1(i,ik),y,dy )
               if ( abs(dy)<dy0 ) then
                  y0=y ; dy0=abs(dy)
!c                  if ( dy0<ep0 ) exit
               end if
            end do
            dviod(i,iorb,ik)=y0
            maxerr=max(maxerr,dy0)
         end do
         write(*,*) "err(dviod)=",maxerr
      end do
      end do

!      do ik=1,MKI
!         rewind iu
!         write(*,*) "viod,dviod(iu,ik,pselect)",iu,ik,pselect
!         MMr=Mr(ik)
!         morb=norb(ik)
!         do i=1,MMr
!            write(iu,'(1x,g16.8,1x,6g16.8)') rad1(i,ik),(viod(i,iorb,ik),dviod(i,iorb,ik),iorb=1,morb)
!         end do
!         iu=iu+1
!      end do

!*
!* "dviod" is modified for force calculation.
!*
      NRc = maxval(NRps)
      lm = 0
      do ik=1,MKI
         lm0=0
         do iorb=1,norb(ik)
            lm0=lm0+3 ; if ( lo(iorb,ik)==0 ) lm0=lm0-2
         end do
         lm=max(lm,lm0)
      end do

      if ( lm>3*lpsmax+1 ) then
         write(*,*) "lpsmax=",lpsmax
         write(*,*) "lm=",lm
         write(*,*) "maxval(lo)=",maxval(lo)
         stop "lpsmax is small."
      end if

!- allocate -------------------------------------------
      allocate( dvork(NRc,lm,MKI) ) ; dvork=0.d0
      mem=mem+bdreal*NRc*lm*MKI ; memax=max(mem,memax)
!------------------------------------------------------

      do ik=1,MKI
         lm=0
         do iorb=1,norb(ik)
            L=lo(iorb,ik)
            NRc=NRps(iorb,ik)
            do J=abs(L-1),L+1
               lm=lm+1
               if ( pselect==1 ) then
                  const=0.5d0*(L*(L+1)-J*(J+1))
                  do i=1,NRc
                     dvork(i,lm,ik)=rad1(i,ik)*dviod(i,iorb,ik)+const*viod(i,iorb,ik)
                  end do
               else if ( pselect==2 ) then
                  const=0.5d0*(2.d0+L*(L+1)-J*(J+1))
                  do i=1,NRc
                     dvork(i,lm,ik)=rad1(i,ik)**2*dviod(i,iorb,ik)+const*rad1(i,ik)*viod(i,iorb,ik)
                  end do
               end if
            end do
         end do
      end do
      const=sqrt(pi4/3.d0)
      do ik=1,MKI
         lm=0
         do iorb=1,norb(ik)
            L=lo(iorb,ik)
            NRc=NRps(iorb,ik)
            do J=abs(L-1),L+1
               lm=lm+1
               do i=1,NRc
                  dviod(i,lm,ik)=const*dvork(i,lm,ik)
               end do
            end do
         end do
      end do

!- deallocate -----------------------------------------
      mem=mem-bdreal*size(dvork) ; deallocate( dvork )
!------------------------------------------------------


! ---------------------------------------------------------------------
! -------------------------- Local potential --------------------------
! ---------------------------------------------------------------------

      write(*,*) "--- Local Potential ---"

!- allocate --------------------------------------------------
      MMr=max( maxval(Mr),maxval(NRps) )
      allocate( vshort(MMr)  ) ; vshort=0.d0
      allocate( vlong(MMr)   ) ; vlong=0.d0
      allocate( tmp1(MMr)    ) ; tmp1=0.d0
      allocate( wtmp(MMr)    ) ; wtmp=0.d0
      allocate( dvshort(MMr) ) ; dvshort=0.d0
      allocate( dvlong(MMr)  ) ; dvlong=0.d0
      allocate( tmp(MMr)     ) ; tmp=0.d0
      allocate( vork(MMr)    ) ; vork=0.d0
      mem=mem+bdreal*MMr*8 ; memax=max(mem,memax)
!-------------------------------------------------------------

      vqls(:,:)  = 0.d0
      dvql(:,:)  = 0.d0
      dvqls(:,:) = 0.d0

      c1=2.d0/sqrt(Pi)
      c2=2.d0/3.d0*c1

      do ik=1,MKI

         MMr        = Mr(ik)
         vlong(:)   = 0.d0
         vshort(:)  = 0.d0
         dvlong(:)  = 0.d0
         dvshort(:) = 0.d0

!         NRc=0
!         do iorb=1,norb(ik)
!            NRc=max( NRc, NRps0(iorb,ik) )
!         end do
         Rc=0.d0
         do iorb=1,norb(ik)
            Rc=max( Rc, Rps0(iorb,ik) )
         end do
!         Rc=1.5d0*Rc
         if ( Rc<1.d-10 ) Rc=5.d0

         do i=1,MMr
            if ( rad(i,ik)>Rc ) then
               NRc=i
               exit
            end if
         end do
         Rc=rad(NRc,ik)

         write(*,*) "ik,NRc,Rc=",ik,NRc,Rc

         call simc(rad(1,ik),vql(1,ik),Rc,Zps(ik),parloc(1,ik),MMr)
         p1=-Zps(ik)*parloc(1,ik) ; p2=sqrt(parloc(2,ik))
         p3=-Zps(ik)*parloc(3,ik) ; p4=sqrt(parloc(4,ik))
         write(*,*) "p1,p2,p3,p4 =",parloc(1,ik),parloc(2,ik)
         write(*,*) "             ",parloc(3,ik),parloc(4,ik)

! --- long range part of the local potential ---

         do i=1,MMr
            r=rad(i,ik)
            if ( r<1.d-9 ) then
               vlong(i)=c1*(p1*p2+p3*p4)
            else
               vlong(i)=( p1*bberf(p2*r)+p3*bberf(p4*r) )/r
            end if
         end do

! --- Derivative (long) ---

         do i=1,MMr
            r=rad(i,ik)
            if ( r<1.d-9 ) then
               dvlong(i)=-c2*r*(p1*p2*p2*p2+p3*p4*p4*p4)
            else
               dvlong(i)=-( p1*bberf(p2*r) + p3*bberf(p4*r) )/(r*r) &
     &              +( p1*c1*p2*exp(-p2*p2*r*r)+p3*c1*p4*exp(-p4*p4*r*r) )/r
            end if
         end do

! --- short range part of the local potential ---

         vshort(1:MMr) = vql(1:MMr,ik) - vlong(1:MMr)

         select case(pselect)
         case default

            Rcloc(ik)=Rc
            NRcloc(ik)=NRc

         case(2)

! --- filtering ---

            const = 2.d0/Pi
            fac   = 1.2d0

            NRc=MMr
            do i=1,MMr
               if ( abs(vshort(i))<ep1 ) then
                  write(*,*) i,vshort(i),rad(i,ik)
                  Rc=fac*rad(i,ik)
                  exit
               end if
            end do
            do i=1,MMr
               if ( rad(i,ik)>=Rc ) then
                  NRc=i
                  exit
               end if
            end do
            Rc=rad(NRc,ik)

            write(*,'(1x,i4,1x,i4,"  Rcloc, vshort =",2g18.10)') ik,NRc,Rc,vshort(NRc)

            Rcloc(ik) =Rc
            NRcloc(ik)=NRc

            call makemaskf

            wtmp(:)=0.d0
            vork(:)=0.d0

            c=1.d0/Rc
            maxerr=0.d0
            do i=1,NRc
               x=rad(i,ik)*c
               if ( x<=dxm ) then
                  y0=1.d0 ; dy0=0.d0
               else
                  m0=int(x/dxm)
                  dy0=1.d10
                  do m=1,20
                     m1=max(m0-m,1)       ; mm=m1-(m0-m)
                     m2=min(m0+m+mm,nmsk) ; mm=(m0+m+mm)-m2
                     call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                     if ( abs(dy)<dy0 ) then
                        y0=y ; dy0=abs(dy)
!c                        if ( dy0<ep0 ) exit
                     end if
                  end do
               end if
!c               if ( y0>1.d0 ) y0=1.d0
               wtmp(i)=y0
               maxerr=max(maxerr,dy0)
            end do
            write(*,'(1x,i4,"  maxerr(maskf,loc) =",g18.10)') ik,maxerr

            do i=1,NRc
               r=rad(i,ik)
               if ( r==0.d0 ) then
                  r1=rad(1,ik)
                  if ( r1==0.d0 ) then
                     tmp1(1)=qc*qc*qc/3.d0
                  else
                     tmp1(1)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
                  end if
                  do j=2,NRc
                     r1=rad(j,ik)
                     tmp1(j)=sin(qc*r1)/(r1*r1*r1)-qc*cos(qc*r1)/(r1*r1)
                  end do
               else
                  do j=1,NRc
                     r1=rad(j,ik)
                     if ( r1==0.d0 ) then
                        tmp1(j)=sin(qc*r)/(r*r*r)-qc*cos(qc*r)/(r*r)
                     else if ( r1==r ) then
                        tmp1(j)=(2*qc*r-sin(2.d0*qc*r))/(4*r*r*r)
                     else
                        tmp1(j)=( sin(qc*(r-r1))/(r-r1)-sin(qc*(r+r1))/(r+r1) )/(2.d0*r*r1)
                     end if
                  end do
               end if
               do j=1,NRc
                  tmp(j)=tmp1(j)*vshort(j)*rab(j,ik)*rad(j,ik)*rad(j,ik)/wtmp(j)
               end do
               call simp(tmp(1:NRc),sum0,2)
               vork(i)=sum0*const*wtmp(i)
            end do

            vshort(1:NRc) = vork(1:NRc)

         end select ! pselect

!
! --- Derivative (short) ---
!

         select case(pselect)
         case default

            maxerr=0.d0
            do i=1,MMr
               r=rad(i,ik)
               dy0=1.d10
               do m=1,20
                  m1=max(i-m,1)
                  m2=min(i+m,MMr)
                  call dpolint(rad(m1,ik),vshort(m1),m2-m1+1,r,y,dy)
                  if ( abs(dy)<abs(dy0) ) then
                     y0=y ; dy0=abs(dy)
!c                     if ( dy0<ep0 ) exit
                  end if
               end do
               dvshort(i)=y0
               maxerr=max(maxerr,dy0)
            end do
            write(*,*) "err(dvshort)=",maxerr

            maxerr=0.d0
            do i=1,MMr
               r=rad(i,ik)
               dy0=1.d10
               do m=1,20
                  m1=max(i-m,1)
                  m2=min(i+m,MMr)
                  call dpolint(rad(m1,ik),vql(m1,ik),m2-m1+1,r,y,dy)
                  if ( abs(dy)<abs(dy0) ) then
                     y0=y ; dy0=abs(dy)
                     if ( dy0<ep0 ) exit
                  end if
               end do
               dvql(i,ik)=y0
               maxerr=max(maxerr,dy0)
            end do
            write(*,*) "err(dvql)=",maxerr

            vqls(1:MMr,ik)=vshort(1:MMr)
            dvqls(1:MMr,ik)=dvshort(1:MMr)

!            vql(1:MMr,ik)=vshort(1:MMr)+vlong(1:MMr)
!            dvql(1:MMr,ik)=dvshort(1:MMr)+dvlong(1:MMr)

#ifdef TEST
            nrtmp=1000
            nqtmp=10000
            allocate( qtmp(nqtmp) )
            allocate( rtmp(nrtmp) )
            allocate( dtmp(nrtmp) )
            dq=qc*10.d0/nqtmp
            drtmp=rad(MMr,ik)/nrtmp
            write(*,*) qc*10.d0,dq,drtmp
            do j=1,nqtmp
               q=dq*(j-1)
               sum0=0.d0
               do i=1,MMr
                  x=rad(i,ik)*q
                  if ( x<=1.d-9 ) then
                     sum0=sum0+rad(i,ik)**2*vshort(i)*rab(i,ik)
                  else
                     sum0=sum0+rad(i,ik)**2*vshort(i)*rab(i,ik)*sin(x)/x
                  end if
                  qtmp(j)=sum0
               end do
            end do
            rewind 80+ik
            do j=1,nqtmp
               write(80+ik,*) dq*(j-1),qtmp(j)
            end do
            do i=1,nrtmp
               r=(i-1)*drtmp
               sum0=0.d0
               sum1=0.d0
               do j=1,nqtmp
                  q=dq*(j-1)
                  x=q*r
                  if ( x<=1.d-9 ) then
                     sum0=sum0+q*q*qtmp(j)
                     sum1=sum1+q*q*q*qtmp(j)*x/3.d0
                  else
                     sum0=sum0+q*q*qtmp(j)*sin(x)/x
                     sum1=sum1+q*q*q*qtmp(j)*(sin(x)/x**2-cos(x)/x)
                 end if
               end do
               rtmp(i)=sum0*dq*2.d0/Pi
               dtmp(i)=-sum1*dq*2.d0/Pi
            end do
            rewind 82+ik
            do i=1,nrtmp
               write(82+ik,*) drtmp*(i-1),rtmp(i),dtmp(i)
            end do
            do j=1,nqtmp
               q=dq*(j-1)
               sum0=0.d0
               do i=1,nrtmp
                  r=drtmp*(i-1)
                  x=r*q
                  if ( x<=1.d-9 ) then
                     sum0=sum0+r*r*dtmp(i)
                  else
                     sum0=sum0+r*r*dtmp(i)*sin(x)/x
                  end if
                  qtmp(j)=sum0*drtmp
               end do
            end do
            rewind 84+ik
            do j=1,nqtmp
               write(84+ik,*) dq*(j-1),qtmp(j)
            end do
            deallocate( dtmp )
            deallocate( rtmp )
            deallocate( qtmp )
#endif

         case(2)

! --- derivative (short-filtered) ---

            maxerr=0.d0
            do i=1,NRc
               r=rad(i,ik)
               dy0=1.d10
               do m=1,20
                  m1=max(i-m,1)
                  m2=min(i+m,NRc)
                  call dpolint(rad(m1,ik),vshort(m1),m2-m1+1,r,y,dy)
                  if ( abs(dy)<abs(dy0) ) then
                     y0=y ; dy0=abs(dy)
!c                     if ( dy0<ep0 ) exit
                  end if
               end do
               dvshort(i)=y0
               maxerr=max(maxerr,dy0)
            end do
            write(*,*) "err(dvshort)=",maxerr

            vqls(1:MMr,ik)=vshort(1:MMr)
            dvqls(1:MMr,ik)=dvshort(1:MMr)

            vql(1:MMr,ik)=vshort(1:MMr)+vlong(1:MMr)
            dvql(1:MMr,ik)=dvshort(1:MMr)+dvlong(1:MMr)

         end select ! pselect

      end do ! ik

!      do ik=1,MKI
!         rewind iu
!         MMr=Mr(ik) ; write(*,*) "loc,short,long,derivative",ik,iu,MMr
!         do i=1,MMr
!            write(iu,'(1x,g14.6,1x,6g14.6)') rad(i,ik),vqls(i,ik),dvqls(i,ik) &
!                                           ,vql(i,ik)-vqls(i,ik),dvql(i,ik)-dvqls(i,ik) &
!                                           ,vql(i,ik),dvql(i,ik)
!         end do
!         iu=iu+1
!      end do

!- deallcoate -----------------------------------------------------
      MMr=max( maxval(Mr),maxval(NRps) )
      deallocate( vshort,vlong,dvlong,dvshort,tmp1,tmp,vork,wtmp )
      mem=mem-bdreal*(MMr*8.d0)
!------------------------------------------------------------------

!
! --- Total # of nonlocal operators ---
!
      Mlma=0
      do a=1,MI
         ik=Kion(a)
         do iorb=1,norb(ik)
            Mlma=Mlma+2*lo(iorb,ik)+1
         end do
      end do
      write(*,*) "Mlma=",Mlma


!c      call watch(ctime1,etime1)
      write(*,*) "TIME(INIT_PS_MOL) =" !,ctime1-ctime0,etime1-etime0
      write(*,*) "MEM (MB) =",mem,memax*B2MB

      return

      CONTAINS

      SUBROUTINE simp(f,s,m)
      implicit none
      integer,intent(IN)  :: m
      real(8),intent(IN)  :: f(:)
      real(8),intent(OUT) :: s
      real(8),allocatable :: g(:)
      integer :: i,n,nmax
      n=size(f) ; nmax=int(n/m)*m
      do i=0,m
         nmax=nmax+i ; if ( nmax>=n ) exit
      end do
      allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
      select case(m)
      case default
         s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
      case(2)
         s=0.d0
         do i=1,nmax-2,2
            s = s + g(i) + 4.d0*g(i+1) + g(i+2)
         end do
         s=s/3.d0
      case(4)
         s=0.d0
         do i=1,nmax-4,4
            s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
         end do
         s=s*2.d0/45.d0
      case(6)
         s=0.d0
         do i=1,nmax-6,6
            s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3)+27*g(i+4)+216*g(i+5)+41*g(i+6)
         end do
         s=s/140.d0
      end select
      deallocate( g )
      return
      END SUBROUTINE simp


      SUBROUTINE makemaskf
      implicit none
      real(8),allocatable :: X12(:,:),iX12(:,:)
      real(8),allocatable :: V(:,:),iV(:,:)
      real(8),allocatable :: fv1(:),fv2(:),fm1(:,:)
      real(8),allocatable :: lambda(:)
      real(8) :: eta_ps,dx,x1,x2,x,c
      real(8),allocatable :: rwork(:)
      integer,parameter :: matz=1
      integer :: ierr,mxtmp,i,j,LWORK,LIWORK
      integer,allocatable :: iwork(:)

!      eta_ps=Rc*qc*qcfac
      eta_ps=etafac

      mxtmp=nmsk-1
      dx=1.d0/mxtmp
      allocate( lambda(mxtmp) )

      allocate( X12(mxtmp,mxtmp) )

      do j=1,mxtmp
         x2=j*dx
         do i=j,mxtmp
            x1=i*dx
            if ( i==j ) then
               X12(i,j)=sin((x1+x2)*eta_ps)/(x1+x2)+mxtmp*Pi-eta_ps
            else
               X12(i,j)=sin((x1+x2)*eta_ps)/(x1+x2)-sin((x1-x2)*eta_ps)/(x1-x2)
               X12(j,i)=X12(i,j)
            end if
         end do
      end do

      LWORK=1+6*mxtmp+2*mxtmp*mxtmp
      LIWORK=3+5*mxtmp
      allocate( rwork(LWORK),iwork(LIWORK) )

      call DSYEVD('V','L',mxtmp,X12,mxtmp,lambda,rwork,LWORK,iwork,LIWORK,ierr)

      deallocate( iwork,rwork )

      c=X12(1,1)/dx
      do i=1,mxtmp
         x=i*dx
         maskr(i)=X12(i,1)/x/c
      end do
      maskr(0)=1.d0
      maskr(nmsk)=0.d0
      dxm=dx
      do i=0,nmsk
         xm(i)=i*dxm
      end do

      lambda0=lambda(1)
      eta0=eta_ps

      deallocate( lambda )
      deallocate( X12 )
      return
      END SUBROUTINE makemaskf


      END SUBROUTINE init_ps_mol
