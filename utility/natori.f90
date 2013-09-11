!--------1---------2---------3---------1---------5---------6---------7--
      implicit none
      real(8),parameter :: HT=27.2116d0, aB=0.529177d0
      real(8),allocatable :: esp(:,:),kbb(:),dedk(:,:),d2edk2(:,:)
      real(8) :: mu_s,mu_d,Vd,kT,eVBM,eCBM,Vgt
      real(8) :: beta,Q,Qf,Qb,fs,fd,ax,epsilon,pi,Cg,Tox,Radius
      real(8) :: dVd,alpha,Q1,mu_s_r,mu_s_l,dk,G0
      real(8) :: cfd0(0:6),bfd0(0:6),cfd(0:6),bfd(0:6)
      real(8) :: Isdf,Isdb,eminf,emaxf,eminb,emaxb
      real(8),allocatable :: tmp(:,:,:),dtmp(:)
      integer :: k,n,idummy,i0,i1,i2,NVd,i,j,m,iloc(1)
      integer :: Md,MB,MBZ,MB1,MB2,MBv,MBc,iflag
      real(8),parameter :: qe=1.602176462d-19 !(C)
      real(8),parameter :: J0=243.413389d0    !(uA)

      cfd0(:)=0.d0
      bfd0(:)=0.d0
!      Md=1 ; bfd0(1)=0.5d0
!      cfd0(0)=-2.d0 ; cfd0(1)=1.d0
      Md=2 ; bfd0(1)=2.d0/3.d0 ; bfd0(2)=-1.d0/12.d0
      cfd0(0)=-2.5d0 ; cfd0(1)=4.d0/3.d0 ; cfd0(2)=-1.d0/12.d0

      read(*,*) MB,MBZ,MB1,MB2,MBv,MBc
      write(*,*) MB,MBZ,MB1,MB2,MBv,MBc

      MBZ=MBZ-1

      write(*,*) "MBZ=",MBZ

      allocate( esp(MB,-MBZ-Md:MBZ+Md)   ) ; esp=0.d0
      allocate( kbb(-MBZ-Md:MBZ+Md)      ) ; kbb=0.d0
      allocate( dedk(MB,-MBZ:MBZ)        ) ; dedk=0.d0
      allocate( d2edk2(MB,-MBZ:MBZ)      ) ; d2edk2=0.d0

      do k=0,MBZ
         write(*,*) k
         do n=1,MB
            read(*,*) kbb(k),idummy,esp(n,k)
         end do
      end do

      do k=1,MBZ
         do n=1,MB
            esp(n,-k)=esp(n,k)
         end do
         kbb(-k)=-kbb(k)
      end do
      do k=MBZ+1,MBZ+Md
         do n=1,MB
            esp(n,k)=esp(n,2*MBZ-k)
         end do
         kbb(k)=kbb(MBZ)+( kbb(MBZ)-kbb(MBZ-1) )
      end do
      do k=-MBZ-1,-MBZ-Md,-1
         do n=1,MB
            esp(n,k)=esp(n,-2*MBZ-k)
         end do
         kbb(k)=kbb(-MBZ)+( kbb(-MBZ)-kbb(-MBZ+1) )
      end do

      if ( 0<MBv .and. MBv<MB ) then
         eVBM=maxval( esp(1:MBv,0:MBZ-1) )
         eCBM=minval( esp(MBc:MB,0:MBZ-1) )
      else
         eVBM=0.d0
         eCBM=0.d0
      end if

      write(*,*) eVBM,eCBM

!      esp(:,:)=esp(:,:)-eCBM
      esp(:,:)=esp(:,:)-eVBM

      rewind 11
      do k=-MBZ-Md,MBZ+Md
         write(11,'(1x,i4,2x,50f12.7)') k,esp(1:MB,k)
      end do

! derivative

      allocate( tmp(MB,-MBZ-Md:MBZ+Md,3) ) ; tmp=0.d0
      tmp(:,:,1)=esp(:,:)

      i0=20
      allocate( dtmp(-i0:i0) ) ; dtmp=0.d0

      rewind 10
      do k=-MBZ,MBZ
         dk=kbb(k+1)-kbb(k)
         bfd(1:Md)=bfd0(1:Md)/dk
         cfd(0:Md)=cfd0(0:Md)/dk**2
         do n=1,1
!            i1=max(n-i0, 1)-n
!            i2=min(n+i0,MB)-n
!            write(*,*) "k=",k,i1,i2
!            do i=i1,i2
!               dtmp(i)=cfd(0)*tmp(n,k,1)
!               do j=1,Md
!                  dtmp(i)=dtmp(i)
!     &                 +cfd(j)*(tmp(n+i,k+j,1)+tmp(n,k-j,1))
!               end do
!               write(*,*) i,n+i,dtmp(i)
!            end do
!            if(k==0)then
!               tmp(n,k,3)=dtmp(0)
!            else
!               iloc(:)=minloc(abs(dtmp(i1:i2)))
!               m=iloc(1)-1+i1
!               tmp(n,k+1,1)=tmp(n+m,k+1,1)
!               tmp(n,k,3)=dtmp(m)
!            end if
            tmp(n,k,2)=0.d0
            do j=1,Md
               tmp(n,k,2)=tmp(n,k,2)+bfd(j)*( tmp(n,k+j,1)-tmp(n,k-j,1) )
            end do
            write(10,'(1x,i4,2x,3g16.8)') k,k*dk,tmp(n,k,1:2)
         end do
      end do
      stop


      Tox     = 10.d0/aB
      Radius  = 80.d0/aB
      pi      = 4.d0*atan(1.d0)
      epsilon = 4.0d0
      Cg      = 2.d0*pi*epsilon/log((Radius+Tox)/Radius)

      write(*,*) "Cg=",Cg,pi

      alpha = 1.d0

      ax   = 10.261d0
      dk   = 2.d0*Pi/ax/(2*MBZ+1)
!      dk   = 2.d0*Pi/(2*MBZ+1)
      kT   = 0.026/HT
      beta = 1.d0/kT
!      mu_s = 1.d0/Ht
!      mu_d = 0.5d0/Ht
!      Vd   = mu_s-mu_d
      Vgt  = 0.5d0/Ht

      NVd = 100
      dVd = 1.5d0/Ht/NVd

      rewind 11

      do i=1,NVd

         Vd = (i-1)*dVd

         mu_s=0.d0

         iflag=0

         do j=1,12

            mu_d=mu_s-Vd

            write(*,'(1x,"mu_s,mu_d,Vd=",3f15.8)') mu_s,mu_d,Vd

            Qf=0.d0
            Qb=0.d0
            do n=MBc,MBc
               do k=-MBZ,MBZ
                  fs=2.d0/( 1.d0+exp( beta*(esp(n,k)-mu_s)) )
                  fd=2.d0/( 1.d0+exp( beta*(esp(n,k)-mu_d)) )
                  if(dedk(n,k)>0.d0)then
                     Qf=Qf+fs
                  else if(dedk(n,k)<0.d0)then
                     Qb=Qb+fd
                  end if
               end do
            end do
            Qf=Qf*dk/(2.d0*Pi)
            Qb=Qb*dk/(2.d0*Pi)
            Q=Qf+Qb

            Q1=Cg*Vgt-Cg*alpha*mu_s

!            write(*,'(1x,"Qf,Qb,Q,Q1",4f15.8)') Qf,Qb,Q,Q1

            if ( Q<Q1 ) then
               if(iflag==2)exit
               iflag=1
               mu_s_l=mu_s
               mu_s_r=mu_s+0.01d0
               mu_s=mu_s_r
            end if
            if ( Q1<Q ) then
               if(iflag==1)exit
               iflag=2
               mu_s_l=mu_s-0.01d0
               mu_s_r=mu_s
               mu_s=mu_s_l
            end if
            if ( Q1==Q ) stop

         end do ! j

         write(*,*)
         write(*,*) "mu_s_l,mu_s_r=",mu_s_l,mu_s_r
         write(*,*)
         if(j>12)stop

         do j=1,120

            mu_s = 0.5d0*( mu_s_l + mu_s_r )

            mu_d=mu_s-Vd

!            write(*,'(1x,"mu_s,mu_d,Vd=",3f18.12)') mu_s,mu_d,Vd

            Qf=0.d0
            Qb=0.d0
            do n=MBc,MBc
               do k=-MBZ,MBZ
                  fs=2.d0/( 1.d0+exp( beta*(esp(n,k)-mu_s)) )
                  fd=2.d0/( 1.d0+exp( beta*(esp(n,k)-mu_d)) )
                  if(dedk(n,k)>0.d0)then
                     Qf=Qf+fs
                  else if(dedk(n,k)<0.d0)then
                     Qb=Qb+fd
                  end if
               end do
            end do
            Qf=Qf*dk/(2.d0*Pi)
            Qb=Qb*dk/(2.d0*Pi)
            Q=Qf+Qb

            Q1=Cg*Vgt-Cg*alpha*mu_s

!            write(*,'(1x,i4,1x,"Qf,Qb,Q,Q1=",4f15.8)') j,Qf,Qb,Q,Q1

            if ( abs(Q-Q1)<1.d-12 ) exit

            if ( Q<Q1 ) then
               mu_s_l=mu_s
            else if ( Q1<Q ) then
               mu_s_r=mu_s
            else
               stop
            end if

         end do ! j
         write(*,'(1x,"mu_s,mu_d,Vd=",3f18.12)') mu_s,mu_d,Vd
         write(*,'(1x,i4,1x,"Qf,Qb,Q,Q1=",4f15.8)') j,Qf,Qb,Q,Q1
         if(j>120)stop
!
!
!
         Isdf=0.d0
         Isdb=0.d0
!         mu_s=mu_s+Vd
!         mu_d=mu_d+Vd
         do n=MBc,MBc
!            do k=-MBZ,MBZ
!               if ( dedk(n,k)>0.d0 ) then
!                  Isdf=Isdf
!     &            +2.d0/(1.d0+exp(beta*(esp(n,k)-mu_s)))*dedk(n,k)
!               else if ( dedk(n,k)<0.d0 ) then
!                  Isdb=Isdb
!     &            +2.d0/(1.d0+exp(beta*(esp(n,k)-mu_d)))*dedk(n,k)
!               end if
!            end do
!            Isdf=Isdf*dk/(2.d0*Pi)
!            Isdb=Isdb*dk/(2.d0*Pi)
            eminf= 10.d10
            emaxf=-10.d10
            eminb= 10.d10
            emaxb=-10.d10
            do k=-MBZ,MBZ
               if ( dedk(n,k)>0.d0 ) then
                  eminf=min(eminf,esp(n,k))
                  emaxf=max(emaxf,esp(n,k))
               else if ( dedk(n,k)<0.d0 ) then
                  eminb=min(eminb,esp(n,k))
                  emaxb=max(emaxb,esp(n,k))
               end if
            end do
            Isdf=Isdf-kT*log( (1.d0+exp(beta*(mu_s-emaxf)))/(1.d0+exp(beta*(mu_s-eminf))) )/Pi
            Isdb=Isdb+kT*log( (1.d0+exp(beta*(mu_d-emaxb)))/(1.d0+exp(beta*(mu_d-eminb))) )/Pi
         end do

!         write(*,*) eminf,emaxf
!         write(*,*) eminb,emaxb
         write(*,*) Isdf,Isdb

         write(11,*) Vd*HT,J0*HT*(Isdf+Isdb)

      end do ! i

      stop
      end
