
      implicit none
      integer,parameter :: ndata=100000
      integer :: n,m,l,i,isign,idummy
      integer :: u1,u2,u3,n1,n2
      real(8),parameter :: HT=27.2116d0
      real(8) :: alpt(ndata),t(ndata),rdum1,rdum2
      real(8) :: alpt_org(ndata)
      real(8) :: dt,gamma,Ecut,dE,E,pi4,d1,A0,D0,Emax
      complex(8),parameter :: zi=(0.d0,1.d0)
      complex(8) :: zum,P,eps,ieps
      character(1) :: indx

! parameters ---------------
!
      dt=0.05d0
      A0=-0.001d0
      gamma=0.005d0
      dE=0.001d0
      Emax=3.d0

      u1=300
      u2=301
      u3=302
!
!---------------------------

      pi4=16.d0*atan(1.d0)
      n1=0
      n2=0
      alpt_org(:)=0.d0
      alpt(:)=0.d0
      t(:)=0.d0

      read(*,*) indx
      write(*,*) indx
      if ( indx=="#" ) then
         read(*,*) indx,dt,A0,n1,n2
         write(*,*) indx,dt,A0,n1,n2
      end if
      do i=1,ndata
         read(*,*,END=999) idummy,t(i),alpt_org(i)
      end do
 999  n=i-1
      
      D0=-A0
      Ecut=4.d0*atan(1.d0)/dt
      if (n1==0) n1=1
      if (n2==0) n2=n

      alpt(:)=alpt_org(:)

      rewind u3
      do i=n1,n2
         alpt(i)=alpt(i)*exp(-t(i)*gamma)
         write(u3,'(1x,i8,2x,f14.8,2x,2f20.12)')
     &        i,t(i),alpt(i),alpt_org(i)
      end do

      m=nint(Emax/dE)

      write(*,*) m,Emax,dE

      rewind u1
      rewind u2

      do l=1,m

         E=dE*l
         zum=(0.d0,0.d0)
         do i=n1,n2
            zum=zum+alpt(i)*exp(zi*E*t(i))
         end do
         zum=zum*dt

         ieps=(D0+zum)/D0
         eps =1.d0/ieps

         E=E*HT
         write(u1,*) E,real(ieps),-aimag(ieps)
         write(u2,*) E,real(eps),aimag(eps)


      end do

      stop
      end

