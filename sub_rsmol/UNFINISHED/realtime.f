!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE realtime
      use global_variables
      implicit none

      integer,parameter :: u1=21,u2=22,u3=23
      integer :: Nhamil,Ntime,dir
      real(8) :: max_etime
!_shigeta      real(8) :: dt,F,dE,t,t1
      real(8) :: F,dE,t,t1
      real(8) :: Hconv
      real(8),allocatable :: polt(:)
      complex(8),parameter :: zi=(0.d0,1.d0),z0=(0.d0,0.d0)
      complex(8),allocatable :: alpha(:)
      complex(8),allocatable :: zc(:),zpsi(:,:),htpsi(:),tpsi(:)
      complex(8) :: zsum
      integer :: Nenergy,MHiter,ix,iy,iz
      integer :: m,n,it,p,q,Ntime0,i,n1,n2,ML0,ierr
      real(8) :: Ds(3,3),Nele,Etot0,hw,Tmax,Emax,s(2),r(2),const
      real(8),allocatable :: rho0(:),Vloc0(:),Vext(:),esp0(:)
      character(6) :: label
      logical :: DISP_SWITCH_LOC

      if (DISP_SWITCH) then
         write(*,'(a60," Real-Time")') repeat("-",60)
      end if

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)

      DISP_SWITCH_LOC = DISP_SWITCH
      DISP_SWITCH     = .false.

*
* --- input ---
*
      if ( myrank==0 ) then
         do i=1,1000
            read(unit,'(a6)') label
            if ( label=='# RSRT' ) exit
            if ( i==1000 ) then
               write(*,*) "Label '# RSRT' is not found."
               goto 900
            end if
         end do
         read(unit,*) dt,Ntime
         read(unit,*) F
         read(unit,*) dir
         read(unit,*) Nhamil
c         read(unit,*) Hconv,MHiter
      end if

      call mpi_bcast(dt    ,1,mpi_real8  ,0,mpi_comm_world,ierr)
      call mpi_bcast(Ntime ,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(F     ,1,mpi_real8  ,0,mpi_comm_world,ierr)
      call mpi_bcast(dir   ,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(Nhamil,1,mpi_integer,0,mpi_comm_world,ierr)

      Tmax = dt*Ntime
      Emax=1.47d0
      dE=0.0003674d0
      Nenergy=Emax/dE
      const=1.d0/F*dt

      if (DISP_SWITCH_LOC) then
         write(*,*)
         write(*,*) "Total time step      =",Ntime
         write(*,*) "Time step            =",dt
         write(*,*) "Field strength       =",F
         write(*,*) "dir of field         =",dir
         write(*,*) "Order of Taylor Exp  =",Nhamil
c         write(*,*) "Hartree(Hconv,MHiter)=",Hconv,MHiter
         write(*,*) "Tmax=",Tmax
         write(*,*) "Emax=",Emax
         write(*,*) "dE=",dE
         write(*,*) "Nenergy=",Nenergy
      end if

      if ( dt<0.d0 ) goto 900
      if ( dir<1 .or. 3<dir ) goto 900

!- allocate --------------------------------------------------------

      allocate( alpha(0:Nenergy) ) ; alpha=z0
      allocate( polt(0:Ntime) )
      allocate( zc(Nhamil) )
      allocate( Vext(n1:n2) )
      allocate( Vloc0(n1:n2),rho0(n1:n2),esp0(MST) )
      allocate( zpsi(n1:n2,MST),tpsi(n1:n2),htpsi(n1:n2) )

      allocate(
     &     www_c(Mx0t-Md:Mx1t+Md,My0t-Md:My1t+Md,Mz0t-Md:Mz1t+Md) )
      www_c=(0.d0,0.d0)

* for eqdiv_mesh
c      m=N_data_max
c      n=N_neighbor_max
* for normal mesh
      ix=Md*(My1t-My0t+1)*(Mz1t-Mz0t+1)
      iy=Md*(Mz1t-Mz0t+1)*(Mx1t-Mx0t+1)
      iz=Md*(Mx1t-Mx0t+1)*(My1t-My0t+1)
      m=max(ix,iy,iz)
      n=maxval(NMBC)

      allocate( sbuf_c(m,n,6),rbuf_c(m,n,6) )
      sbuf_c=(0.d0,0.d0)
      rbuf_c=(0.d0,0.d0)

      m=0
      do n=0,nprocs-1
         if ( n/=myrank ) m=max(m,lma_nsend(n))
      end do
      allocate( sbufnl_c(m,0:nprocs-1),rbufnl_c(m,0:nprocs-1) )
      sbufnl_c=(0.d0,0.d0)
      rbufnl_c=(0.d0,0.d0)

!-------------------------------------------------------------------

*
* --- external field ---
*
      Vext(n1:n2)=F*LL(dir,n1:n2)*H

*
* --- Static dipole moment ---
*
      Ds(1,1)=sum(LL(1,n1:n2)*rho(n1:n2))*dV*H
      Ds(2,1)=sum(LL(2,n1:n2)*rho(n1:n2))*dV*H
      Ds(3,1)=sum(LL(3,n1:n2)*rho(n1:n2))*dV*H

      call mpi_allreduce
     &     (Ds(1,1),Ds(1,2),3,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      Ds(:,3)=0.d0
      do i=1,MI
         Ds(1,3)=Ds(1,3)+Zps(Kion(i))*Rion(1,i)
         Ds(2,3)=Ds(2,3)+Zps(Kion(i))*Rion(2,i)
         Ds(3,3)=Ds(3,3)+Zps(Kion(i))*Rion(3,i)
      end do

      Ds(1:3,1)=Ds(1:3,3)-Ds(1:3,2)

      if (DISP_SWITCH_LOC) then
         write(*,'(1x,"Static dipole moment =",3f15.8)') Ds(1,1:3)
         write(*,'(1x,"                     =",3f15.8)') Ds(2,1:3)
         write(*,'(1x,"                     =",3f15.8)') Ds(3,1:3)
      end if

*
* --- Taylor expansion coefficient ---
*
      do n=1,Nhamil
         zc(n)=(-zi*dt)**n
         do m=1,n
            zc(n)=zc(n)/m
         end do
      end do

*
* --- Initial wf ---
*
      do p=1,MST
         zpsi(n1:n2,p)=exp(zi*Vext)*psi(n1:n2,p)
      end do

      polt(0)=0.d0

      rho0(n1:n2)=rho(n1:n2)
      Vloc0(n1:n2)=Vloc(n1:n2)
      esp0(1:MST)=esp(1:MST)
      Etot0=Etot
      Nele=sum(rho)*dV

      rewind u1
      rewind u2

      if (DISP_SWITCH_LOC) then
         write(*,'(1x,i6,2x,f12.8,3f18.12)') 0,0.d0,polt(0),Nele,Etot
      end if

      if ( myrank==0 ) then
         write(u1,*) 0,0.d0,polt(0)
         write(u2,*) 0.d0,Nele,Etot
      end if

      Ntime0=1

*
* --- Time evolution ---
*
      Time_evolution : do it=Ntime0,Ntime

         t=it*dt

         do p=1,MST
            tpsi=zpsi(:,p)
            do n=1,Nhamil
               call hpsi_c(tpsi,htpsi)
               zpsi(:,p)=zpsi(:,p)+zc(n)*htpsi
               tpsi=htpsi
            end do
         end do

         rho(n1:n2)=occ(1)*abs(zpsi(n1:n2,1))**2
         do p=2,MST
            rho(n1:n2)=rho(n1:n2)+occ(p)*abs(zpsi(n1:n2,p))**2
         end do

         call Hartree(rho,Vh,1.d-25,2000)
         call Exc_Cor
         Vloc(n1:n2)=Vion(n1:n2)+Vh(n1:n2)+Vxc(n1:n2) !+Vext(n1:n2)

*
* result
*
         s(1) = sum(rho(n1:n2))*dV
         s(2) = sum(LL(dir,n1:n2)*rho(n1:n2))*dV*H
         call mpi_allreduce(s,r,2,mpi_real8
     &        ,mpi_sum,mpi_comm_world,ierr)
         Nele     = r(1)
         polt(it) = r(2) - Ds(dir,1)

         call etot_rt
*
* --- polarizability ---
*
         do m=0,Nenergy
            hw=m*dE
            zsum=exp(zi*hw*t)*polt(it)*(1-3*(t/Tmax)**2+2*(t/Tmax)**3)
            alpha(m)=alpha(m)+zsum/F*dt
         end do


         if ( DISP_SWITCH_LOC ) then
            write(*,'(1x,i5,2x,f7.3,5f15.8)') it,t,polt(it),Nele,Etot
     &           ,real(alpha(0))
         end if

         if ( myrank==0 ) then
            write(u1,*) it,t,polt(it)
            write(u2,*) t,Nele,Etot
            rewind u3
            do n=0,Nenergy
               hw=n*dE
               write(u3,*) hw,real(alpha(n)),aimag(alpha(n))
            end do
         end if

      end do Time_evolution

!----------------------------------------------------- Fourier transform

c      Tmax = dt*Ntime
c      Nenergy = Ntime
c      Emax = Pi/dt
c      dE = Emax/Nenergy
c      alpha=z0
c      do m=0,Nenergy
c         hw=m*dE
c         zsum=z0
c         do n=1,Ntime
c            t=n*dt
c            zsum=zsum+exp(zi*hw*t)*polt(n)
c     &           *(1-3*(t/Tmax)**2+2*(t/Tmax)**3)
c         end do
c         alpha(m)=zsum/F*dt
c      end do
c      rewind u3
c      do n=0,Nenergy
c         hw=n*dE
c         write(u3,*) hw,real(alpha(n)),aimag(alpha(n))
c      end do

      deallocate( rbufnl_c,sbufnl_c )
      deallocate( rbuf_c,sbuf_c )
      deallocate( www_c )

      deallocate( htpsi,tpsi,zpsi )
      deallocate( esp0,rho0,Vloc0 )
      deallocate( Vext )
      deallocate( zc )
      deallocate( polt )
      deallocate( alpha )

      return

 900  call stop_program1("real_time")

      CONTAINS

         subroutine etot_rt
         implicit none
         integer :: p,ierr
         real(8),allocatable :: esp_tmp(:)
         real(8) :: s(2),r(2)

         allocate( esp_tmp(MST) ) ; esp_tmp=0.d0

         do p=1,MST
            if ( occ(p)<1.d-10 ) cycle
            call hpsi_c(zpsi(:,p),htpsi)
            esp_tmp(p)=sum( conjg(zpsi(:,p))*htpsi(:) )*dV
         end do

         call mpi_allreduce(esp_tmp,esp,MST,mpi_real8
     &        ,mpi_sum,mpi_comm_world,ierr)

         s(1)=sum( Vh(n1:n2)*rho(n1:n2) )*dV
         s(2)=sum( Vxc(n1:n2)*rho(n1:n2) )*dV
         call mpi_allreduce
     &        (s,r,2,mpi_real8,mpi_sum,mpi_comm_world,ierr)

         Etot=sum(occ*esp)-0.5d0*r(1)-r(2)+Exc+Eion

         deallocate( esp_tmp )

         end subroutine etot_rt

      END SUBROUTINE realtime
