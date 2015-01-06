!
! Born Oppenheimer MD
! written by Yasuteru Shigeta
! Version 1.0.1
! 
!
      subroutine bomd(lmnhc)
      use global_variables
      implicit none
      integer i,ierr,itime,j,k,n1,n2
      real(8) kine,tote,tote0,pm,ltemp
      real(8) tott,bathe,rescale,dif,dthf
      logical disp_switch_loc,inivel,lmnhc
      real(8) :: ctime0,ctime1,etime0,etime1
c
      if(myrank==0) write(*,*) "start md"
c
c     setup variables
c
      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)
      disp_switch_loc=.true.
      inivel =.true.
      dt=dt*FS_TO_AU
      dthf=dt*0.5d0
      temp=temp*KB/TOJOUL
c
c     allocate and send local variables
c
      if(lmnhc) allocate(xmns(nns),vmns(nns),gmns(nns),qmns(nns))
      allocate(Velocity(3,Mi)) ; Velocity=0.0d0
      allocate(   Force(3,Mi)) ; Force   =0.0d0
      allocate(iatom(Mi))
      if(myrank==0) then
         do i=1,Mi
            do k=1,20
               if(Aion(i).eq.batm(k)) then
                  iatom(i)=k
                  exit
               endif
            enddo
         enddo
      end if
      call mpi_bcast(iatom,Mi,mpi_integer,0,mpi_comm_world,ierr)

c
c     setup initial force
c
      Econst=0.0d0
      call force_ion
      call constraint_rg
c
c     setup initial velocity
c
      if(myrank==0) then
         write(*,*) "initial velocity"
         if(inivel) then
            call setv
         else
            call mdio(0)
         endif
         if(lmnhc) then
            call setmnhc(Mi)
            call emnhc(bathe)
         else
            bathe=0.0d0
         endif
         call calkin(kine)
         tote0=kine+Etot+bathe+Econst
         if(DISP_SWITCH_LOC) then
            call mdio(2)
            write(*,'(a)') "initial energy"
            write(*,'(4f15.8)') tote0,Etot,kine,bathe
            write(*,'(a)') "initial force"
            do i=1,Mi
               write(*,'(1x,i4,i3,3g21.12)') i,Kion(i),Force(:,i)
            end do
         endif
         open(3,file='traj.dat',status='unknown')
         open(4,file='info.dat',status='unknown')
      endif
      if(lmnhc) call sicoef
c
c broadcast initial data
c
      call mpi_bcast(Velocity,Mi*3,mpi_real8,0,mpi_comm_world,ierr)
      if(lmnhc) then
         call mpi_bcast(xmns,nns,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(vmns,nns,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(gmns,nns,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(qmns,nns,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(fkt ,  1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(fnkt,  1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(kine,  1,mpi_real8,0,mpi_comm_world,ierr)
      endif
c
c write initial condition 
c
      if(myrank==0) then
         dif =0.0d0
         tott=0.0d0
         ltemp=kine/KB*TOJOUL*2.0d0/3.0d0/dble(Mi)
         write(*,*) "step=",0,"t=",tott,"Eel=",Etot
         write(4,'(7f15.8)') tott,tote0,dif,Etot,kine,bathe,ltemp
         do i=1,Mi
            write(3,'(9f12.8)') (Rion(j,i),j=1,3),
     &                          (Velocity(j,i),j=1,3),
     &                          (Force(j,i),j=1,3)
         enddo
      endif
c
c loop start
c
      do itime=1,nstep
         call watch(ctime0,etime0)
         if(lmnhc) then
            call mnhc(rescale,kine)
            Velocity=Velocity*rescale
            kine    =kine*rescale*rescale
         endif
         do i=1,Mi
            pm=dthf/(pmass(iatom(i))*AMU)
            Velocity(:,i)=Velocity(:,i)+Force(:,i)*pm
         enddo
         do i=1,Mi
            Rion(:,i)=Rion(:,i)+Velocity(:,i)*dt
         enddo
         Econst=0.0d0
         call eforce
         call constraint_rg
         do i=1,Mi
            pm=dthf/(pmass(iatom(i))*AMU)
            Velocity(:,i)=Velocity(:,i)+Force(:,i)*pm
         enddo
         call calkin(kine)
         if(lmnhc) then
            call mnhc(rescale,kine)
            Velocity=Velocity*rescale
            kine    =kine*rescale*rescale
            call emnhc(bathe)
         endif
         tote=kine+Etot+bathe
c
         if(myrank==0) then
            dif=abs(tote-tote0)
            ltemp=kine/KB*TOJOUL*2.0d0/3.0d0/dble(Mi)
            tott=dt*itime/FS_TO_AU
            write(*,*) "step=",itime,"t=",tott,"Eel=",Etot
            write(4,'(7f15.8)') tott,tote,dif,Etot,kine,bathe,ltemp
            do i=1,Mi
               write(3,'(9f12.8)') (Rion(j,i),j=1,3),
     &                             (Velocity(j,i),j=1,3),
     &                             (Force(j,i),j=1,3)
            enddo
         endif
         call watch(ctime1,etime1)
         write(14,'(i6,2f20.5)') itime,ctime1-ctime0,etime1-etime0
      enddo
c
c loop end
c
      if(DISP_SWITCH_LOC) call mdio(1)
      if(lmnhc) call unsetmnhc(bathe)
 9999 deallocate(Velocity,Force)
      return
      end subroutine bomd
c-----------------------------------------------------------------------
      subroutine calkin(kine)
      use global_variables
      implicit none
      real(8),intent(inout) :: kine
      integer i
      real(8) pm
      kine=0.0d0
      do i=1,Mi
         pm=pmass(iatom(i))
         kine=kine
     &       +(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
      enddo
      kine=kine*0.5d0*AMU
      return
      end
c-----------------------------------------------------------------------
      subroutine mdio(io)
      use global_variables
      implicit none
      integer,intent(in)    :: io
      real(8) dummy(3)
      character*50 aline
      integer i,j,k
      if(io==0) then
         open(2,file='initial.dat',status='unknown')
         read(2,'(a)') aline
            do i=1,Mi
            read(2,'(3f15.8)') (dummy(k),k=1,3)
         enddo
         read(2,'(a)') aline
         do i=1,Mi
            read(2,'(3f15.8)') (Velocity(k,i),k=1,3)
         enddo
         read(2,'(a)') aline
         do i=1,Mi
            read(2,'(3f15.8)') (dummy(k),k=1,3)
         enddo
         close(2)
      elseif(io==1) then
         open(2,file='final.dat',status='unknown')
         write(2,'(a)') "final coordinate"
         do i=1,Mi
            write(2,'(3f15.8)') (Rion(k,i),k=1,3)
         enddo
         write(2,'(a)') "final velocity"
         do i=1,Mi
            write(2,'(3f15.8)') (Velocity(k,i),k=1,3)
         enddo
         write(2,'(a)') "final force"
         do i=1,Mi
            write(2,'(3f15.8)') (Force(k,i),k=1,3)
         enddo
         close(2)
      elseif(io==2) then
         write(*,'(a)') "initial coordinate"
         do i=1,Mi
            write(*,'(3f15.8)') (Rion(k,i),k=1,3)
         enddo
         write(*,'(a)') "initial velocity"
         do i=1,Mi
            write(*,'(3f15.8)') (Velocity(k,i),k=1,3)
         enddo
         write(*,'(a)') "initial force"
         do i=1,Mi
            write(*,'(3f15.8)') (Force(k,i),k=1,3)
         enddo
      endif
      return
      end subroutine mdio
c-----------------------------------------------------------------------
      subroutine setv
      use global_variables
      implicit none
      real(8),allocatable   :: r(:)
      integer i,k
      real(8) pm,kine
c
      allocate(r(6*Mi))
      call rnum(r,Mi*6)
c
      do i=1,Mi
         Velocity(1,i)=r(6*i-5)*sin(pi2*r(6*i-4))
         Velocity(2,i)=r(6*i-3)*sin(pi2*r(6*i-2))
         Velocity(3,i)=r(6*i-1)*sin(pi2*r(6*i  ))
      enddo
      do i=1,Mi
         pm=pmass(iatom(i))*amu
         do k=1,3
            vw(k)=vw(k)+Velocity(k,i)*pm
         enddo
         vw(4)=vw(4)+pm
      enddo
      do i=1,Mi
         do k=1,3
            Velocity(k,i)=Velocity(k,i)-vw(k)/vw(4)
         enddo
      enddo
      call calkin(kine)
      Velocity=Velocity*sqrt(1.5d0*Mi*temp/kine)
      deallocate(r)
      return
      end subroutine setv
c-----------------------------------------------------------------------
      subroutine rnum(randv,N)
      use global_variables
      implicit none
      real(8),intent(inout) :: randv(N)
      integer,intent(in)    :: N
      real(8),parameter :: a=1953125.0d0
      real(8),parameter :: b=33554432.0d0
      real(8),parameter :: r0=177147.0d0
      real(8) r
      integer i
      r=r0
      do i=1,N
        r=mod(a*r,b)
        randv(i)=r/b
      enddo
      return
      end subroutine rnum
c-----------------------------------------------------------------------
      subroutine eforce
      use global_variables
      implicit none
      real(8) err,err2
c
      Miter=0
      Force=0.0d0
      call Etot_ion
      call prep_ps
      select case(MEO)
      case(1)
         call prep_hartree2
      case(2)
         call prep_hartree1
      case default
         write(*,*) "MEO is invalid! : MEO=",MEO
         call stop_program
      end select
      if(ifac_vloc>0) call calc_vloc_rho
      call scf(err,err2)
      call force_ion
c
      return
      end subroutine eforce
c-----------------------------------------------------------------------
      subroutine sicoef
      use global_variables
      implicit none
c
      real(8) :: anc,a2,a4,a8
      real(8) :: w1,w2,w3,w4,w5,c0,c1
      real(8) :: wt(7)
      integer  i
c
      anc=dfloat(nc)
      a2=dt/(anc*2.0d0)
      a4=dt/(anc*4.0d0)
      a8=dt/(anc*8.0d0)
c
      if(nys.eq.1) then
        wt2(1)=a2
        wt4(1)=a4
        wt8(1)=a8
      endif
c
      if(nys.eq.3) then
        c0   =2.0d0**(1.0d0/3.0d0)
        w1   =1.0d0/(2.0d0-c0)
        w2   =1.0d0-2.0d0*w1
        w3   =w1
        wt2(1)=w1*a2
        wt2(2)=w2*a2
        wt2(3)=w3*a2
        wt4(1)=w1*a4
        wt4(2)=w2*a4
        wt4(3)=w3*a4
        wt8(1)=w1*a8
        wt8(2)=w2*a8
        wt8(3)=w3*a8
      endif
c
      if(nys.eq.5) then
        c1   =4.0d0**(1.0d0/3.0d0)
        w1   =1.0d0/(4.0d0-c1)
        w2   =w1
        w3   =1.0d0-4.0d0*w1
        w4   =w1
        w5   =w1
        wt2(1)=w1 * a2
        wt2(2)=w2 * a2
        wt2(3)=w3 * a2
        wt2(4)=w4 * a2
        wt2(5)=w5 * a2
        wt4(1)=w1 * a4
        wt4(2)=w2 * a4
        wt4(3)=w3 * a4
        wt4(4)=w4 * a4
        wt4(5)=w5 * a4
        wt8(1)=w1 * a8
        wt8(2)=w2 * a8
        wt8(3)=w3 * a8
        wt8(4)=w4 * a8
        wt8(5)=w5 * a8
      endif
c
      if(nys.eq.7) then
          wt(1)= 0.784513610477560d0
          wt(2)= 0.235573213359357d0
          wt(3)=-1.17767998417887d0
          wt(4)= 1.0d0 - 2.0d0*(wt(1)+wt(2)+wt(3))
          wt(5)=-1.17767998417887d0
          wt(6)= 0.235573213359357d0
          wt(7)= 0.784513610477560d0
        do i=1,7
          wt2(i)=wt(i)*a2
          wt4(i)=wt(i)*a4
          wt8(i)=wt(i)*a8
        enddo
      endif
      return
      end subroutine sicoef
c-----------------------------------------------------------------------
      subroutine mnhc(sc,ekt)
      use global_variables
      implicit none
c
      real(8),intent(inout) :: sc
      real(8),intent(inout) :: ekt
      real(8) :: aa
      integer jc,iys,ins,jns,nns1
      nns1=nns+1
c
      sc=1.0d0
      gmns(1)=(ekt*2.0d0-fnkt)/qmns(1)
      do jc=1,nc
         do iys=1,nys
            vmns(nns)=vmns(nns)+gmns(nns)*wt4(iys)
            do ins=1,nns-1
               aa=exp(-wt8(iys)*vmns(nns1-ins))
               jns=nns-ins
               vmns(jns)=vmns(jns)*aa*aa+gmns(jns)*aa*wt4(iys)
            enddo
            aa=exp(-wt2(iys)*vmns(1))
            sc=sc*aa 
            gmns(1)=(sc**2*ekt*2.0d0-fnkt)/qmns(1)
            do ins=1,nns
               xmns(ins)=xmns(ins)+vmns(ins)*wt2(iys) 
            enddo
            do ins=1,nns-1
               aa=exp(-wt8(iys)*vmns(ins+1))
               vmns(ins)  = vmns(ins)*aa*aa+gmns(ins)*aa*wt4(iys)
               gmns(ins+1)=(qmns(ins)*vmns(ins)**2-fkt)/qmns(ins+1)
            enddo
            vmns(nns)     =vmns(nns)+gmns(nns)*wt4(iys)
         enddo
      enddo
      return
      end subroutine mnhc
c-----------------------------------------------------------------------
      subroutine setmnhc(np)
      use global_variables
      implicit none
c
      integer,intent(IN) :: np
      real(8),parameter :: wave =2000.0d0       !(1/cm)
      real(8) :: randv(nns*6)
      real(8) :: tunit,fkt2m
      integer l,ins
c
      call rnum(randv,nns*6)
c
      fkt =temp
      fnkt=fkt * np *3.0d0
      fkt2m=-2.0d0 * fkt
      tunit=tojoul / (clight * planck * wave)
c
      write(*,1000) fkt,fnkt,fkt * tunit * tunit
      write(*,1001) tunit,1.0d0/tunit
 1000 format("fkt =",f10.6,2x,"fnkt=",f10.6,2x,"qmns=",f10.6)
 1001 format("tunit=",f10.6,2x,"De =",f10.6)
c
      l=0
      do ins=1,nns
         l=l+1
         qmns(ins)=fkt*tunit*tunit
         xmns(ins)=1.0d0
         vmns(ins)=sqrt(fkt2m*dlog(randv(2*l-1))/qmns(ins))
     *             * sin(pi2*randv(2*l))
         gmns(ins)=0.0d0
      enddo
      qmns(1)=fnkt*tunit*tunit
      do ins=1,nns-1
         gmns(ins+1)=(qmns(ins)*vmns(ins)**2-fkt)/qmns(ins+1)
      enddo
      return
      end subroutine setmnhc
c-----------------------------------------------------------------------
      subroutine unsetmnhc(bathe)
      use global_variables
      implicit none
c
      real(8),intent(in) :: bathe
      integer :: k
c
      write(*,'(2f15.8)') bathe
      write(*,'(8f15.8)') (xmns(k),k=1,nns)
      write(*,'(8f15.8)') (vmns(k),k=1,nns)
      write(*,'(8f15.8)') (gmns(k),k=1,nns)
      write(*,'(8f15.8)') (qmns(k),k=1,nns)
      deallocate(xmns,vmns,gmns,qmns)
c
      return
      end subroutine unsetmnhc
c-----------------------------------------------------------------------
      subroutine emnhc(bathe)
      use global_variables
      implicit none
      real(8),intent(inout) :: bathe
c
      real(8) :: bathk,bathp
      integer ins
c
      bathk=0.0d0
      do ins=1,nns
         bathk=bathk+vmns(ins)**2*qmns(ins)
      enddo
c
      bathp=xmns(1)*fnkt
      do ins=2,nns
         bathp=bathp+xmns(ins)*fkt
      enddo
c
      bathk=bathk * 0.5d0
      bathe=bathk+bathp
c
      return
      end subroutine emnhc
c-----------------------------------------------------------------------
      subroutine constraint_rg_old
      use global_variables
      implicit none
      real(8) disrg,krg,rg0,pm,gcent(3)
      integer i,j,k
      krg=1.0d3
      rg0=0.0d0
      gcent=0.0d0
      vw=0.0d0
      do i=1,Mi
         pm=pmass(iatom(i))*amu
         do k=1,3
            vw(k)=vw(k)+Rion(k,i)*pm
         enddo
         vw(4)=vw(4)+pm
      enddo
      do k=1,3
         vw(k)=vw(k)/vw(4)
      enddo
      if(myrank==0) write(*,*) vw
      disrg=(vw(1)-gcent(1))**2+(vw(2)-gcent(2))**2+(vw(3)-gcent(3))**2
      Econst=Econst+krg*disrg*0.5d0
      disrg=sqrt(disrg)
      if(disrg.ne.rg0) then
         do i=1,Mi
            pm=pmass(iatom(i))*amu/vw(4)
            do k=1,3
               Force(k,i)=Force(k,i)
     &                   -krg*(1.0d0-rg0/disrg)*(pm*vw(k)-gcent(k))
            enddo
         enddo
      endif
      return
      end subroutine constraint_rg_old
c-----------------------------------------------------------------------
      subroutine constraint_rg
      use global_variables
      implicit none
      real(8) disrg,krg,pm
      integer i,j,k
      krg=1.0d2
      vw=0.0d0
      do i=1,Mi
         pm=pmass(iatom(i))*amu
         do k=1,3
            vw(k)=vw(k)+Rion(k,i)*pm
         enddo
         vw(4)=vw(4)+pm
      enddo
      do k=1,3
         vw(k)=vw(k)/vw(4)
      enddo
      disrg=vw(1)**2+vw(2)**2+vw(3)**2
      Econst=Econst+krg*disrg*0.5d0
      disrg=sqrt(disrg)
      if(disrg.ne.0.0d0) then
         do i=1,Mi
            pm=pmass(iatom(i))*amu/vw(4)
            do k=1,3
               Force(k,i)=Force(k,i)-krg*pm*vw(k)
            enddo
         enddo
      endif
      return
      end subroutine constraint_rg
