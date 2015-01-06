c-----------------------------------------------------------------------
c Driver rootine for Verlet CPMD with two thermostadt
c-----------------------------------------------------------------------
      subroutine cpmd(l2therm)
      use global_variables
      implicit none
      integer i,ierr,itime,j,k,n1,n2
      real(8) kine,tote,tote0,pm,ltemp,tott,bathe,dif,dt2,dthf,dtinvhf
      real(8) fke
      real(8), allocatable :: Rnew(:,:),Rold(:,:)
      real(8) ae,ar,co,cc,cf,tunit
      real(8) barmass,cormass,rmass,smass,effmass,aeff
      logical disp_switch_loc,inivel,l2therm
      real(8) :: ctime0,ctime1,etime0,etime1
c
      if(myrank==0) write(*,*) "start cpmd"
c
c     setup variables
c
      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)
      disp_switch_loc=.true.
      inivel =.true.
c
      dt =dt*FS_TO_AU
      dt2=dt*dt
      dthf=dt*0.5d0
      dtinvhf=0.5d0/dt
      cf=dt2/emass*2.0d0
      temp=temp*KB/TOJOUL
      fnkt=temp*3.0d0*Mi
      MSTOCC =aint(sum(occ))/2
      MSTOCC2=(MSTOCC+1)*MSTOCC/2
c
c     allocate and send local variables
c
      if(l2therm) allocate(xe(-3:1),xr(-3:1))
      allocate(Rnew(3,Mi))
      allocate(Rold(3,Mi))
      allocate(Velocity(3,Mi))
      allocate(   Force(3,Mi))
      allocate(psi_c(n1:n2,MSTOCC)); psi_c=psi(n1:n2,1:MSTOCC)
      allocate(psi_v(n1:n2,MSTOCC)); psi_v=0.0d0
      allocate(psi_n(n1:n2,MSTOCC))
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
      call Total_Energy_cpmd
      call force_ion
      call constraint_rg
      call calfke(fke)
      psi_v=psi_c
      do i=1,MSTOCC
         call hpsi(psi_c(n1,i),psi_n(n1,i))
      end do
c
c     setup initial velocity
c
      if(myrank==0) then
         if(inivel) then
            call setv
         else
            call mdio(0)
         endif
         call calkin(kine)
         if(l2therm) then
            call settherm2(bathe,kine,fke,dt2)
         else
            bathe=0.0d0
         endif
         tote0=kine+Etot+bathe+fke+Econst
         tote =tote0
         if(DISP_SWITCH_LOC) then
            write(*,'(a)') "initial energy"
            write(*,'(4f15.8)') tote0,Etot,kine+fke,bathe
            write(*,'(a)') "initial force"
            do i=1,Mi
               write(*,'(1x,i4,i3,3g21.12)') i,Kion(i),Force(:,i)
            end do
         endif
         open(3,file='traj.dat',status='unknown')
         open(4,file='info.dat',status='unknown')
      endif
c
c broadcast initial data
c
      call mpi_bcast(Velocity,Mi*3,mpi_real8,0,mpi_comm_world,ierr)
      if(l2therm) then
         call mpi_bcast(xe,5,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(xr,5,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(ve,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(vr,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(qe,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(qr,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(fe,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(fr,1,mpi_real8,0,mpi_comm_world,ierr)
         call mpi_bcast(bathe,1,mpi_real8,0,mpi_comm_world,ierr)
      endif
c
c write initial condition 
c
      if(myrank==0) then
         dif =0.0d0
         tott=0.0d0
         ltemp=kine/KB*TOJOUL*2.0d0/3.0d0/dble(Mi)
         write(4,111) tott,tote,dif,Etot,kine,fke,bathe,ltemp
         write(*,111) tott,tote,dif,Etot,kine,fke,bathe,ltemp
         do i=1,Mi
            write(3,'(9f12.8)') (Rion(j,i),j=1,3),
     &                          (Velocity(j,i),j=1,3),
     &                          (Force(j,i),j=1,3)
         enddo
      endif
c
c set back to previous time step
c
      ae=ve*dthf
      ar=vr*dthf
      do i=1,Mi
         barmass=pmass(iatom(i))*AMU
         cormass=kmass(iatom(i))*emass*4.0d0/3.0d0
         effmass=barmass-cormass
         rmass=barmass/effmass
         smass=cormass/effmass
         aeff=ar*rmass-ae*smass
         aeff=(1.0d0+aeff)*dt
         pm=dt2/effmass*0.5d0
         Rold(:,i)=Rion(:,i)-aeff*Velocity(:,i)+Force(:,i)*pm
      enddo

      if(l2therm) bathe=bathe-0.5d0*(ve*ve*qe+vr*vr*qr)
c
c loop start
c
      do itime=0,nstep
         call watch(ctime0,etime0)
         if(l2therm) then
            ve=(4.0d0*(xe(0)+xe(-2))-7.0d0*xe(-1)-xe(-3))*dtinvhf
            vr=(4.0d0*(xr(0)+xr(-2))-7.0d0*xr(-1)-xr(-3))*dtinvhf
            ae=ve*dthf
            ar=vr*dthf
            do i=1,Mi
               barmass=pmass(iatom(i))*AMU
               cormass=kmass(iatom(i))*emass*4.0d0/3.0d0
               effmass=barmass-cormass
               rmass=barmass/effmass
               smass=cormass/effmass
               aeff=ar*rmass-ae*smass
               pm=dt2/effmass
               cc=2.0d0/(1.0d0+aeff)
               co=(1.0d0-aeff)/(1.0d0+aeff)
               cf=pm/(1.0d0+aeff)
               Rnew(:,i)=cc*Rion(:,i)-co*Rold(:,i)+cf*Force(:,i)
            enddo
            cc=2.0d0/(1.0d0+ae)
            co=(1.0d0-ae)/(1.0d0+ae)
            cf=dt2/emass*2.0d0/(1.0d0+ae)
            psi_n=cc*psi_c-co*psi_v-cf*psi_n
         else
            do i=1,Mi
               pm=dt2/(pmass(iatom(i))*AMU)
               Rnew(:,i)=2.0d0*Rion(:,i)-Rold(:,i)+pm*Force(:,i)
            enddo
            psi_n=2.0d0*psi_c-psi_v-cf*psi_n
         endif
         call rotorb
c
         Velocity=(Rnew-Rold)*dtinvhf
         call calkin(kine)
         Rold=Rion
         Rion=Rnew
         psi_v=(psi_n-psi_v)*dtinvhf
         call calfke(fke)
         psi_v=psi_c
         psi_c=psi_n
         if(l2therm) then
            if(ve.lt.0.0d0) then
               xe(1)=xe(0)
            else
               xe(1)=2.0d0*xe(0)-xe(-1)+fe*dt2
            endif
            xr(1)=2.0d0*xr(0)-xr(-1)+fr*dt2
            ve=(xe(1)-xe(-1))*dtinvhf
            vr=(xr(1)-xr(-1))*dtinvhf
            do i=-3,0
               xe(i)=xe(i+1)
               xr(i)=xr(i+1)
            enddo
            bathe=bathe+0.5d0*(ve*ve*qe+vr*vr*qr)
         endif
         tote=kine+Etot+bathe+fke+Econst
c
         if(myrank==0) then
            dif=abs(tote-tote0)
            ltemp=kine/KB*TOJOUL*2.0d0/3.0d0/dble(Mi)
            tott=dt*itime/FS_TO_AU
            write(4,111) tott,tote,dif,Etot,kine,fke,bathe,ltemp
            write(*,111) tott,tote,dif,Etot,kine,fke,bathe,ltemp
            do i=1,Mi
               write(3,'(9f12.8)') (Rion(j,i),j=1,3),
     &                             (Velocity(j,i),j=1,3),
     &                             (Force(j,i),j=1,3)
            enddo
         endif
!Force
         Econst=0.0d0
         call eforce_cpmd(.true.)
         call constraint_rg
         do i=1,MSTOCC
            call hpsi(psi_c(n1,i),psi_n(n1,i))
         end do
         if(l2therm) then
            call calkin2(etemp)
            fe=2.0d0*(fke-etemp)/qe
            fr=(2.0d0*kine-fnkt)/qr
            if(ve.lt.0.0d0) then
               fe=0.0d0
            elseif(ve.eq.0.0d0) then
               fe=fe*0.5d0
            endif
            bathe=xr(0)*fnkt-etemp
         endif
         call watch(ctime1,etime1)
         write(14,'(i6,2f20.5)') itime,ctime1-ctime0,etime1-etime0
!Force
      enddo
c
c loop end
c
      if(DISP_SWITCH_LOC) call mdio(1)
 9999 deallocate(Velocity,Force)
      deallocate(psi_c,psi_n,psi_v)
      deallocate(Rnew,Rold,iatom)
      if(l2therm) deallocate(xe,xr)
 111  format(8(f15.8,","))
      return
      end subroutine cpmd
c-----------------------------------------------------------------------
c Energy and Force for CPMD
c-----------------------------------------------------------------------
      subroutine eforce_cpmd(lcpmd)
      use global_variables
      implicit none
      integer n1,n2
      real(8) err,err2
      logical lcpmd
      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)
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
      call calc_vloc_rho
      if(lcpmd) then
         psi(n1:n2,1:MSTOCC)=psi_c(n1:n2,1:MSTOCC)
         call Total_Energy_cpmd
      else
         call scf(err,err2)
         psi_c(n1:n2,1:MSTOCC)=psi(n1:n2,1:MSTOCC)
         call Total_Energy
      endif
      call force_ion
      if(myrank==0) then
         write(*,'("Forcesum=",3F15.11)') sum(Force(1,1:Mi)),
     &                 sum(Force(2,1:Mi)),sum(Force(3,1:Mi))
      endif
      return
      end subroutine eforce_cpmd
c-----------------------------------------------------------------------
c orbital rotation after coordinate update
c-----------------------------------------------------------------------
      subroutine rotorb
      use global_variables
      implicit none
      real(8),parameter   :: eps=1.0d-8
      integer,parameter   :: maxit=100
      real(8) error
      integer i,j,ii,it,ierr
      integer n1,n2,ML0
      allocate(tau(MSTOCC,MSTOCC),sig(MSTOCC,MSTOCC))
      allocate(gam(MSTOCC,MSTOCC),gamn(MSTOCC,MSTOCC))
      allocate(scr(MSTOCC,MSTOCC))
      tau=0.0d0; sig=0.0d0; gam=0.0d0; gamn=0.0d0; scr=0.0d0
c
      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)
c
      scr=0.0d0
      call dsyr2k('u','t',MSTOCC,ML0,-dV*0.5d0,psi_c(n1,1),ML0,
     &            psi_n(n1,1),ML0,0.0d0,scr,MSTOCC)
      call mpi_allreduce
     &     (scr,tau,MSTOCC**2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      scr=0.0d0
      call dsyrk
     &     ('u','t',MSTOCC,ML0,-dV,psi_n(n1,1),ML0,0.0d0,scr,MSTOCC)
      call mpi_allreduce
     &     (scr,sig,MSTOCC**2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      scr=0.0d0
      if(myrank==0) then
         do i=1,MSTOCC
            do j=i+1,MSTOCC
               tau(j,i)=tau(i,j)
            enddo
            tau(i,i)=1.0d0+tau(i,i)
         enddo
         do i=1,MSTOCC
            do j=i+1,MSTOCC
               sig(j,i)=sig(i,j)
            enddo
            sig(i,i)=1.0d0+sig(i,i)
         enddo
         gam=sig*0.5d0
         call dgemm('t','n',MSTOCC,MSTOCC,MSTOCC,0.5d0,tau,
     &                      MSTOCC,tau,MSTOCC,0.5d0,sig,MSTOCC)
         it=0
  100    continue
         it=it+1
         gamn=sig
         scr=tau-gam
         call dgemm('t','n',MSTOCC,MSTOCC,MSTOCC,-0.5d0,scr,
     &                      MSTOCC,scr,MSTOCC,1.0d0,gamn,MSTOCC)
         error=0.0d0
         scr=gamn-gam
         gam=gamn
         error=sum(scr(1:MSTOCC,1:MSTOCC)**2)
         do i=1,MSTOCC
            do j=i+1,MSTOCC
               gam(i,j)=0.5d0*(gam(i,j)+gam(j,i))
               gam(j,i)=gam(i,j)
            enddo
         enddo
         error=sqrt(error)
         if(error.gt.eps.and.it.le.maxit) goto 100
      endif
      call mpi_bcast(gam,MSTOCC**2,mpi_real8,0,mpi_comm_world,ierr)
      do i=1,MSTOCC
         do j=1,MSTOCC
            psi_n(n1:n2,i)=psi_n(n1:n2,i)+psi_c(n1:n2,j)*gam(j,i)
         end do
      end do
c
      deallocate(tau,sig,gam,gamn,scr)
      return
      end subroutine rotorb
c-----------------------------------------------------------------------
c fictious kinetic energy
c-----------------------------------------------------------------------
      subroutine calfke(fke)
      use global_variables
      implicit none
      integer :: ierr,n1,n2
      real(8) :: fke,fke0
c
      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      fke0=0.0d0 ; fke=0.0d0
      fke0=sum(psi_v(n1:n2,1:MSTOCC)**2)
      call mpi_allreduce(fke0,fke,1,mpi_real8,mpi_sum,mpi_comm_world,
     &                   ierr)
      fke=fke*dV*emass
      if(myrank==0) write(*,*) "fke=",fke
c
      return
      end subroutine calfke
c-----------------------------------------------------------------------
c Total energy for CPMD
c-----------------------------------------------------------------------
      subroutine Total_Energy_cpmd
      use global_variables
      implicit none
      integer :: i,ii,ia,ib,a,n1,n2,ierr
      real(8) :: c,s(5),r(5),e,Eloc,Eloc_ion
      real(8) :: ctime0,ctime1,etime0,etime1
c
      call watch(ctime0,etime0)
      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)
c
      s=0.0d0; r=0.0d0
c
c --- density update ---
c
      rho(n1:n2)=0.0d0
      do i=1,MSTOCC
         rho(n1:n2)=rho(n1:n2)+psi_c(n1:n2,i)*psi_c(n1:n2,i)*2.0d0
      end do
      call mpi_allgatherv(rho(n1),ircnt(myrank),mpi_real8,rho
     &     ,ircnt,idisp,mpi_real8,mpi_comm_world,ierr)
c
c --- Hartree energy ---
c
      call Hartree(rho,Vh,1.d-25,2000,2*ifac_vloc)
c
c --- Exchange-Correlation Energy ---
c
      call Exc_Cor
c
c --- Local potential ---
c
      Vloc(n1:n2)=Vion(n1:n2)+Vh(n1:n2)+Vxc(n1:n2)
      call mpi_allgatherv(Vloc(n1),ircnt(myrank),mpi_real8
     &    ,Vloc,ircnt,idisp,mpi_real8,mpi_comm_world,ierr)
c
c --- Kinetic, Non-local, etc. ---
c
      Ekin=0.d0 ; Enl=0.d0 ; Eloc=0.d0 ; Eloc_ion=0.d0

      do i=1,MSTOCC
         e      = 0.d0
         ekinsp = 0.d0
         enlsp  = 0.d0
         call hpsi_spe(i,e,n1,n2,10)
         Ekin = Ekin + ekinsp*2.0d0
         Enl  = Enl  + enlsp*2.0d0
      end do
c
c --- parts of total energy ---
c
      s(1)=sum(rho(n1:n2)*Vion(n1:n2))*dV
      s(2)=sum(rho(n1:n2)*Vloc(n1:n2))*dV
      s(3)=Ekin
      s(4)=Enl

      call mpi_allreduce(s,r,4,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      Eloc_ion=r(1)
      Eloc    =r(2)
      Ekin    =r(3)
      Enl     =r(4)

      Etot = Ekin + Enl + Eloc_ion + E_Hartree + Ex + Ec + Eion

      call watch(ctime1,etime1)

      if(myrank==0) then
!      if( myrank==0.AND.DISP_SWITCH )then
         write(*,*) "TIME(Total_Energy)=",ctime1-ctime0,etime1-etime0
         write(*,*) "(II) ",Eion,"(LI) ",Eloc_ion
         write(*,*) "(LO) ",Eloc,"(NL) ",Enl
         write(*,*) "(KI) ",Ekin,"(HT) ",E_Hartree
         write(*,*) "(EX) ",Ex,  "(CO) ",Ec
         write(*,*) "(TO) ",Etot
      end if
c
      return
      end subroutine Total_Energy_cpmd
c-----------------------------------------------------------------------
      subroutine calkin2(ekin2)
      use global_variables
      implicit none
      real(8),intent(inout) :: ekin2
      integer i
      real(8) pm
      ekin2=0.0d0
      do i=1,Mi
         pm=kmass(iatom(i))
         ekin2=ekin2
     &         +(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
      enddo
      ekin2=ekin2*0.5d0*emass*4.0d0/3.0d0
      return
      end
c-----------------------------------------------------------------------
      subroutine settherm2(bathe,kine,fke,dt2)
      use global_variables
      implicit none
      real(8) tunit,bathe,kine,fke,dt2
      real(8),parameter :: waver=1000.0d0
      real(8),parameter :: wavee=10000.0d0
      call calkin2(etemp)
      tunit=tojoul/(clight*planck*wavee)
      qe=sum(occ)*tunit*tunit
      xe=0.0d0
!      ve=sqrt(etemp/qe)
      ve=0.0d0
      fe=2.0d0*(fke-etemp)/qe
      if(ve.lt.0.0d0) then
         fe=0.0d0
      elseif(ve.eq.0.0d0) then
         fe=fe*0.5d0
      endif
      xe(-1)=xe(0)-ve*dt      +fe*dt2*0.5d0
      xe(-2)=xe(0)-ve*dt*2.0d0+fe*dt2*0.5d0*4.0d0
      xe(-3)=xe(0)-ve*dt*3.0d0+fe*dt2*0.5d0*9.0d0
      tunit=tojoul/(clight*planck*waver)
      qr=3.0d0*Mi*tunit*tunit
      xr=0.0d0
      vr=sqrt(fnkt/qr)
      fr=(2.0d0*kine-fnkt)/qr
      xr(-1)=xr(0)-vr*dt      +fr*dt2*0.5d0
      xr(-2)=xr(0)-vr*dt*2.0d0+fr*dt2*0.5d0*4.0d0
      xr(-3)=xr(0)-vr*dt*3.0d0+fr*dt2*0.5d0*9.0d0
      bathe=xr(0)*fnkt-etemp
      bathe=bathe+0.5d0*(ve*ve*qe+vr*vr*qr)
      write(*,'(a,2f15.8,a,f8.5)') "qe,qr=",qe,qr,"bathe=",bathe
      write(*,'(a,2f15.8,a,f8.5)') "fe,fr=",fe,fr,"bathe=",-etemp
      return
      end subroutine settherm2
