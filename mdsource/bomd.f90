!-----------------------------------------------------------------------
!  Main routine for molecular dynamics
!-----------------------------------------------------------------------
SUBROUTINE bomd

  use atom_module, only: Natom,ki_atom
  use cpmd_variables, only: disp_switch,Velocity,psi_v,psi_n,KB,dt &
                           ,Force,lcpmd,lquench,lbere,lbathnew,nstep &
                           ,inivel,linitnose,lmeta,dtsuz,lbath,lbathnewe &
                           ,lscaleele,Rion,Rion0,lscale,linitnosee,lblue,AMU,pmass &
                           ,iatom,TOJOUL,deltat,FS_TO_AU,dsettemp,temp,MI &
                           ,MB_0_CPMD,MB_1_CPMD,MB_0_SCF,MB_1_SCF,batm,itime
  use parallel_module
  use total_energy_module
  use wf_module
  use watch_module
  use array_bound_module, only: MB_0,MB_1
  use cpmdio2_module

  use rotorb_module
  use blue_moon_module

  implicit none

  integer :: i,j,k,n,s,ierr,ib1,ib2
  real(8) :: tott,dif,kine,tote,tote0,dt2,ltemp
  real(8) :: fke,bathe,engconf,enose
  real(8) :: ctime0,ctime1,etime0,etime1
  real(8) :: ctime_cpmd(0:11),etime_cpmd(0:11)
  real(8),allocatable :: mu(:)
  logical :: ltime,flag_etlimit
  real(8) :: sctot
  real(8) :: vcmio(4)
!------parameter for new nose
  integer,parameter :: nit=1
  integer,parameter :: ncalls=9
!------------------------------------
  integer :: MI3

  integer,parameter :: unit_trjxyz = 90

  lblue = .false.

  if ( myrank == 0 ) write(*,'(1x,a60," CPMD")') repeat("-",60)

  if ( myrank == 0 ) then
     call read_cpmd_variables ! in 'alloc_cpmd.f90'
     if ( lblue ) call read_blue
  end if
  call send_cpmd_variables ! in 'alloc_cpmd.f90'
  if ( lblue ) call bcast_blue_data

  tote0       = 0.d0
  enose       = 0.d0
  dsettemp    = temp
  MI          = Natom
  MI3         = 3*Natom
  disp_switch = (myrank==0)
  ltime       = .true.
  MB_0_SCF    = MB_0
  MB_1_SCF    = MB_1
  flag_etlimit= .false.

  if ( myrank == 0 ) then
     open(unit_trjxyz,file='TRAJECTORY.xyz')
     if ( lcpmd ) then
        if ( lbath ) then
           write(*,*) "start NVT cpmd"
        else if ( lbathnew ) then
           write(*,*) "start NVT with new Nose-Hoover"
        else
           write(*,*) "start NVE cpmd"
        endif
     else
        if ( lbath ) then
           write(*,*) "start NVT bomd"
        else
           write(*,*) "start NVE bomd"
        endif
     endif
  endif

!
!     setup variables
!
  dt        = deltat*FS_TO_AU ! time in AU
  dt2       = dt*0.50d0       ! dt/2
  temp      = temp*KB/TOJOUL  ! temperature in AU

!
!     allocate MD variables
!
  call alloc_md   ! in 'alloc_cpmd.f90'
  call asign_atom ! in 'asign_atom.f90'

!
!     initial velocity
!
  if ( myrank == 0 ) then
     if ( inivel ) then
        write(*,*) "generate initial velocity"
        call setv ! in 'setv.f90'
        call vcom(kine) ! in 'vcom.f90'
     else
        write(*,*) "read initial coordinate and velocity"
        call mdio(0,tote0) ! in 'mdio.f90'
        call calkin(kine) ! in 'calkin.f90'
     endif
     call mdio(2,tote0)
  endif

!
!     broadcast initial data
!
  if ( .not.inivel ) then
     call mpi_bcast(Rion,MI3,mpi_real8,0,mpi_comm_world,ierr)
     call mpi_bcast(tote0,1,mpi_real8,0,mpi_comm_world,ierr)
  endif
  call mpi_bcast(Velocity,MI3,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(kine,1,mpi_real8,0,mpi_comm_world,ierr)

!
!     initial wave function
!
  if ( lcpmd ) then
     call active_band ! in 'active_band.f90'
     MB_0 = MB_0_CPMD
     MB_1 = MB_1_CPMD
     ib1  = MB_0_CPMD
     ib2  = MB_1_CPMD
     call alloc_cpmd ! in 'bomd.f90' (this file)
     if ( lquench ) then
        if ( myrank == 0 ) write(*,*) "Quench system to BO surface!!!!!!!!"
        call getforce ! in 'getforce.f90'
        psi_v(:,:,:,:)=0.0d0
        psi_n(:,:,:,:)=0.0d0
        fke=0.0d0
     endif
     if(.not. inivel) then
        call read_data_cpmd_k_para
!        call gather_wf_data(1)   ! In the case of band para
     else
         continue
!        write(*,*) "point1!!!!!!!!!"
!        call gather_wf_data(0)   ! In the case of band para
!        write(*,*) "point2!!!!!!!!!"
     endif
     call calfke(fke) ! in 'calfke.f90'
     if ( myrank == 0 ) write(*,*) "calfke=",fke
  else
     fke=0.0d0
  endif

!
!     initial force
!
  if ( lcpmd ) then
     if ( inivel ) then
        call getforce ! in 'getforce.f90'
        do i=1,Natom
           Force(:,i)=Force(:,i)/(pmass(iatom(i))*AMU)
        enddo
        call wf_force ! in 'wf_force.f90'
     else
        if ( myrank == 0 ) then
           write(*,*) "energy and force evaluation form previous data"
        end if
        disp_switch=(myrank==0)
        call getforce_cpmd(.true.,ltime) ! in 'getforce_cpmd.f90'
        call wf_force ! in 'wf_force.f90'
        call calc_total_energy(.false.,disp_switch)
        disp_switch=.false.
     endif
  else
     call getforce ! in 'getforce.f90'
  endif

!
!     initial bath parameters
!
  if ( lbath ) then
     call sicoef ! in 'mnhc.f90'
     call setmnhc(kine) ! in 'mnhc.f90'
     call emnhc(bathe) ! in 'mnhc.f90'
  else
     bathe=0.0d0
  endif

!------------------------------------------------------

  if ( lbathnewe .and. inivel )  then
     call nosepae(dt) ! in 'bomd.f90'
  endif

  if ( linitnosee ) then
     call nosepae(dt) ! in 'bomd.f90'
  endif

  if ( lbathnew .and. inivel ) then
     call nosepa(dt) ! in 'bomd.f90'
  endif

  if ( linitnose ) then
     call nosepa(dt) ! in 'bomd.f90'
  endif

!------------------------------------------------------

  if ( (lbathnew.and.(.not.inivel)) .and. (.not.linitnose) ) then
     call read_nose_data ! in 'bomd.f90'
  endif
  if ( (lbathnewe.and.(.not.inivel)) .and. (.not.linitnosee) ) then
     call read_nosee_data ! in 'bomd.f90'
  endif
  if ( (lbathnew.or.lbathnewe) .and. (.not.inivel) ) then
     call make_nose_time(dt) ! in 'bomd.f90'
  endif

!
!     initialize Metadynamics (interface to PLUMED)
!
  if ( lmeta ) then
!     allocate(mu(Natom))
!     mem_new = mem_new + bdreal*Natom
!     do i=1,Natom
!        mu(i)=pmass(iatom(i))*amu
!     enddo
!     CALL init_metadyn(Natom,dt,mu,mu,1,1.0D0,"plumed.dat"//char(0))
     stop "metadynamics is temporary unavailable:stop@bomd"
  else
     engconf=0.0d0
  endif

!
!     initial data output
!
  if ( inivel ) tote0=kine+Etot+fke+bathe+engconf
  if ( myrank == 0 ) then
     write(*,'(a)') "initial energy"
     write(*,'(5f15.8)') tote0,Etot,kine,fke,bathe
     open( 3,file='traj.dat',status='unknown')
     open( 4,file='info.dat',status='unknown')
     open(15,file='time.dat',status='unknown')
     if ( ltime ) then
        open(16,file='time_cpmd_loop.dat',status='unknown')
        write(16,'(6x,9(1x,a9))') "PSI_V","PSI_N","ROTORB","GETFORCE_CPMD" &
                   ,"WF_FORCE","PSI_V","ROTORB2","CALFKE","TOTAL_ENERGY"
        open(17,file='time_force_once.dat',status='unknown')
        write(17,'(9a10)') "EWALD","LOCAL&PCC","PS_NLOC","PSI_RHO" &
                   ,"HARTREE ","EXC","VLOC","FORCE","FORCE"
     end if
     dif   = 0.0d0
     tott  = 0.0d0
     ltemp = kine*TOJOUL/KB/1.5d0/dble(Natom)
     write(4,10) tott,tote0,dif,Etot,kine,fke,bathe,ltemp
  endif

!
! loop start
!
  disp_switch = .false.

  vcmio(1:4)=0.0d0

  do itime=1,nstep

     call watch(ctime0,etime0)

     if ( lbere ) then
        call berendsen(dt2) ! in 'bomd.f90'
     endif

     if ( lbathnewe ) then
        sctot=1.0d0
        call calfke(fke) ! in 'calfke.f90'
        do i=1,nit
           do j=1,ncalls
              call enosmove(fke,dtsuz(j),sctot) ! in 'bomd.f90'
           enddo
        enddo
        psi_v(:,ib1:ib2,:,:)=psi_v(:,ib1:ib2,:,:)*sctot
        if ( myrank == 0 ) write(*,*) "sctot2------>", sctot
     endif

     if ( lbathnew ) then
        do i=1,nit
           do j=1,ncalls
              call nose(dtsuz(j))
           enddo
        enddo
     endif

     if ( lbath ) call mnhc

!-------subtruct center of mass
     call comvel(vcmio,.true.) ! in 'bomd.f90'
!------------------------------

     Rion0(:,:)    = Rion(:,:)
     Velocity(:,:) = Velocity(:,:) + Force(:,:)*dt2
     Rion(:,:)     = Rion(:,:)     + Velocity(:,:)*dt

     if ( lblue ) call cpmdshake

     if ( lcpmd ) then
        call watch(ctime_cpmd(0),etime_cpmd(0))
        psi_v(:,ib1:ib2,:,:)=psi_v(:,ib1:ib2,:,:)+psi_n(:,ib1:ib2,:,:)*dt2
        call watch(ctime_cpmd(1),etime_cpmd(1))
        psi_n(:,ib1:ib2,:,:)=unk(:,ib1:ib2,:,:)+psi_v(:,ib1:ib2,:,:)*dt
        call watch(ctime_cpmd(2),etime_cpmd(2))
        call rotorb
        call watch(ctime_cpmd(3),etime_cpmd(3))
        call getforce_cpmd(.true.,ltime)
        call watch(ctime_cpmd(4),etime_cpmd(4))
        call wf_force
        call watch(ctime_cpmd(5),etime_cpmd(5))
        psi_v(:,ib1:ib2,:,:)=psi_v(:,ib1:ib2,:,:)+psi_n(:,ib1:ib2,:,:)*dt2
        call watch(ctime_cpmd(6),etime_cpmd(6))
        if ( lscaleele ) call rscve(fke)
        call rotorb2
        call watch(ctime_cpmd(7),etime_cpmd(7))
        call calfke(fke)
        call watch(ctime_cpmd(8),etime_cpmd(8))
        call calc_total_energy(.false.,disp_switch)
        call watch(ctime_cpmd(9),etime_cpmd(9))
     else
        call getforce
        call calc_total_energy(.false.,disp_switch)
     endif

     if ( lmeta ) then
!         do i=1,Natom
!            Force(:,i)=Force(:,i)*(pmass(iatom(i))*amu)
!         enddo
!         call meta_force_calculation(aaL,itime,Rion,0,0,Force,0,0,engconf)
!         do i=1,Natom
!            Force(:,i)=Force(:,i)/(pmass(iatom(i))*amu)
!         enddo
        stop "metadynamics is temporary unavailable:stop@bomd"
     endif

     Velocity(:,:) = Velocity(:,:) + Force(:,:)*dt2

!---------------------for constraint(2)----
     if ( lblue ) then
        call rattle
        call write_blue_data(itime)
     end if
!------------------------------------------

!-----------------subtruct center of mass (Nessesary for RS-CPMD)
     call comvel(vcmio,.false.)
!-----------------------------------------

     if ( lbath ) then
        call mnhc
        call emnhc(bathe)
     endif

     if ( lbathnew ) then
        do i=1,nit
           do j=1,ncalls
              call nose(dtsuz(j))
           enddo
        enddo
        call noseene(bathe)
     endif

     if ( lbere ) then
        call berendsen(dt2)
     endif

     if ( lscale ) then
        write(*,*) "lscale------->", lscale
        call scaling_velo
     endif

     call calkin(kine)

     if ( lbathnewe ) then
        sctot=1.0d0
        call calfke(fke)
        do i=1,nit
           do j=1,ncalls
              call enosmove(fke,dtsuz(j),sctot)
           enddo
        enddo
        psi_v(:,:,:,:)=psi_v(:,:,:,:)*sctot
        if ( myrank == 0 ) write(*,*) "sctot2------>", sctot
        call noseeneele(enose)
        bathe=bathe+enose
     endif

     call watch(ctime1,etime1)

     if ( myrank == 0 ) then

        tott  = deltat*itime
        tote  = kine+Etot+fke+bathe+engconf
        ltemp = kine/KB*TOJOUL/1.5d0/dble(Natom)
        dif   = abs(tote-tote0)
!-------------------------------------------trajectry
        if ( mod(itime,1) == 0 ) then
           do i=1, Natom
              write(3,'(6f18.9)') Rion(1:3,i),Velocity(1:3,i)
           enddo
        endif
!-------------------------------------------
        if ( mod(itime,1) == 0 ) then
           write(unit_trjxyz,*) Natom
           write(unit_trjxyz,*) "CPMD on RSDFT STEP->",itime
           do i=1,Natom
              write(unit_trjxyz,'(a2,3f14.6)') &
                   batm(iatom(i)),Rion(1:3,i)*0.529177210d0
           enddo
        endif
!-------------------------------------------
        write(*,'(a,x,f8.2,x,a,x,f15.8,x,a,x,f6.1,x,a,x,f15.8)') &
             "t=",tott,"dif=",dif,"ltemp=",ltemp,"fke=",fke
        write(*,10) tott,tote,dif,Etot,kine,fke,bathe,ltemp,Sum(esp),enose
        write(*,*) Sum(occ), Sum(esp)

        write(4,10) tott,tote,dif,Etot,kine,fke,bathe,ltemp,Sum(esp),enose
        write(15,'(i6,2f20.5)') itime,ctime1-ctime0,etime1-etime0
        if ( ltime ) then
           write(16,'(i6,9f10.5)') itime,(etime_cpmd(k+1)-etime_cpmd(k),k=0,8)
        end if
     endif

     call global_watch(.false.,flag_etlimit)
     if ( flag_etlimit ) then
        if ( myrank == 0 ) write(*,*) "elapsed time limit exceeded : flag_etlimit=",flag_etlimit
        exit
     end if

  enddo ! itime
!
! loop end
!

  if ( myrank==0 ) call mdio(1,tote0)
  if ( lcpmd ) then
!     call write_data2(0)
     if (.not. allocated(unk)) write(*,*) "UNK IS NOT ALLOCATED!!!!!"
     if ( lbathnew  ) call write_nose_data
     if ( lbathnewe ) call write_nosee_data
     call write_data_cpmd_k_para
     call dealloc_cpmd
  endif
  if ( lbath ) call unsetmnhc(bathe)
  if ( lmeta ) then
     deallocate(mu)
  endif
  call dealloc_md
  close(3)
  close(4)
  close(15)
  close(unit_trjxyz)
  if ( lblue ) close(889)
  if ( ltime ) then
     close(16)
     close(17)
  end if

98 continue

  return

99 stop "stop@bomd(99)"

10 format(10f15.8)

END SUBROUTINE bomd



!kk Rescaling temp.-------06-06-2012-----------------------------

!----------------------------------------------------------------
subroutine scaling_velo
!   use global_variables
  use cpmd_variables
  use atom_module, only: Natom
  implicit none
  real(8) :: kine,settmp,pm,scale,velosum
  integer i

  kine=0.0d0
  do i=1,Natom
     pm=pmass(iatom(i))*amu
     kine=kine+(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
  enddo
  kine=kine*0.5d0

  kine=kine/KB*TOJOUL/1.5d0/dble(Natom)

  if ( myrank == 0 ) write(*,*) "SCALE KINETIC ENERGY----temp---->",kine

  if(((kine.gt.dsettemp+trange).or.&
       &(kine.lt.(dsettemp-trange))).and.(kine.ne.0.0d0)) then
     scale=DSQRT(dsettemp/kine)
     do i=1,Natom
        Velocity(1,i)=Velocity(1,i)*scale
        Velocity(2,i)=Velocity(2,i)*scale
        Velocity(3,i)=Velocity(3,i)*scale
     enddo
  endif

  velosum=sum(Velocity)

  return
end subroutine scaling_velo

!---------------------------------------------------------------
!---------------set Nose-Hoover parameters for ION  KK-----------
subroutine nosepa(delt_elec)
!        use global_variables
  use cpmd_variables, only: omegan,dtsuz,dsettemp,Mi,nch,gkt,qnospc,etap1,pi,etap1dot
  use atom_module, only: Natom
  implicit none
  real(8) :: w_yosh7_1,w_yosh7_2,w_yosh7_3,w_yosh7_4,w_yosh7_0
  real(8) :: rnr(3)
  real(8) :: wntau
  parameter (wntau=7.26d-7)
  real(8),parameter :: seed=135792468D0
  real(8),parameter :: ry=13.60569193D0
  real(8),parameter :: evtokel=11604.505D0
  integer,parameter :: ip=1
  integer,parameter :: nit=1
  real(8) :: wnosp0,wnosep,wnosep2,factem
  integer :: ipp,l
  real(8) :: tempw,alfa1,sigma,dofst,delt_elec,NOSPT0
  real*8  xranf
  external xranf



!-------wnosp0 Nose freq. ------------------------
  WNOSP0=omegan
!-------------------------------------------------

  wnosep=wnosp0*wntau
  wnosep2=wnosep*wnosep

  w_yosh7_1 = 0.192d0
  w_yosh7_2 =  0.554910818409783619692725006662999D0
  w_yosh7_3 =  0.124659619941888644216504240951585D0
  w_yosh7_4 =  -0.843182063596933505315033808282941D0
  w_yosh7_0 = 1.0d0 - 2.0d0*(w_yosh7_1+w_yosh7_2+  &
       &             w_yosh7_3+w_yosh7_4)
  dtsuz(1) = w_yosh7_4*delt_elec/dble(nit)
  dtsuz(2) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(3) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(4) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(5) = w_yosh7_0*delt_elec/dble(nit)
  dtsuz(6) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(7) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(8) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(9) = w_yosh7_4*delt_elec/dble(nit)

!----------------------------------
  factem=2.d0*ry*evtokel
  ipp=1
  tempw=dsettemp
  NOSPT0=TEMPW
!----------------------------------
!-------------for test-------------
  DOFST=3.0d0*(Natom-1)
!----------------------------------

  if(ipp.gt.1) then
     gkt(ipp,1) = 3.d0*(Mi-1)*tempw/factem
     qnospc(1,ip) = gkt(ipp,1)/wnosep2
     do l=2,nch
        qnospc(l,ip) = gkt(ipp,1)/wnosep2/(3.d0*Mi)
     enddo
  else
     gkt(ipp,1) = DOFST*tempw/factem
     qnospc(1,ip) = gkt(ipp,1)/wnosep2
     do l=2,nch
        qnospc(l,ip) = gkt(ipp,1)/wnosep2/DOFST
     enddo
  endif

  do l=1,nch
     etap1(l,ip) = 0.0d0
     sigma=dsqrt(NOSPT0/qnospc(l,ip)/factem)
     call mprand(seed,3,rnr)
     alfa1=2.0d0*pi*rnr(1)
     etap1dot(l,ip)=sqrt(log(xranf())*(-2.0d0))*  &
     &                               cos(alfa1)*sigma
  enddo
  return
end subroutine nosepa

!---------------------------------------------------------

subroutine mprand(seed,n,a)
  implicit none
  integer n
  real*8  a(n),seed
  real*8  xranf
  external xranf
  integer i

  do i=1,n
     a(i)=xranf()
  enddo

  return
end subroutine mprand

function xranf()
  implicit none
  real*8 xranf
  integer m,konst
  data    m/100001/,konst/125/
  save    m
  m=m*konst
  m=m-2796203*(m/2796203)
  xranf=dble(m)/2796203.D0
  return
end function xranf

!----------------------------------------------------------

!----------------------------------------------------------
!--------------Nose'-Hoover chain for ION (same as CPMD)

subroutine nose(step)
!        use global_variables
  use cpmd_variables, only: dsettemp,Mi,gkt,fetapv,nch,qnospc,etap1dot,Velocity,etap1
  implicit none
  integer :: l,ipp
  real(8),parameter :: ry=13.60569193D0
  real(8),parameter :: evtokel=11604.505D0
  integer,parameter :: ip=1
  real(8) :: step,factem,facton,tempw
  real(8) :: pkewant,ekinp,aaex
  real(8) :: DOFST
!----------
  ipp=1
  factem=2.d0*ry*evtokel
  tempw=dsettemp

  facton=0.5d0/factem
  DOFST=3.0d0*(Mi-1)
!----------

  ekinp=0.d0
  call calkin(ekinp)

  if(ipp.gt.1) then
     gkt(ipp,1) = 3.d0*(Mi-1)*tempw/factem
  else
     gkt(ipp,1) = DOFST*tempw/factem
  endif

  pkewant = 0.5d0*gkt(ipp,1)
  fetapv(1) = 2.0d0*(ekinp - pkewant)
  do l=2,nch
     ekinp = 0.5D0*qnospc(l-1,ip)*etap1dot(l-1,ip)* &
          &         etap1dot(l-1,ip)
     pkewant = tempw*facton
     fetapv(l) = 2.0d0*(ekinp - pkewant)
  enddo
!     ==------------------------------------------------------------==
!     == MOVE LAST CHAIN ELEMENT                                    ==
!     ==------------------------------------------------------------==
  etap1dot(nch,ip) = etap1dot(nch,ip) +   &
       &                  0.25d0*step*fetapv(nch)/qnospc(nch,ip)
!     ==------------------------------------------------------------==
!     ==  LOOP OVER REST OF NOSE-HOOVER CHAIN                       ==
!     ==------------------------------------------------------------==
  do l=1,nch-1
     aaex = DEXP(-0.125d0*step*etap1dot(nch+1-l,ip))

     etap1dot(nch-l,ip) = etap1dot(nch-l,ip)*aaex*aaex + &
          & 0.25d0*step*fetapv(nch-l)*aaex/qnospc(nch-l,ip)

  enddo
!     ==--------------------------------------------------------------==
!     ==  DO VELOCITY SCALING AND PROPAGATE THERMOSTAT POSITIONS      ==
!     ==--------------------------------------------------------------==
  aaex = DEXP(-0.5d0*step*etap1dot(1,ip))

  Velocity=Velocity*aaex

  do l=1,nch
     etap1(l,ip) = etap1(l,ip) + 0.5d0*step*etap1dot(l,ip)
  enddo
!     ==--------------------------------------------------------------==
!     == RECALCULATE CHAIN FORCE 1                                    ==
!     ==--------------------------------------------------------------==
  ekinp=0.d0
  call calkin(ekinp)

  pkewant = 0.5d0*gkt(ipp,1)
  fetapv(1) = 2.0d0*(ekinp - pkewant)
!     ==-------------------------------------------------------------==
!     ==  LOOP OVER CHAIN AGAIN                                      ==
!     ==-------------------------------------------------------------==
  do l=1,nch-1
     aaex = DEXP(-0.125d0*step*etap1dot(l+1,ip))

     etap1dot(l,ip) = etap1dot(l,ip)*aaex*aaex + &
          & 0.25d0*step*fetapv(l)*aaex/qnospc(l,ip)

     ekinp = 0.5d0*qnospc(l,ip)*etap1dot(l,ip)*etap1dot(l,ip)
     pkewant = tempw*facton
     fetapv(l+1) = 2.0d0*(ekinp - pkewant)
  enddo
!     ==------------------------------------------------------------==
!     == MOVE LAST CHAIN ELEMENT                                    ==
!     ==------------------------------------------------------------==
  etap1dot(nch,ip) = etap1dot(nch,ip) +  &
       &                  0.25d0*step*fetapv(nch)/qnospc(nch,ip)
!     ==------------------------------------------------------------==
  return
end subroutine nose
!----------------------------------------------------------------------------

!---------------------------------------------------------------------------
!---------------------calculate energy of bath for ION 9/25 KK   -----------
!---------------------------------------------------------------------------
subroutine noseene(bathe)
!        use global_variables
  use cpmd_variables
  implicit none
  real(8),parameter :: ry=13.60569193D0
  real(8),parameter :: evtokel=11604.505D0
  integer :: ipp,ip,l
  real(8) :: bathe,factem,tempw


  ipp=1
  factem=2.d0*ry*evtokel
  tempw=dsettemp

  bathe=0.0d0

  do ip=1,1
     bathe = 0.5d0*qnospc(1,ip)*etap1dot(1,ip)*etap1dot(1,ip) &
          &  + gkt(ipp,1)*etap1(1,ip)
     do l=2,nch
        bathe = bathe + 0.5d0*qnospc(l,ip)*etap1dot(l,ip)*     &
             &   etap1dot(l,ip) + tempw*etap1(l,ip)/factem
     enddo
  enddo
  return
end subroutine noseene

!-------------------------------------------------------------------------
!----------------Restart Nose'--------------------------------------------
!-------------------------------------------------------------------------
subroutine write_nose_data
  !        use global_variables
  use cpmd_variables
  implicit none
  integer :: i
  integer,parameter :: ip=1 

  if(myrank.eq.0) then
     open(999,file="Bathdata.dat1")
     write(999,*)(etap1(i,ip),etap1dot(i,ip),i=1,nch) 
     write(999,*) gkt(1,1)
     write(999,*)(qnospc(i,ip),i=1,nch) 
     write(999,*)(fetapv(i),i=1,nch)
     close(999)
  endif
  return
end subroutine write_nose_data

subroutine read_nose_data
!        use global_variables
  use cpmd_variables
  implicit none
  integer :: i,ierr
  integer,parameter :: ip=1

  if(myrank.eq.0) then
     open(999,file="Bathdata.dat")
     read(999,*)(etap1(i,ip),etap1dot(i,ip),i=1,nch) 
     read(999,*) gkt(1,1)
     read(999,*)(qnospc(i,ip),i=1,nch)
     read(999,*)(fetapv(i),i=1,nch)
     close(999)
  endif

  call mpi_bcast(etap1(1,ip),nch,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(etap1dot(1,ip),nch,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(gkt(1,1),1,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(qnospc(1,ip),nch,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(fetapv(1),nch,mpi_real8,0,mpi_comm_world,ierr)

  return
end subroutine read_nose_data


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!-------Berendsen thermostat for ions 08/16 KK---------------------------------

subroutine berendsen(dt2)
!        use global_variables
  use cpmd_variables
  implicit none
  real(8) :: kine,settmp,pm,lambda,dt2
  integer i
  real(8) :: thresh,taubp
  parameter(thresh=1.0d-6)  ! don't scale if we are too "cold"
  parameter(taubp=2.0d0)    ! time constant

  kine=0.0d0
  do i=1,Mi
     pm=pmass(iatom(i))*amu
     kine=kine+(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
  enddo
  kine=kine*0.5d0

  kine=kine/KB*TOJOUL/1.5d0/dble(Mi)

  if((kine/dsettemp).gt.thresh) then
     lambda=DSQRT(1.0d0+dt2/taubp*((dsettemp/kine)-1.0d0))
     do i=1,Mi
        Velocity(1,i)=Velocity(1,i)*lambda
        Velocity(2,i)=Velocity(2,i)*lambda
        Velocity(3,i)=Velocity(3,i)*lambda
     enddo
  endif

  return
end subroutine berendsen

!-----------------------------------------------------------------------------
!---------------Scaling fictitious kinetic energy KK--------------------------
subroutine rscve(fke)
  use cpmd_variables
  implicit none
  real(8) :: alfae,fke
  !       if(fke.gt.ekin1.or.fke.lt.ekin2.and.fke.ne.0.d0) then
  if ( fke > ekin1 .and. fke /= 0.d0 ) then
  !   write(*,*) "SCALE ELEC. KIN. ENE."
     alfae=dsqrt(ekinw/fke)
     psi_v=psi_v*alfae
  endif
  return
end subroutine rscve


!--------------------------------------------------------------------
!-------------determine Nose parameters for electrons---10/1 KK------
!----------------under construction----------------------------------
subroutine nosepae(delt_elec)
!        use global_variables
  use cpmd_variables
  implicit none
  real(8) :: wntau
  real(8) :: omega,wnosee2,delt_elec 
  integer :: ngws, i, n
  parameter (wntau=7.26d-7)
  real(8) :: w_yosh7_1,w_yosh7_2,w_yosh7_3,w_yosh7_4,w_yosh7_0
  integer,parameter :: nit=1
  integer,parameter :: ip=1
  integer,parameter :: nedof0=6

!        write(*,*) "wnose0-->", wnose0
!        write(*,*) "ekinw--->", ekinw

  wnosee=wnose0*wntau
  n=mstocck(1,1)

!        omega=dV*ML1*ML2*ML3
!        ngws=int((1.0d0/(2.0d0*3.1415d0*3.1415d0))*omega*(Ecut**(3/2)))
!        write(*,*) "plane wave number->",ngws

!        n1=n*(2*ngws-1)
!        nc=n*n
!        ntrue=(n1-nc)

!       if NEDOF<0 use true number of DOFs
!       We must determine only nedof !!!!!!!!!!
!        nedof=n*ML1*ML2*ML3

!        nedof=ngws
!        nedof=2592

  nedof=int(n*nedof0)

!        write(*,*) "n------->", n
!        write(*,*) "nedof--->", nedof


  do i=1,nche
     qnosee(i)=0.0d0
  enddo

  wnosee2 = wnosee*wnosee

  qnosee(1) = 2.0d0*ekinw/wnosee2
  do i=2,nche
     qnosee(i) = 2.d0*ekinw/wnosee2/dble(nedof)
  enddo

  do i=1,nche
     etadot(i,ip)=0.0d0
  enddo

  etadot(1,ip) = dsqrt(2.d0*ekinw/qnosee(1))
  do i=2,nche
     etadot(i,ip) = dsqrt(2.d0*ekinw/qnosee(i)/dble(nedof))
  enddo

  w_yosh7_1 = 0.192d0
  w_yosh7_2 =  0.554910818409783619692725006662999D0
  w_yosh7_3 =  0.124659619941888644216504240951585D0
  w_yosh7_4 =  -0.843182063596933505315033808282941D0
  w_yosh7_0 = 1.0d0 - 2.0d0*(w_yosh7_1+w_yosh7_2+  &
       &             w_yosh7_3+w_yosh7_4)
  dtsuz(1) = w_yosh7_4*delt_elec/dble(nit)
  dtsuz(2) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(3) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(4) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(5) = w_yosh7_0*delt_elec/dble(nit)
  dtsuz(6) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(7) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(8) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(9) = w_yosh7_4*delt_elec/dble(nit)

  return
end subroutine nosepae
!-----------------------------------------------------------------------------
!----------------Nose'-Hoover chain for electrons (same as CPMD) 10/1 KK------
!-----------------------------------------------------------------------------
subroutine enosmove(ekinc,step,sctot)
!      use global_variables
  use cpmd_variables
  implicit none
  real(8) ::  ekinc,step,sctot
  real(8) ::  ckewant,ckine,aae,f1,f2
  integer ::  l
  integer,parameter :: ip=1
!     ==--------------------------------------------------------------==
!     ==  COMPUTE ALL THE THERMOSTAT FORCES                           ==
!     ==--------------------------------------------------------------==
!------------------------find bug--------------------
!      write(*,*) "ekinc, ekinw----->", ekinc, ekinw

  ckewant = ekinw/dble(nedof)
  feta(1) = 2.d0*(ekinc - ekinw)
  do l=2,nche+1
     ckine = 0.5d0*qnosee(l-1)*etadot(l-1,ip)*etadot(l-1,ip)
     feta(l) = 2.d0*(ckine - ckewant)
  enddo
!     ==--------------------------------------------------------------==
!     ==  ETADOT(NCHE) IS TREATED SEPARATELY                          ==
!     ==--------------------------------------------------------------==
  aae = exp(-0.125d0*step*etadot(nche-1,ip))
  etadot(nche,ip) = etadot(nche,ip)*aae*aae + &
       &               0.25d0*step*feta(nche)*aae/qnosee(nche)
  ckine = 0.5d0*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
  feta(nche+1) = 2.d0*(ckine - ckewant)
  f1 = feta(nche-1)
  f2 = feta(nche+1)
  feta(nche-1) = f1 + f2
!     ==--------------------------------------------------------------==
!     ==  LOOP BACKWARDS OVER CHAIN                                   ==
!     ==--------------------------------------------------------------==
  do l=1,nche-1
     aae = exp(-0.125d0*step*etadot(nche+1-l,ip))
     etadot(nche-l,ip) = etadot(nche-l,ip)*aae*aae + &
          &                  0.25d0*step*feta(nche-l)*aae/qnosee(nche-l)
  enddo
!     ==--------------------------------------------------------------==
!     ==  CALCULATE SCALING FACTOR AND APPLY TO ELECTRON KE           ==
!     ==  AND ACCUMULATE TOTAL SCALING FACTOR IN SCTOT                ==
!     ==--------------------------------------------------------------==
  aae = exp(-0.5d0*step*etadot(1,ip))
  sctot = sctot*aae
  ekinc = ekinc*aae*aae
!     ==--------------------------------------------------------------==
!     ==  PROPAGATE THERMOSTAT POSITIONS                              ==
!     ==--------------------------------------------------------------==
  do l=1,nche
     etap(l,ip) = etap(l,ip) + 0.5d0*step*etadot(l,ip)
     feta(l) = 0.d0
  enddo
!     ==--------------------------------------------------------------==
!     ==  COMPUTE THERMOSTAT FORCES FOR 1 AND NCHE-1                  ==
!     ==  (THE REST ARE COMPUTED IN THE PROPAGATION LOOP)             ==
!     ==--------------------------------------------------------------==
  feta(1) = 2.d0*(ekinc - ekinw)
  ckine = 0.5d0*qnosee(nche)*etadot(nche,ip)*etadot(nche,ip)
  feta(nche-1) = 2.d0*(ckine - ckewant)
!     ==--------------------------------------------------------------==
!     ==  LOOP FORWARDS OVER CHAIN                                    ==
!     ==--------------------------------------------------------------==
  do l=1,nche-1
     aae = exp(-0.125d0*step*etadot(l+1,ip))
     etadot(l,ip) = etadot(l,ip)*aae*aae + &
          &                                 0.25d0*step*feta(l)*aae/qnosee(l)
     ckine = 0.5d0*qnosee(l)*etadot(l,ip)*etadot(l,ip)
     feta(l+1) = feta(l+1) + 2.d0*(ckine - ckewant)
  enddo
  aae = exp(-0.125d0*step*etadot(nche-1,ip))
  etadot(nche,ip) = etadot(nche,ip)*aae*aae +  &
       &               0.25d0*step*feta(nche)*aae/qnosee(nche)
!     ==-------------------------------------------------------------==
  return
end subroutine enosmove

!----------------------------------------------------------------------
subroutine noseeneele(enose)
  !      use global_variables
  use cpmd_variables
  implicit none
  real(8) :: enose
  integer :: i
  integer,parameter :: ip=1

  enose=0.0d0
  enose = enose + 0.5d0*qnosee(1)*etadot(1,ip)*etadot(1,ip) + &
       &                   2.0d0*ekinw*etap(1,ip)
  do i=2,nche
     enose = enose + 0.5d0*qnosee(i)*etadot(i,ip)*etadot(i,ip) + &
          &                   2.0d0*ekinw*etap(i,ip)/dble(nedof)
  enddo
  enose = enose + 2.d0*ekinw*etap(nche-1,ip)/dble(nedof)
  return
end subroutine noseeneele
!----------------------------------------------------------------------
subroutine write_nosee_data
  !      use global_variables
  use cpmd_variables
  implicit none
  integer :: i
  integer,parameter :: ip=1

  if(myrank.eq.0) then
     open(999,file="Bathdata_ele.dat1")
     write(999,*)(etap(i,ip),etadot(i,ip),i=1,nche) 
     write(999,*)(feta(i),i=1,nche+1) 
     write(999,*)(qnosee(i),i=1,nche) 
     write(999,*) nedof
     close(999)
  endif
  return
end subroutine write_nosee_data


subroutine read_nosee_data
!      use global_variables
  use cpmd_variables
  implicit none
  integer :: i, ierr
  integer,parameter :: ip=1

  if(myrank.eq.0) then
     open(999,file="Bathdata_ele.dat")
     read(999,*)(etap(i,ip),etadot(i,ip),i=1,nche) 
     read(999,*)(feta(i),i=1,nche+1)
     read(999,*)(qnosee(i),i=1,nche) 
     read(999,*) nedof
     close(999)
  endif

  call mpi_bcast(etap(1,ip),nche,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(etadot(1,ip),nche,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(feta(1),nche+1,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(qnosee(1),nche,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(nedof,1,mpi_integer,0,mpi_comm_world,ierr)
      
  return
end subroutine read_nosee_data
!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------
subroutine make_nose_time(delt_elec)
!        use global_variables
  use cpmd_variables
  implicit none
  real(8) :: delt_elec
  real(8) :: w_yosh7_1,w_yosh7_2,w_yosh7_3,w_yosh7_4,w_yosh7_0
  integer,parameter :: nit=1

  w_yosh7_1 = 0.192d0
  w_yosh7_2 =  0.554910818409783619692725006662999D0
  w_yosh7_3 =  0.124659619941888644216504240951585D0
  w_yosh7_4 =  -0.843182063596933505315033808282941D0
  w_yosh7_0 = 1.0d0 - 2.0d0*(w_yosh7_1+w_yosh7_2+  &
       &             w_yosh7_3+w_yosh7_4)
  dtsuz(1) = w_yosh7_4*delt_elec/dble(nit)
  dtsuz(2) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(3) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(4) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(5) = w_yosh7_0*delt_elec/dble(nit)
  dtsuz(6) = w_yosh7_1*delt_elec/dble(nit)
  dtsuz(7) = w_yosh7_2*delt_elec/dble(nit)
  dtsuz(8) = w_yosh7_3*delt_elec/dble(nit)
  dtsuz(9) = w_yosh7_4*delt_elec/dble(nit)

  return
end subroutine make_nose_time
!-------------------------------------------------------------------
subroutine comvel(vcmio,inout)
  !      use global_variables
  use cpmd_variables
  implicit none
  real(8) :: kine,pm,temp1,temp2
  real(8) :: tscal,vscale,totmass
  real(8) :: vcmio(4),vcm(3)
  integer :: i,j
  logical :: inout
!calculate first temp.
  kine=0.0d0
  do i=1,Mi
     pm=pmass(iatom(i))*amu
     kine=kine+(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
  enddo
  temp1=kine*0.5d0
  if(inout) then

     vcm(1)=0.0d0
     vcm(2)=0.0d0
     vcm(3)=0.0d0
     totmass=0.0d0

     do i=1,Mi
        vcm(1)=vcm(1)+pmass(iatom(i))*Velocity(1,i)
        vcm(2)=vcm(2)+pmass(iatom(i))*Velocity(2,i)
        vcm(3)=vcm(3)+pmass(iatom(i))*Velocity(3,i)
        totmass=totmass+pmass(iatom(i))
     enddo

     do j=1,3
        vcmio(j)=vcm(j)
     enddo
     vcmio(4)=totmass

  else

     do j=1,3
        vcm(j)=vcmio(j)
     enddo
     totmass=vcmio(4)

  endif

  do i=1,Mi
     do j=1,3
        Velocity(j,i)=Velocity(j,i)-0.5d0*vcm(j)/totmass
     enddo
  enddo
!----------------debug--------
!      if(inout) then
!         write(*,"(a3,4f15.8,i6)") 'in ',vcm(1),vcm(2),vcm(3),totmass,itime
!      else
!         write(*,"(a3,4f15.8,i6)") 'out',vcm(1),vcm(2),vcm(3),totmass,itime
!      endif
!----------------------------

!calculate second temp.
  kine=0.0d0
  do i=1,Mi
     pm=pmass(iatom(i))*amu
     kine=kine+(Velocity(1,i)**2+Velocity(2,i)**2+Velocity(3,i)**2)*pm
  enddo
  temp2=kine*0.5d0

  if(temp2.gt.1.d-5) then
     tscal=temp1/temp2
     vscale=DSQRT(tscal)
     do i=1,Mi
        Velocity(1,i)=Velocity(1,i)*vscale
        Velocity(2,i)=Velocity(2,i)*vscale
        Velocity(3,i)=Velocity(3,i)*vscale
     enddo
  endif
 
      
  return
end subroutine comvel
!-------------------------------------------------------------------







