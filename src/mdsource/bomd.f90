!-----------------------------------------------------------------------
!  Main routine for molecular dynamics
!-----------------------------------------------------------------------
SUBROUTINE bomd

  use lattice_module, only: lattice, get_aa_lattice
  use atom_module, only: Natom,ki_atom,zn_atom,aa_atom &
                        ,convert_to_xyz_coordinates_atom
  use cpmd_variables, only: disp_switch,Velocity,psi_v,psi_n,dt &
                           ,Force,lcpmd,lquench,lbere,lbathnew,nstep &
                           ,inivel,dtsuz,lbathnewe &
                           ,lscaleele,Rion,lscale,lblue &
                           ,AMU,pmass,Etot,trjstep,Ndof,omegan,ekinw,wnose0 &
                           ,deltat,FS_TO_AU,temp,MI &
                           ,MB_0_CPMD,MB_1_CPMD,MB_0_SCF,MB_1_SCF,batm,itime &
                           ,wrtstep, all_traj
  use parallel_module
  use total_energy_module
  use wf_module
  use watch_module
  use array_bound_module, only: MB_0,MB_1
  use cpmdio_module
  use rotorb_module
  use blue_moon_module
  use calc_kine_temp_module
  use velocity_scaling_module
  use velocity_scaling_ele_module
  use berendsen_module
  use nose_hoover_chain_module
  use nose_hoover_chain_ele_module

  implicit none

  integer :: i,j,k,n,s,ierr,ib1,ib2
  real(8) :: tott,dif,kine,tote,tote0,dt2,ltemp
  real(8) :: fke,Ebath_ion,Ebath_ele,Ebath
  real(8) :: ctime0,ctime1,etime0,etime1
  real(8) :: ctime_cpmd(0:11),etime_cpmd(0:11)
  real(8),allocatable :: mu(:)
  logical :: ltime,flag_etlimit
  integer,parameter :: unit_trjxyz = 90
  logical,external :: exit_program
  type(lattice) :: aa_obj

  call write_border( 0, "" )
  call write_border( 0, " CPMD START -----------" )

#ifndef _DRSDFT_
  write(*,*) "RS-CPMD is not available for COMPLEX16 WFs"
  write(*,*) "Please re-compile the program"
  call stop_program( "" )
#endif

  call check_disp_switch( .false., 1 )

  lblue = .false.

  call read_cpmd_variables ! in 'alloc_cpmd.f90'
  if ( lblue ) call read_blue

  tote0       = 0.d0
  MI          = Natom
  disp_switch = (myrank==0)
  ltime       = .true.
  MB_0_SCF    = MB_0
  MB_1_SCF    = MB_1
  flag_etlimit= .false.

  if ( wrtstep == 0 ) wrtstep=nstep+1 

  if ( myrank == 0 ) then
     open(unit_trjxyz,file='TRAJECTORY.xyz',status="replace")
     close(unit_trjxyz)
     if ( all_traj > 0 ) open(3,file='traj.dat',status='replace')
     !close(3)
     open( 4,file='info.dat',status='replace')
     open(15,file='time.dat',status='replace')
     if ( ltime ) then
        open(16,file='time_cpmd_loop.dat',status='replace')
        write(16,'(6x,9(1x,a9))') "PSI_V","PSI_N","ROTORB","GETFORCE_CPMD" &
                   ,"WF_FORCE","PSI_V","ROTORB2","CALFKE","TOTAL_ENERGY"
        open(17,file='time_force_once.dat',status='replace')
        write(17,'(9a10)') "EWALD","LOCAL&PCC","PS_NLOC","PSI_RHO" &
                   ,"HARTREE ","EXC","VLOC","FORCE","FORCE"
     end if
  end if

  if ( lbathnew ) then
     if ( disp_switch ) write(*,*) "start NVT with new Nose-Hoover"
  else
     if ( disp_switch ) write(*,*) "start NVE cpmd"
  end if

! --- setup variables

  dt  = deltat*FS_TO_AU ! time in AU
  dt2 = dt*0.50d0       ! dt/2

! --- allocate MD variables

  call alloc_md   ! Rion,Rion0,Velocity,Force

! --- Degree of fredom

  call calc_Ndof

! --- initial coordinates & velocity

  if ( inivel ) then
     if ( disp_switch ) write(*,*) "generate initial velocity"
     if ( myrank == 0 ) call setv( temp, Velocity )
     call get_aa_lattice( aa_obj )
     Rion(:,:)=aa_atom(:,:)
     call convert_to_xyz_coordinates_atom( aa_obj, Rion )
  else
     if ( disp_switch ) write(*,*) "read initial coordinate and velocity"
     if ( myrank == 0 ) call mdio( 0, tote0 )
  end if

  call mpi_bcast(Rion,size(Rion),mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(tote0,1,mpi_real8,0,mpi_comm_world,ierr)
  call mpi_bcast(Velocity,size(Velocity),mpi_real8,0,mpi_comm_world,ierr)

! --- initial wave function

  fke = 0.0d0

  if ( lcpmd ) then
     call active_band
     MB_0 = MB_0_CPMD
     MB_1 = MB_1_CPMD
     ib1  = MB_0_CPMD
     ib2  = MB_1_CPMD
     call alloc_cpmd
     if ( .not.inivel ) call read_data_cpmdio
     if ( lquench ) then
        call getforce
        psi_v=0.0d0
     else
        call getforce_cpmd( ltime )
     end if
     call wf_force
     call calc_total_energy( .false., Etot )
     call calfke( fke )
     if ( disp_switch ) write(*,*) "fictitious kinetic energy (init)=",fke
  else
     call getforce
  end if

! --- initial bath parameters

  Ebath_ele = 0.0d0
  if ( lbathnewe ) then
     call init_nose_hoover_chain_ele(dt,ekinw,wnose0,nint(sum(occ)),Ebath_ele)
  end if

  Ebath_ion = 0.0d0
  if ( lbathnew  ) then
     call init_nose_hoover_chain( dt, temp, omegan, Ebath_ion )
  end if

  Ebath = Ebath_ion + Ebath_ele

! --- initial data output

  call calc_kine( Velocity, kine, ltemp )

  if ( inivel ) tote0 = kine + Etot + fke + Ebath

  if ( myrank == 0 ) then
     write(*,'(a)') "initial energy"
     write(*,'(5a15)') "Etotal","Etot_DFT","kine","fke","Ebath"
     write(*,'(5f15.8)') tote0,Etot,kine,fke,Ebath
     dif  = 0.0d0
     tott = 0.0d0
     if ( inivel ) write(4,10) tott,tote0,dif,Etot,kine,fke,Ebath,ltemp,sum(esp),Ebath_ele
  end if

! --- loop start

  disp_switch = .false.

  do itime=1,nstep

     call watch(ctime0,etime0)

     if ( lbathnewe ) then
        call calfke( fke )
        call nose_hoover_chain_ele( fke, psi_v, MB_0_CPMD,MB_1_CPMD )
     end if

     if ( lscale ) call velocity_scaling( temp, Velocity )

     if ( lbere ) call berendsen( temp, dt2, Velocity )

     if ( lbathnew ) call nose_hoover_chain( Velocity )

     Velocity(:,:) = Velocity(:,:) + Force(:,:)*dt2

     call vcom( Velocity ) ! center of mass motion off

     if ( lblue ) then
        call shake( Rion, Velocity )
     else
        Rion(:,:) = Rion(:,:) + Velocity(:,:)*dt
     end if

     if ( lcpmd ) then

        call watch(ctime_cpmd(0),etime_cpmd(0))

        psi_v(:,ib1:ib2,:,:)=psi_v(:,ib1:ib2,:,:)+psi_n(:,ib1:ib2,:,:)*dt2

        call watch(ctime_cpmd(1),etime_cpmd(1))

        psi_n(:,ib1:ib2,:,:)=unk(:,ib1:ib2,:,:)+psi_v(:,ib1:ib2,:,:)*dt

        call watch(ctime_cpmd(2),etime_cpmd(2))

        call rotorb

        call watch(ctime_cpmd(3),etime_cpmd(3))

        call getforce_cpmd( ltime ) ! band-parallel WFs(unk) are gathered here

        call watch(ctime_cpmd(4),etime_cpmd(4))

        call wf_force

        call watch(ctime_cpmd(5),etime_cpmd(5))

        psi_v(:,ib1:ib2,:,:)=psi_v(:,ib1:ib2,:,:)+psi_n(:,ib1:ib2,:,:)*dt2

        call watch(ctime_cpmd(6),etime_cpmd(6))

        if ( lscaleele ) call velocity_scaling_ele( fke, psi_v )

        call rotorb2
        call watch(ctime_cpmd(7),etime_cpmd(7))

        call calfke(fke)

        call watch(ctime_cpmd(8),etime_cpmd(8))
        call calc_total_energy( .false., Etot )
        call watch(ctime_cpmd(9),etime_cpmd(9))

     else ! BOMD

        call getforce
        call calc_total_energy( .false., Etot )

     end if

     Velocity(:,:) = Velocity(:,:) + Force(:,:)*dt2

     if ( lblue ) then ! Blue-Moon Method
        call rattle( Rion, Velocity )
        if ( mod(itime-1,trjstep)==0 ) call write_blue_data(itime,myrank==0)
     end if

     call vcom( Velocity ) ! center of mass motion off

     if ( lbathnew ) then
        call nose_hoover_chain( Velocity, Ebath_ion )
     end if

     if ( lbere ) call berendsen( temp, dt2, Velocity )

     if ( lscale ) call velocity_scaling( temp, Velocity )

     if ( lbathnewe ) then
        call calfke( fke )
        call nose_hoover_chain_ele( fke,psi_v,MB_0_CPMD,MB_1_CPMD,Ebath_ele )
     end if
     Ebath = Ebath_ion + Ebath_ele

     call watch(ctime1,etime1)

     call calc_kine( Velocity, kine, ltemp )

     if ( myrank == 0 ) then

        tott  = deltat*itime
        tote  = kine+Etot+fke+Ebath
        dif   = abs(tote-tote0)

        write(*,'(1x,f10.3,9f15.8)') tott,tote,Etot,kine,fke,ltemp
        write(4,10) tott,tote,dif,Etot,kine,fke,Ebath,ltemp,sum(esp),Ebath_ele
        write(15,'(i6,2f20.5)') itime,ctime1-ctime0,etime1-etime0
        if ( ltime ) then
           write(16,'(i6,9f10.5)') itime,(etime_cpmd(k+1)-etime_cpmd(k),k=0,8)
        end if

        if ( mod(itime-1,trjstep) == 0 ) then
           open(unit_trjxyz,file="TRAJECTORY.xyz",position="append")
           write(unit_trjxyz,*) Natom
           write(unit_trjxyz,*) "CPMD on RSDFT STEP->",itime
           do i=1,Natom
              write(unit_trjxyz,'(a2,3f14.6)') &
                   batm(zn_atom(ki_atom(i))),Rion(1:3,i)*0.529177210d0
           end do
           close(unit_trjxyz)
        end if

        if ( all_traj > 0 ) then
           !open(3,file='traj.dat',position="append")
           do i=1,Natom
              write(3,'(3f24.16," R",i8,f10.3)') Rion(1:3,i),i,tott
              write(3,'(3f24.16," V",i8,f10.3)') Velocity(1:3,i),i,tott
              write(3,'(3f24.16," F",i8,f10.3)') Force(1:3,i),i,tott
           end do
           !close(3)
           if ( all_traj /= 1 .and. mod(itime,all_traj) == 0 ) then
              close(3)
              open(3,file='traj.dat',position="append")
           end if
        end if

     end if

     if ( mod(itime,wrtstep) == 0 ) then
        if ( myrank == 0 ) call mdio( 1, tote0 )
        if ( lcpmd ) then
           if ( lbathnew  ) call write_nose_data
           if ( lbathnewe ) call write_nosee_data
           call write_data_cpmdio
        end if
     end if

     call global_watch(.false.,flag_etlimit)
     if ( flag_etlimit ) then
        if ( myrank == 0 ) write(*,*) &
             "elapsed time limit exceeded : flag_etlimit=",flag_etlimit
        exit
     end if

     if ( exit_program() ) exit

  end do ! itime

!
! --- loop end
!

  if ( myrank == 0 ) call mdio( 1,tote0 )
  if ( lcpmd ) then
     if ( lbathnew  ) call write_nose_data
     if ( lbathnewe ) call write_nosee_data
     call write_data_cpmdio
  end if

  if ( lcpmd ) call dealloc_cpmd
  call dealloc_md

  if ( myrank == 0 ) then
     if ( all_traj > 0 ) close(3)
     close(4)  ! info.dat
     close(15) ! time.dat
     if ( ltime ) then
        close(16) ! time_cpmd_loop.dat
        close(17) ! time_force_once.dat
     end if
  end if

98 continue

  call write_border( 0, " CPMD END -----------" )
  call write_border( 0, "" )

  return

99 stop "stop@bomd(99)"

10 format(10f15.8)

END SUBROUTINE bomd
