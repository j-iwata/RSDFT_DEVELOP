!-----------------------------------------------------------------------
!     allocation and deallocation
!-----------------------------------------------------------------------

subroutine alloc_md
  use cpmd_variables, only: Rion,Velocity,Force,MI
  implicit none 
  if ( .not.allocated(Rion) ) allocate( Rion(3,MI) )
  Rion=0.0d0
  allocate( Velocity(3,MI) ); Velocity=0.0d0
  allocate(    Force(3,MI) ); Force=0.0d0
  return
end subroutine alloc_md



subroutine dealloc_md
  use cpmd_variables
  implicit none 
  deallocate(Rion,Velocity,Force)
  return
end subroutine dealloc_md



subroutine alloc_cpmd
  use array_bound_module, only: MB,MBZ_0,MBZ_1,MSP_0,MSP_1,ML_0,ML_1
  use cpmd_variables, only: tau,sig,gam,gamn,scr,wrk,psi_n,psi_v,MBC &
                           ,nprocs &
                           ,MB_0_CPMD,MB_1_CPMD
  implicit none 
  ! integer,allocatable :: ir_i(:),id_i(:)
  ! integer ls,le,li,i,k,n
  logical :: lblas

  lblas=.true.

  if ( lblas ) then
    allocate( tau(MBC,MBC), sig(MBC,MBC) )
    allocate( gam(MBC,MBC),gamn(MBC,MBC) )
    allocate( scr(MBC,MBC), wrk(MBC,MBC) )
  else
    ! allocate( id_i(0:nprocs-1),ir_i(0:nprocs-1) )
    ! ir_i(:)=0
    ! do i=1,MBC
    !   k=mod(i-1,nprocs)
    !   ir_i(k)=ir_i(k)+1
    ! end do
    ! do k=0,nprocs-1
    !   id_i(k)=sum(ir_i(0:k))-ir_i(k)
    ! end do
    ! ls=id_i(myrank)+1
    ! le=id_i(myrank)+ir_i(myrank)
    ! li=ir_i(myrank)
    ! deallocate( id_i,ir_i )
    ! allocate( tau(MBC,ls:le), sig(MBC,ls:le) )
    ! allocate( gam(MBC,ls:le),gamn(MBC,ls:le) )
    ! allocate( scr(MBC,ls:le), wrk(MBC,MBC)   )
  endif
  tau=0.0d0 ; sig =0.0d0
  gam=0.0d0 ; gamn=0.0d0
  scr=0.0d0 ; wrk =0.0d0
!  allocate( psi_v(ML_0:ML_1,MB,MBZ_0:MBZ_1,MSP_0:MSP_1) )
!  allocate( psi_n(ML_0:ML_1,MB,MBZ_0:MBZ_1,MSP_0:MSP_1) )
  allocate( psi_v(ML_0:ML_1,MB_0_CPMD:MB_1_CPMD,MBZ_0:MBZ_1,MSP_0:MSP_1) )
  allocate( psi_n(ML_0:ML_1,MB_0_CPMD:MB_1_CPMD,MBZ_0:MBZ_1,MSP_0:MSP_1) )
  psi_v=0.0d0
  psi_n=0.0d0
  return
end subroutine alloc_cpmd



subroutine dealloc_cpmd
  use cpmd_variables
  implicit none 
  integer,allocatable :: ir_i(:),id_i(:)
  integer ls,le,li,i,k,n,n1,n2
  logical :: lblas
  lblas=.false.
  deallocate(tau,sig,gam,gamn,scr,wrk)
  deallocate(psi_v,psi_n,mstocck)
  return
end subroutine dealloc_cpmd



!-----------------------------------------------------------------------
!     read and send cpmd input variables
!-----------------------------------------------------------------------
subroutine read_cpmd_variables
   use cpmd_variables
   use io_tools_module
   implicit none 
   call write_border( 0, 'read_cpmd_variables(start)' )
   if ( myrank == 0 ) then
      open(2,file='cpmd_var.dat',status='old')
      read(2,*) nstep
      read(2,*) deltat
      read(2,*) temp
      read(2,*) omegan
      read(2,*) emass
      read(2,*) trange
      read(2,*) ekinw, ekin1, ekin2
      read(2,*) wnose0, ekinw
      read(2,*) lcpmd
      read(2,*) lbath
      read(2,*) inivel
      read(2,*) lmeta
      read(2,*) lbere
      read(2,*) lbathnew
      read(2,*) lscale
      read(2,*) lscaleele
      read(2,*) lquench
      read(2,*) lforce_fast
      read(2,*) lbathnewe
      read(2,*) linitnose
      read(2,*) linitnosee
      read(2,*) lblue
   end if
   call IOTools_readIntegerKeyword( "TRJSTEP", trjstep , 2 )
   call IOTools_readIntegerKeyword( "WRTSTEP", wrtstep , 2 )
   call IOTools_readIntegerKeyword( "ALLTRAJ", all_traj, 2 )
   call IOTools_readIntegerKeyword( "CPMDIO"   , ctrl_cpmdio, 2 )
   call IOTools_readIntegerKeyword( "CPMDWRITE", ctrl_cpmdio, 2 )
   if ( myrank == 0 ) then
      close(2)
      write(*,*) "nstep =",nstep
      write(*,*) "deltat=",deltat
      write(*,*) "temp  =",temp
      write(*,*) "omegan=",omegan
      write(*,*) "emass =",emass
      write(*,*) "trange=",trange
      write(*,*) "ekin1,ekin2=",ekin1,ekin2
      write(*,*) "wnose0=",wnose0
      write(*,*) "ekinw =",ekinw
      write(*,*) "lcpmd =",lcpmd
      write(*,*) "lbath =",lbath
      write(*,*) "inivel=",inivel
      write(*,*) "lmeta =",lmeta
      write(*,*) "lbere =",lbere
      write(*,*) "lbathnew=",lbathnew
      write(*,*) "lscale=",lscale
      write(*,*) "lscalee=",lscaleele
      write(*,*) "lquench=",lquench
      write(*,*) "lforce_fast",lforce_fast
      write(*,*) "linitnose=",linitnose
      write(*,*) "linitnosee=",linitnosee
      write(*,*) "lblue=",lblue
      if ( lblue ) write(*,*) "Constraint ON"
      write(*,*) "TRJSTEP=",trjstep
      write(*,*) "WRTSTEP=",wrtstep
      write(*,*) "all_traj=",all_traj
   end if
   call send_cpmd_variables
   call write_border( 0, 'read_cpmd_variables(end)' )
   return
end subroutine read_cpmd_variables



subroutine send_cpmd_variables
!   use global_variables
   use cpmd_variables
   implicit none 
   integer :: ierr
   call write_border( 0, 'send_cpmd_variables(start)' )
   call mpi_bcast(nstep,1,mpi_integer,0,mpi_comm_world,ierr)
   call mpi_bcast(deltat,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(trange,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(ekinw,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(ekin1,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(ekin2,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(temp,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(omegan,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(emass,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(wnose0,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(ekinw,1,mpi_real8,0,mpi_comm_world,ierr)
   call mpi_bcast(lcpmd,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lbath,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(inivel,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lmeta,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lbere,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lbathnew,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lbathnewe,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lscale,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lscaleele,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lquench,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lforce_fast,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lrestart_e,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lrestart_p,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(linitnose,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(linitnosee,1,mpi_logical,0,mpi_comm_world,ierr)
   call mpi_bcast(lblue,1,mpi_logical,0,mpi_comm_world,ierr)
   call write_border( 0, 'send_cpmd_variables(end)' )
   return
end subroutine send_cpmd_variables
