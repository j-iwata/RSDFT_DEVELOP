!-----------------------------------------------------------------------
!     calculate fictitious kinetic energy
!-----------------------------------------------------------------------
SUBROUTINE calfke(fke)
  use cpmd_variables, only: DISP_SWITCH,psi_v,emass,MB_0_CPMD,MB_1_CPMD
  use electron_module, only: occ
  use array_bound_module
  use rgrid_module, only: dV
  use watch_module
  implicit none
  real(8),intent(INOUT) :: fke
  integer :: n,k,s,ierr
  real(8) :: ctime0,ctime1,etime0,etime1,fke0,fke1
  include 'mpif.h'

  call watch(ctime0,etime0)

  fke0=0.d0
  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1
  do n=MB_0_CPMD,MB_1_CPMD
     if ( abs(occ(n,k,s))<1.d-10 ) cycle
     fke0=fke0+sum(abs(psi_v(:,n,k,s))**2)*occ(n,k,s)
  end do
  end do
  end do

  call mpi_allreduce(fke0,fke1,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
  fke=fke1*emass*dV

  call watch(ctime1,etime1)

  if (DISP_SWITCH) then
     write(*,*) "TIME(calfke)=",ctime1-ctime0,etime1-etime0
  end if

   return
END SUBROUTINE calfke
