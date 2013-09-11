!-----------------------------------------------------------------------
!     Evaluate wave function force
!-----------------------------------------------------------------------
SUBROUTINE wf_force
  use rgrid_module, only: dV
  use array_bound_module
  use wf_module
  use parallel_module, only: MB_d
  use cpmd_variables, only: psi_n,mstocck,emass,MB_0_CPMD,MB_1_CPMD
  use watch_module
  use hamiltonian_module
  implicit none
  include 'mpif.h'
  integer :: i,s,k,ns,ne,MBT,ierr
  real(8) :: c
  real(8),allocatable :: esp0(:,:,:)
  real(8) :: ctime0,ctime1,etime0,etime1

   call watch(ctime0,etime0)

   allocate( esp0(MB,MBZ,MSP) ) ; esp0=0.d0

   do s=MSP_0,MSP_1
   do k=MBZ_0,MBZ_1
      MBT=mstocck(k,s)
      do ns=MB_0_CPMD,MB_1_CPMD,MB_d
         ne=min(ns+MB_d-1,MB_1_CPMD)
#ifdef _DRSDFT_
         call hamiltonian &
              (k,s,unk(ML_0,ns,k,s),psi_n(ML_0,ns,k,s),ML_0,ML_1,ns,ne)
#endif
         esp0(ns,k,s)=sum( unk(:,ns,k,s)*psi_n(:,ns,k,s) )*dV
      enddo  ! band-loop
   enddo  ! k-loop
   enddo  ! s-loop

   call mpi_allreduce &
        (esp0,esp,MB*MBZ*MSP,mpi_real8,mpi_sum,mpi_comm_world,ierr)

   deallocate( esp0 )

   c=-1.d0/emass
   psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:) = c*psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)

   call watch(ctime1,etime1)

   return
END SUBROUTINE wf_force
