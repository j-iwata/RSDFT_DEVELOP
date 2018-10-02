!-----------------------------------------------------------------------
!     Evaluate wave function force
!-----------------------------------------------------------------------
SUBROUTINE wf_force
  use rgrid_module, only: dV
  use array_bound_module
  use wf_module
  use parallel_module, only: MB_d, myrank
  use cpmd_variables, only: psi_n,mstocck,emass,MB_0_CPMD,MB_1_CPMD
  use watch_module
  use hamiltonian_module
  implicit none
  include 'mpif.h'
  integer :: i,s,k,ns,ne,n,MBT,ierr
  real(8) :: c
  real(8) :: ctime0,ctime1,etime0,etime1
  real(8) :: ct(0:7),et(0:7),ct0,ct1,ct2,et0,et1,et2

  ct=0.0d0
  et=0.0d0

   call watch(ct(0),et(0))

   allocate( esp0(MB,MBZ,MSP) ) ; esp0=0.d0

   call watch(ct(1),et(1))

   do s=MSP_0,MSP_1
   do k=MBZ_0,MBZ_1
      MBT=mstocck(k,s)
      do ns=MB_0_CPMD,MB_1_CPMD,MB_d
         ne=min(ns+MB_d-1,MB_1_CPMD)
         call watch(ct0,et0)
#ifdef _DRSDFT_
         call hamiltonian &
              (k,s,unk(:,ns:ne,k,s),psi_n(:,ns:ne,k,s),ML_0,ML_1,ns,ne)
#endif
         call watch(ct1,et1)
         do n=ns,ne
            esp0(n,k,s)=sum( unk(:,n,k,s)*psi_n(:,n,k,s) )*dV
         end do
         call watch(ct2,et2)
         ct(2)=ct(2)+ct1-ct0 ; et(2)=et(2)+et1-et0
         ct(3)=ct(3)+ct2-ct1 ; et(3)=et(3)+et2-et1
      enddo  ! band-loop
   enddo  ! k-loop
   enddo  ! s-loop

   call watch(ct(4),et(4))

   call mpi_allreduce &
        (esp0,esp,MB*MBZ*MSP,mpi_real8,mpi_sum,mpi_comm_world,ierr)

   call watch(ct(5),et(5))

   deallocate( esp0 )

   call watch(ct(6),et(6))

   c=-1.d0/emass
   psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:) = c*psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)

   call watch(ct(7),et(7))

!   if ( myrank == 0 ) then
!      write(*,*) "time(wf_fore1)=",ct(1)-ct(0),et(1)-et(0)
!      write(*,*) "time(wf_fore2)=",ct(2),et(2)
!      write(*,*) "time(wf_fore3)=",ct(3),et(3)
!      write(*,*) "time(wf_fore5)=",ct(5)-ct(4),et(5)-et(4)
!      write(*,*) "time(wf_fore6)=",ct(6)-ct(5),et(6)-et(5)
!      write(*,*) "time(wf_fore7)=",ct(7)-ct(6),et(7)-et(6)
!   end if

   return
END SUBROUTINE wf_force
