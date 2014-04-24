!-----------------------------------------------------------------------
!     orthogonalization of wave function velocity
!-----------------------------------------------------------------------
subroutine rotorb2

  use parallel_module
  use wf_module
  use electron_module, only: occ
  use array_bound_module, only: MSP_0,MSP_1,MBZ_0,MBZ_1
  use cpmd_variables, only: MBC,MBT,mstocck,psi_v,wrk,MB_0_CPMD,MB_1_CPMD
  use overlap_cpmd_module

  implicit none

  integer :: k,s,n1,n2,ML0,MB0
  real(8),parameter :: on= 1.0d0
  real(8),parameter :: om=-1.0d0

  n1  = idisp(myrank)+1
  n2  = idisp(myrank)+ircnt(myrank)
  ML0 = ircnt(myrank)
  MB0 = MB_1_CPMD-MB_0_CPMD+1

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1

     MBT = mstocck(k,s)

     call overlap5(s,k)

     call dgemm('n','n',ML0,MB0,MBT,om,unk(n1,1,k,s),ML0 &
          ,wrk(1,MB_0_CPMD),MBC,on,psi_v(n1,MB_0_CPMD,k,s),ML0)

  enddo ! k
  enddo ! s

  return

end subroutine rotorb2
