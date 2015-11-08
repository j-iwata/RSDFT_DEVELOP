MODULE ps_getDij_module

  use parallel_module, only: myrank,nprocs_g ! no need of nprocs_g
  use array_bound_module, only: MSP_0,MSP_1
  use VarPSMemberG, only: Dij00,Dij,N_nzqr,QRij
  use VarParaPSnonLocG, only: MJJ_Q,JJP_Q,nrqr_xyz,num_2_rank_Q &
                             ,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q &
                             ,rbufnl_Q,nl_max_send_Q
  use para_rgrid_comm, only: do3StepComm_real
  use localpot_module, only: Vloc
  use rgrid_module, only: dV
  use pseudopot_module, only: pselect

  implicit none

  PRIVATE
  PUBLIC :: getDij

CONTAINS

  SUBROUTINE getDij

    implicit none
    integer :: s,kk1,i,j
    real(8) :: IntQV(N_nzqr,MSP_0:MSP_1)
    real(8) :: pIntQV(N_nzqr,MSP_0:MSP_1)

    select case ( pselect )
    case default

       return

    case ( 102 )

       pIntQV = 0.0d0
       IntQV  = 0.0d0

       do s=MSP_0,MSP_1
          do kk1=1,N_nzqr
             do j=1,MJJ_Q(kk1)
                i=JJP_Q(j,kk1)
                pIntQV(kk1,s)=pIntQV(kk1,s)+QRij(j,kk1)*Vloc(i,s)
             end do
             pIntQV(kk1,s)=dV*pIntQV(kk1,s)
          end do
       end do

       call do3StepComm_real( nrqr_xyz,num_2_rank_Q, &
                            & sendmap_Q,recvmap_Q,qr_nsend, &
                            & N_nzqr,MSP_0,MSP_1,pIntQV )
       IntQV(:,:)=pIntQV(:,:)

       do s=MSP_0,MSP_1
          Dij(1:N_nzqr,s) = Dij00(1:N_nzqr) + IntQV(1:N_nzqr,s)
       end do

    end select

  END SUBROUTINE getDij

END MODULE ps_getDij_module
