MODULE PSnonLocDij
  use parallel_module, only: myrank
  use array_bound_module, only: MSP_0,MSP_1
  use VarPSMemberG, only: Dij00,Dij,N_nzqr,QRij
  use VarParaPSnonLocG, only: MJJ_Q,JJP_Q,nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q
  use ParaRGridComm, only: threeWayComm
  use localpot_module, only: Vloc
  use rgrid_module, only: dV
  use pseudopot_module, only: pselect
  implicit none
CONTAINS
  SUBROUTINE getDij
    implicit none
    integer :: s,kk1,i,j
    real(8) :: IntQV(N_nzqr,MSP_0:MSP_1)
#ifdef _DRSDFT_
    real(8) :: pIntQV(N_nzqr,MSP_0:MSP_1),IntQV(N_nzqr,MSP_0:MSP_1)
#else
    complex(8) :: pIntQV(N_nzqr,MSP_0:MSP_1)
#endif
    pIntQV  =0.d0
    IntQV   =0.d0

    select case ( pselect )
    case ( 2 )
      stop 'not implemented yet'
    case ( 102 )
      do s=MSP_0,MSP_1
        do kk1=1,N_nzqr
          do j=1,MJJ_Q(kk1)
            i=JJP_Q(j,kk1)
            pIntQV(kk1,s)=pIntQV(kk1,s)+QRij(j,kk1)*Vloc(i,s)
          end do
          pIntQV(kk1,s)=dV*pIntQV(kk1,s)
        end do
      end do
      call threeWayComm( nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,N_nzqr,MSP_0,MSP_1,pIntQV,0 )
      IntQV(:,:)=real(pIntQV(:,:))
      do s=MSP_0,MSP_1
        Dij(1:N_nzqr,s)=Dij00(1:N_nzqr)+IntQV(1:N_nzqr,s)
      end do
    case default
      stop 'pselect\=2,102 is not availiable'
    end select
    return
  END SUBROUTINE getDij

END MODULE PSnonLocDij
