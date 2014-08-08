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
write(400+myrank,*) ">>>> getDij"
    pIntQV  =0.d0
    IntQV   =0.d0

    select case ( pselect )
    case ( 2 )
      stop 'not implemented yet'
    case ( 102 )
write(501,*) 'getDij'
      do s=MSP_0,MSP_1
        do kk1=1,N_nzqr
write(501,*) 's,kk1,MJJ_Q(kk1)'
write(501,'(3I5)') s,kk1,MJJ_Q(kk1)
write(501,*) 's,kk1,i,j,QRij(j,kk1),Vloc(i,s)'
          do j=1,MJJ_Q(kk1)
            i=JJP_Q(j,kk1)
            pIntQV(kk1,s)=pIntQV(kk1,s)+QRij(j,kk1)*Vloc(i,s)
write(501,'(4I5,2g20.7)') s,kk1,i,j,QRij(j,kk1),Vloc(i,s)
          end do
          pIntQV(kk1,s)=dV*pIntQV(kk1,s)
write(501,*) 's,kk1,pIntQV(kk1,s)'
write(501,'(2I5,2g20.7)') s,kk1,pIntQV(kk1,s)
        end do
      end do
      call threeWayComm( nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,N_nzqr,MSP_0,MSP_1,pIntQV )
      IntQV(:,:)=real(pIntQV(:,:))
      do s=MSP_0,MSP_1
        Dij(1:N_nzqr,s)=Dij00(1:N_nzqr)+IntQV(1:N_nzqr,s)
do i=1,N_nzqr
write(540,'(2I5,g20.7)') s,i,Dij(i,s)
enddo
      end do
    case default
      stop 'pselect\=2,102 is not availiable'
    end select
write(400+myrank,*) "<<<< getDij"
    return
  END SUBROUTINE getDij

END MODULE PSnonLocDij
