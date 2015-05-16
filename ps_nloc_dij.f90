MODULE PSnonLocDij
  use parallel_module, only: myrank,nprocs_g ! no need of nprocs_g
  use array_bound_module, only: MSP_0,MSP_1
  use VarPSMemberG, only: Dij00,Dij,N_nzqr,QRij
  use VarParaPSnonLocG, only: MJJ_Q,JJP_Q,nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,nl_max_send_Q
  use ParaRGridComm, only: do3StepComm_real
  use localpot_module, only: Vloc
  use rgrid_module, only: dV
  use pseudopot_module, only: pselect
  implicit none
CONTAINS
  SUBROUTINE getDij
    implicit none
    integer :: s,kk1,i,j
    real(8) :: IntQV(N_nzqr,MSP_0:MSP_1)
    real(8) :: pIntQV(N_nzqr,MSP_0:MSP_1)

    pIntQV  =0.d0
    IntQV   =0.d0

    select case ( pselect )
    case ( 2 )
       return
!      stop 'not implemented yet'
    case ( 102 )
!!write(8600+myrank,*) repeat('-',60)
!write(8700+myrank,*) repeat('-',60)
!write(8800+myrank,*) repeat('-',60)
      do s=MSP_0,MSP_1
        do kk1=1,N_nzqr
          do j=1,MJJ_Q(kk1)
            i=JJP_Q(j,kk1)
            pIntQV(kk1,s)=pIntQV(kk1,s)+QRij(j,kk1)*Vloc(i,s)
!write(8600+myrank,'(4I6,2G20.7)') s,kk1,j,i,QRij(j,kk1),Vloc(i,s)
          end do
          pIntQV(kk1,s)=dV*pIntQV(kk1,s)
!write(8700+myrank,'(2I6,3G20.7)') s,kk1,pIntQV(kk1,s),dV
        end do
      end do
!write(5000+myrank,*) repeat('=',60)
!write(5000+myrank,'(6I6)') nrqr_xyz(1:6)
!write(5000+myrank,'(I6)') nl_max_send_Q
!write(5000+myrank,*) repeat('-',20)
do i=0,nprocs_g-1
do j=1,nl_max_send_Q
!write(5000+myrank,'(3I6)') i,j,sendmap_Q(j,i)
enddo
enddo
!write(5000+myrank,*) repeat('-',20)
      call do3StepComm_real( nrqr_xyz,num_2_rank_Q, &
                          & sendmap_Q,recvmap_Q,qr_nsend, &
                          & N_nzqr,MSP_0,MSP_1,pIntQV )
      IntQV(:,:)=pIntQV(:,:)
do s=MSP_0,MSP_1
do kk1=1,N_nzqr
!write(8800+myrank,'(2I6,1G20.7)') s,kk1,IntQV(kk1,s)
enddo
enddo
      do s=MSP_0,MSP_1
        Dij(1:N_nzqr,s)=Dij00(1:N_nzqr)+IntQV(1:N_nzqr,s)
      end do
do s=MSP_0,MSP_1
do kk1=1,N_nzqr
!write(8820+myrank,'(2I6,1G20.7)') s,kk1,Dij(kk1,s)
enddo
enddo
    case default
      stop 'pselect\=2,102 is not availiable'
    end select
    return
  END SUBROUTINE getDij

END MODULE PSnonLocDij
