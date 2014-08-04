Module PSnonLocOpG2
use parallel_module, only: myrank
  use VarPSMember
  use VarPSMemberG
  use ps_nloc2_variables, only: Mlma,nzlma,MJJ,JJP,nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,zero,uVk,iuV
  use rgrid_module, only: dV
  use ParaRGridComm, only: threeWayComm
  implicit none
  include 'mpif.h'
  private
  public :: op_ps_nloc2_uspp
CONTAINS

  SUBROUTINE op_ps_nloc2_uspp(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#endif
    integer :: i,ib,j,i1,i2,m,lma,nb,ierr,nreq
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    complex(8) :: zc

    integer :: iqr,lma1,lma2
#ifdef _SHOWALL_OP_
write(200+myrank,*) ">>>>op_ps_nloc2_uspp"
write(200+myrank,*) '----------------------------op_ps_nloc2_uspp'
write(200+myrank,'(A17,6I7)') 'k,s,n1,n2,ib1,ib2',k,s,n1,n2,ib1,ib2
write(200+myrank,*) '----------------------------op_ps_nloc2_uspp'
#endif
    nb = ib2-ib1+1

    if ( Mlma <= 0 ) return

    allocate( uVunk(nzlma,ib1:ib2),uVunk0(nzlma,ib1:ib2) )

!!$OMP parallel

    do ib=ib1,ib2
!!$OMP do
       do lma=1,nzlma
          uVunk(lma,ib)=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)
#ifdef _DRSDFT_
             uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*tpsi(i,ib)
#else
             uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*tpsi(i,ib)
#endif
          end do
          uVunk(lma,ib)=iuV(lma)*dV*uVunk(lma,ib)
       end do
!!$OMP end do
    end do

!    select case( iswitch_eqdiv )
!    case default
    
      call threeWayComm( nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,uVunk )

!    case( 2 )
!       call comm_eqdiv_ps_nloc2_mol(nzlma,ib1,ib2,uVunk)
!    end select

    do ib=ib1,ib2
      do iqr=1,N_nzqr
        lma1=nzqr_pair(iqr,1)
        lma2=nzqr_pair(iqr,2)
        if (lma1<lma2) stop 'NZQR_PAIR is strange'
        if (lma1==lma2) then
          do j=1,MJJ(lma1)
            i=JJP(j,lma1)
            htpsi(i,ib)=htpsi(i,ib)+Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
          end do
        else
          do j=1,MJJ(lma1)
            i=JJP(j,lma1)
            htpsi(i,ib)=htpsi(i,ib)+Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
          end do
          do j=1,MJJ(lma2)
            i=JJP(j,lma2)
            htpsi(i,ib)=htpsi(i,ib)+Dij(iqr,s)*uVk(j,lma2,k)*uVunk(lma1,ib)
          end do
        end if
      end do
    end do
!!$OMP end parallel

    deallocate( uVunk0,uVunk )
#ifdef _SHOWALL_OP_
write(200+myrank,*) ">>>>op_ps_nloc2_uspp"
#endif


    return
  END SUBROUTINE op_ps_nloc2_uspp
End Module PSnonLocOpG2
