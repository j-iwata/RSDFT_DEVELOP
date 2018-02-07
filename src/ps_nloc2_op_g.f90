MODULE PSnonLocOpG2

  use parallel_module, only: myrank
  use var_ps_member
  use var_ps_member_g
  use ps_nloc2_variables, only: Mlma,nzlma,MJJ,JJP,nrlma_xyz &
       ,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,zero,uVk,iuV
  use rgrid_module, only: dV, Igrid
  use para_rgrid_comm, only: do3StepComm

  implicit none

  include 'mpif.h'

  private
  public :: op_ps_nloc2_uspp

CONTAINS

  SUBROUTINE op_ps_nloc2_uspp( tpsi, htpsi, htpsi00, k_in, s_in )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8),intent(INOUT),optional :: htpsi00(:,:)
    real(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#else
    complex(8),intent(IN)  :: tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8),intent(INOUT),optional :: htpsi00(:,:)
    complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
#endif
    integer,optional,intent(IN) :: k_in, s_in
    integer :: i,ib,j,i1,i2,m,lma,nb,ierr,nreq,k,s,i0
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    complex(8) :: zc

    integer :: iqr,lma1,lma2

    i0 = Igrid(1,0)
    nb = size( tpsi, 2 )
    k  = 1 ; if ( present(k_in) ) k=k_in
    s  = 1 ; if ( present(s_in) ) s=s_in

    if ( Mlma <= 0 ) return

    allocate( uVunk(nzlma,nb),uVunk0(nzlma,nb) )

!!$OMP parallel

    do ib=1,nb
!!$OMP do
       do lma=1,nzlma
          uVunk(lma,ib)=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)-i0+1
#ifdef _DRSDFT_
             uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*tpsi(i,ib)
#else
             uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*tpsi(i,ib)
#endif
          end do
          uVunk(lma,ib)=dV*uVunk(lma,ib)
       end do
!!$OMP end do
    end do

    call do3StepComm( nrlma_xyz,num_2_rank,sendmap,recvmap &
                     ,lma_nsend,sbufnl,rbufnl,nzlma,1,nb,uVunk )

    if ( present(htpsi00) ) then

       do ib=1,nb
          do iqr=1,N_nzqr
             lma1=nzqr_pair(iqr,1)
             lma2=nzqr_pair(iqr,2)
             if ( lma1 < lma2 ) stop 'NZQR_PAIR is strange'
             if ( lma1 == lma2 ) then
                do j=1,MJJ(lma1)
                   i=JJP(j,lma1)-i0+1
                   htpsi00(i,ib) = htpsi00(i,ib) &
                        + Dij00(iqr)*uVk(j,lma1,k)*uVunk(lma2,ib)
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
                end do
             else
                do j=1,MJJ(lma1)
                   i=JJP(j,lma1)-i0+1
                   htpsi00(i,ib) = htpsi00(i,ib) &
                        + Dij00(iqr)*uVk(j,lma1,k)*uVunk(lma2,ib)
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
                end do
                do j=1,MJJ(lma2)
                   i=JJP(j,lma2)-i0+1
                   htpsi00(i,ib) = htpsi00(i,ib) &
                        + Dij00(iqr)*uVk(j,lma2,k)*uVunk(lma1,ib)
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma2,k)*uVunk(lma1,ib)
                end do
             end if
          end do
       end do

    else

       do ib=1,nb
          do iqr=1,N_nzqr
             lma1=nzqr_pair(iqr,1)
             lma2=nzqr_pair(iqr,2)
             if ( lma1 < lma2 ) stop 'NZQR_PAIR is strange'
             if ( lma1 == lma2 ) then
                do j=1,MJJ(lma1)
                   i=JJP(j,lma1)-i0+1
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
                end do
             else
                do j=1,MJJ(lma1)
                   i=JJP(j,lma1)-i0+1
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma1,k)*uVunk(lma2,ib)
                end do
                do j=1,MJJ(lma2)
                   i=JJP(j,lma2)-i0+1
                   htpsi(i,ib) = htpsi(i,ib) &
                        + Dij(iqr,s)*uVk(j,lma2,k)*uVunk(lma1,ib)
                end do
             end if
          end do
       end do

    end if
!!$OMP end parallel

    deallocate( uVunk0,uVunk )

  END SUBROUTINE op_ps_nloc2_uspp

END MODULE PSnonLocOpG2
