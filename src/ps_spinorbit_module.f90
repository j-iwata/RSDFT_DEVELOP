MODULE ps_spinorbit_module

  use ps_nloc2_variables
  use pseudopot_module, only: flag_so
  use grid_module, only: get_info_rs_grid

  implicit none

  PRIVATE
  PUBLIC :: op_ps_spinorbit

  integer :: comm_grid
  real(8) :: dV=0.0d0

CONTAINS


  SUBROUTINE op_ps_spinorbit( k, n1,n2, tpsi, htpsi )

    implicit none
    integer,intent(IN) :: k,n1,n2
    complex(8),intent(IN)  :: tpsi(n1:n2,2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,2)
    complex(8),allocatable :: uVunk(:,:),uVunk0(:,:)
    include 'mpif.h'
    integer :: i,is,j,i1,i2,m,lma,ns,ierr,nreq,lma1,lma2
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    complex(8) :: zc1,zc2
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: work(:,:)

    if ( .not.flag_so ) return

    if ( Mlma <= 0 ) return

!    call write_border( 1, " op_ps_spin_orbit(start)" )

    if ( dV == 0.0d0 ) call get_info_rs_grid( dV_rs=dV, comm_rs=comm_grid )

    ns = size( tpsi, 2 )

    allocate( uVunk(nzlma,ns)  ) ; uVunk=zero
    allocate( uVunk0(nzlma,ns) ) ; uVunk0=zero

!$OMP parallel

    do is=1,ns
!$OMP do
       do lma=1,nzlma
          uVunk(lma,is)=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             uVunk(lma,is)=uVunk(lma,is)+conjg(uVk_so00(j,lma,k))*tpsi(i,is)
          end do
          uVunk(lma,is)=uVunk(lma,is)*dV
       end do
!$OMP end do
    end do

!$OMP single
    do i=1,6
       select case(i)
       case(1,3,5)
!!$OMP single
          j=i+1
!!$OMP end single
!!$OMP workshare
          uVunk0(:,:)=uVunk(:,:)
!!$OMP end workshare
       case(2,4,6)
!!$OMP single
          j=i-1
!!$OMP end single
       end select
!!$OMP single
       do m=1,nrlma_xyz(i)
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if( irank >= 0 )then
             i2=0
             do is=1,ns
                do i1=1,lma_nsend(irank)
                   i2=i2+1
                   sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),is)
                end do
             end do
             nreq=nreq+1
             call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*ns &
                  ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
          end if
          if( jrank >= 0 )then
             nreq=nreq+1
             call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*ns &
                  ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          call mpi_waitall(nreq,ireq,istatus,ierr)
          if( jrank >= 0 )then
             i2=0
             do is=1,ns
                do i1=1,lma_nsend(jrank)
                   i2=i2+1
                   uVunk(recvmap(i1,jrank),is) &
                        =uVunk(recvmap(i1,jrank),is)+rbufnl(i2,jrank)
                end do
             end do
          end if
       end do
!!$OMP end single
    end do
!$OMP end single

!    allocate( work(size(htpsi,1),2) ) ; work=zero

    do m=1,N_nzqr
       lma1=nzqr_pair(1,m)
       lma2=nzqr_pair(2,m)
       zc1=Kij00(m)*uVunk(lma2,1)
       zc2=Kij00(m)*uVunk(lma2,2)
!$OMP do
       do j=1,MJJ(lma1)
          i=JJP(j,lma1)
          htpsi(i,1)=htpsi(i,1)+uVk_so11(j,lma1,k)*zc1+uVk_so12(j,lma1,k)*zc2
          htpsi(i,2)=htpsi(i,2)+uVk_so21(j,lma1,k)*zc1+uVk_so22(j,lma1,k)*zc2
!          work(i,1)=work(i,1)+uVk_so11(j,lma1,k)*zc1+uVk_so12(j,lma1,k)*zc2
!          work(i,2)=work(i,2)+uVk_so21(j,lma1,k)*zc1+uVk_so22(j,lma1,k)*zc2
       end do
!$OMP end do
    end do ! m

!$OMP end parallel

!    write(*,*) m,sum(conjg(tpsi(:,1))*work(:,1)+conjg(tpsi(:,2))*work(:,2))*dV
!    htpsi=htpsi+work
!    deallocate( work )

    deallocate( uVunk0,uVunk )

!    call write_border( 1, " op_ps_spin_orbit(end)" )

  END SUBROUTINE op_ps_spinorbit


END MODULE ps_spinorbit_module
