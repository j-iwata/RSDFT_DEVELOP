MODULE overlap_cpmd_mbp_module

  use wf_module, only: unk
  use rgrid_module, only: dV
! use parallel_module
  use parallel_module, only: idisp, ircnt,id_class  &
                            ,comm_grid, comm_band, myrank, myrank_g, myrank_b, np_band
  use cpmd_variables, only: wrk,tau,sig,scr,MBC,MBT,psi_n,psi_v &
                           ,ir_band_cpmd,id_band_cpmd,MB_0_CPMD,MB_1_CPMD
  use watch_module
  use calc_overlap_mbp_module, only: calc_overlap_mbp
  use rsdft_mpi_module, only: rsdft_allgatherv, rsdft_allreduce_sum

  implicit none

  PRIVATE
  PUBLIC :: overlap2_mbp, overlap4_mbp, overlap5_mbp

  real(8),allocatable :: w1(:),w2(:)

CONTAINS


  SUBROUTINE overlap2_mbp(s,k,ii0,jj0)
    implicit none
    integer,intent(IN) :: s,k
    integer,intent(OUT) :: ii0,jj0
    integer :: i,j,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir_loc(:),id_loc(:)
    integer :: ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    integer :: MB0
    logical :: iflag_diag_add

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    call watcht(myrank==0,"",0)

    call watcht(myrank==0,"overlap(2-1)",1)

    iflag_diag_add=.true.
    MB0=size(unk,2)
!   call calc_overlap_mbp(ML0,MBT,MB0,psi_n(n1,1,k,s),psi_n(n1,1,k,s),-dV,sig(:,:),iflag_diag_add,ii0,jj0)
    call calc_overlap_mbp(ML0,MBT,MB0,psi_n(n1,MB_0_CPMD,k,s),psi_n(n1,MB_0_CPMD,k,s),-dV,sig(:,:),iflag_diag_add,ii0,jj0)

!   call dsyrk('L','t',MBT,ML0,-dV,psi_n(n1,1,k,s),ML0,zero,sig,MBC)

    call watcht(myrank==0,"overlap(2-2)",1)

!   n = (MBC*(MBC+1))/2
!   if ( .not.allocated(w1) ) then
!      allocate( w1(n) ) ; w1=0.0d0
!   end if

!   do j=1,MBC
!   do i=j,MBC
!      m=(j-1)*MBC-(j*(j-1))/2+i
!      w1(m)=sig(i,j)
!   end do
!   end do

!   !call mpi_allreduce(MPI_IN_PLACE,w1,n,mpi_real8,mpi_sum,comm_grid,ierr)
!   call rsdft_allreduce_sum( w1(1:n), comm_grid )

!   do j=1,MBC
!   do i=j,MBC
!      m=(j-1)*MBC-(j*(j-1))/2+i
!      sig(i,j)=w1(m)
!   end do
!   end do

    call watcht(myrank==0,"overlap(2-3)",1)

!!$OMP parallel do
!   do j=1  ,MBC
!      !do i=j+1,MBC
!      !   sig(j,i) = sig(i,j)
!      do i=1,j-1
!         sig(i,j) = sig(j,i)
!      end do
!      sig(j,j)=sig(j,j)+1.0d0
!   end do
!!$OMP end parallel do

    call watcht(myrank==0,"overlap(2-4)",1)

    return

  END SUBROUTINE overlap2_mbp


  SUBROUTINE overlap4_mbp(s,k,ii0,jj0)
    implicit none
    integer,intent(IN) :: s,k
    integer,intent(OUT) :: ii0,jj0
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir_loc(:),id_loc(:)
    integer ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    integer :: MB0
    logical :: iflag_diag_add

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    call watcht(myrank==0,"",0)

    call watcht(myrank==0,"overlap(4-1)",1)


#ifdef _DRSDFT_
!   call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_n(n1,1,k,s),-dV,tau)
    iflag_diag_add=.true.
    MB0=size(unk,2)
!   call calc_overlap_mbp(ML0,MBT,MB0,unk(n1,1,k,s),psi_n(n1,1,k,s),-dV,tau(:,:),iflag_diag_add,ii0,jj0)
    call calc_overlap_mbp(ML0,MBT,MB0,unk(n1,MB_0_CPMD,k,s),psi_n(n1,MB_0_CPMD,k,s),-dV,tau(:,:),iflag_diag_add,ii0,jj0)
#endif

    call watcht(myrank==0,"overlap(4-2)",1)

!!$OMP parallel do
!   do j=1,MBC
!      !do i=j+1,MBC
!      !   tau(j,i) = tau(i,j)
!      do i=1,j-1
!         tau(i,j) = tau(j,i)
!      end do
!      tau(j,j)=tau(j,j)+1.0d0
!   end do
!!$OMP end parallel do

    call watcht(myrank==0,"overlap(4-3)",1)


    return

  END SUBROUTINE overlap4_mbp


  SUBROUTINE overlap5_mbp(s,k,ii0,jj0)
    implicit none
    integer,intent(IN) :: s,k
    integer,intent(OUT) :: ii0,jj0
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
    integer,allocatable :: ir(:),id(:)
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir_loc(:),id_loc(:)
    integer ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    integer :: MB0
    logical :: iflag_diag_add

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    call watcht(myrank==0,"",0)

    call watcht(myrank==0,"overlap(5-1)",1)

#ifdef _DRSDFT_
!   call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_v(n1,1,k,s),dV,wrk)
    iflag_diag_add=.false.
    MB0=size(unk,2)
!   call calc_overlap_mbp(ML0,MBT,MB0,unk(n1,1,k,s),psi_v(n1,1,k,s),dV,wrk(:,:),iflag_diag_add,ii0,jj0)
    call calc_overlap_mbp(ML0,MBT,MB0,unk(n1,MB_0_CPMD,k,s),psi_v(n1,MB_0_CPMD,k,s),dV,wrk(:,:),iflag_diag_add,ii0,jj0)
#endif

    call watcht(myrank==0,"overlap(5-2)",1)

!!$OMP parallel do
!   do j=1,MBC
!   !do i=j+1,MBC
!   !   wrk(j,i) = wrk(i,j)
!   do i=1,j-1
!      wrk(i,j) = wrk(j,i)
!   end do
!   end do
!!$OMP end parallel do 

    call watcht(myrank==0,"overlap(5-3)",1)

    return

  END SUBROUTINE overlap5_mbp


END MODULE overlap_cpmd_mbp_module
