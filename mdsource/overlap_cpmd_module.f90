MODULE overlap_cpmd_module

  use wf_module, only: unk
  use rgrid_module, only: dV
  use parallel_module
  use cpmd_variables, only: wrk,tau,sig,scr,MBC,MBT,psi_n,psi_v &
       ,ir_band_cpmd,id_band_cpmd,MB_0_CPMD,MB_1_CPMD
  use watch_module
  use overlap_module

  implicit none

  PRIVATE
  PUBLIC :: overlap2, overlap4, overlap5

CONTAINS


  SUBROUTINE overlap2(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr
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

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)

    allocate( ir(0:np_band-1) ) ; ir=0
    allocate( id(0:np_band-1) ) ; id=0
    ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
    id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0

    call watcht(myrank==0,"",0)

    call mpi_allgatherv(psi_n(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
         ,psi_n(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)

    call watcht(myrank==0,"overlap(2-1)",1)

    call dsyrk('u','t',MBT,ML0,-dV,psi_n(n1,1,k,s),ML0,zero,sig,MBC)

    call watcht(myrank==0,"overlap(2-2)",1)

    call mpi_allreduce &
         (MPI_IN_PLACE,sig,MBC*MBC,mpi_real8,mpi_sum,comm_grid,ierr)

    call watcht(myrank==0,"overlap(2-3)",1)

!$OMP parallel do
    do i=1,MBC
    do j=1,i-1
       sig(i,j) = sig(j,i)
    end do
    end do
!$OMP end parallel do

    call watcht(myrank==0,"overlap(2-4)",1)

    deallocate( id,ir )

    return

  END SUBROUTINE overlap2


  SUBROUTINE overlap4(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr
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

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)

    allocate( ir(0:np_band-1) ) ; ir=0
    allocate( id(0:np_band-1) ) ; id=0
    ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
    id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0

    call watcht(myrank==0,"",0)

    call mpi_allgatherv(  unk(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
         ,  unk(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
    call mpi_allgatherv(psi_n(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
         ,psi_n(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)

    call watcht(myrank==0,"overlap(4-1)",1)

#ifdef _DRSDFT_
    call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_n(n1,1,k,s),-dV,tau)
#endif

    call watcht(myrank==0,"overlap(4-2)",1)

!$OMP parallel do
    do i=1,MBC
    do j=1,i-1
       tau(j,i) = tau(i,j)
    end do
    end do
!$OMP end parallel do

    call watcht(myrank==0,"overlap(4-3)",1)

    deallocate( id,ir )

    return

  END SUBROUTINE overlap4


  SUBROUTINE overlap5(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr
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

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)

    allocate( ir(0:np_band-1) ) ; ir=0
    allocate( id(0:np_band-1) ) ; id=0
    ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
    id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0

    call watcht(myrank==0,"",0)

    call mpi_allgatherv(  unk(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
         ,  unk(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
    call mpi_allgatherv(psi_v(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
         ,psi_v(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)

    call watcht(myrank==0,"overlap(5-1)",1)

#ifdef _DRSDFT_
    call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_v(n1,1,k,s),dV,wrk)
#endif

    call watcht(myrank==0,"overlap(5-2)",1)

!$OMP parallel do
    do i=1,MBC
    do j=1,i-1
       wrk(j,i) = wrk(i,j)
    end do
    end do
!$OMP end parallel do 

    call watcht(myrank==0,"overlap(5-3)",1)

    deallocate( id,ir )

    return

  END SUBROUTINE overlap5


END MODULE overlap_cpmd_module
