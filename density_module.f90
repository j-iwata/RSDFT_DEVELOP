MODULE density_module

  use rgrid_module, only: dV
  use electron_module, only: Nelectron
  use wf_module
  use parallel_module,only: comm_grid,comm_band,comm_bzsm,comm_spin
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1,MSP

!--------------------------------------------------------------------------- 20140507 HL
#ifdef _USPP_
  use WFDensityG, only: get_rhonks
  use VarSysParameter, only: pp_kind
#endif
!=========================================================================== 20140507 HL

  implicit none

  PRIVATE
  PUBLIC :: rho,init_density,normalize_density,calc_density

  real(8),allocatable :: rho(:,:)

CONTAINS

!---------------------------------------------------------------------------------------
  SUBROUTINE init_density
    implicit none
    real(8) :: c,d
    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,MSP) )
       call random_number(rho)
       call normalize_density
    end if
  END SUBROUTINE init_density

!---------------------------------------------------------------------------------------
  SUBROUTINE normalize_density
    implicit none
    real(8) :: c,d
    integer :: ierr
    include 'mpif.h'
    c=sum(rho)*dV
    call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    c=Nelectron/d
    rho=c*rho
  END SUBROUTINE normalize_density

!---------------------------------------------------------------------------------------
!--------------------------------------------------------------------------- 20140507 HL
#ifdef _USPP_
  SUBROUTINE calc_density
    ! IN:	unk(:.n.k.s),ML_0,ML_1,MSP_0,MSP_1
    ! OUT:	rho(n1:n2,s)
    implicit none
    integer :: n,k,s
    integer :: n1,n2,n0
    real(8),allocatable :: rhonks(:)

    select case ( pp_kind )
    case ( 'NCPP' )
       rho(:,:)=0.d0
       do s=MSP_0,MSP_1
          do k=MBZ_0,MBZ_1
             do n=MB_0 ,MB_1
                rho(:,s)=rho(:,s)+occ(n,k,s)*abs( unk(:,n,k,s) )**2
             end do
          end do
       end do
    case ( 'USPP' )
       n1=ML_0
       n2=ML_1
       n0 = ML_1 - ML_0 + 1

       allocate( rhonks(n1:n2) ) ; rhonks(:)=0.d0

       rho(:,:)=0.d0
       do s=MSP_0,MSP_1
          do k=MBZ_0,MBZ_1
             do n=MB_0,MB_1
                rhonks(:)=0.d0
                call get_rhonks( rhonks,n1,n2,n,k,s )
                rho(:,s) = rho(:,s) + occ(n,k,s)*rhonks(:)
             end do
          end do
       end do
    end select

    call reduce_and_gather

    deallocate( rhonks )

  END SUBROUTINE calc_density
#else
  SUBROUTINE calc_density
    implicit none
    integer :: n,k,s
    rho(:,:)=0.d0
    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do n=MB_0 ,MB_1
             rho(:,s)=rho(:,s)+occ(n,k,s)*abs( unk(:,n,k,s) )**2
          end do
       end do
    end do
    call reduce_and_gather
  END SUBROUTINE calc_density
#endif
!=========================================================================== 20140507 HL

!---------------------------------------------------------------------------------------
  SUBROUTINE reduce_and_gather
    implicit none
    real(8),allocatable :: w(:)
    integer :: n,k,s,m,ierr
    include 'mpif.h'
    m=ML_1-ML_0+1
    allocate( w(m) )
    do s=MSP_0,MSP_1
       call mpi_allreduce(rho(ML_0,s),w,m,mpi_real8,mpi_sum,comm_band,ierr)
       call mpi_allreduce(w,rho(ML_0,s),m,mpi_real8,mpi_sum,comm_bzsm,ierr)
    end do
    ! The following assumes all 'MSP_1-MSP_0+1' are the same
    m=m*(MSP_1-MSP_0+1)
    call mpi_allgather(rho(ML_0,MSP_0),m,mpi_real8,rho,m,mpi_real8 &
         ,comm_spin,ierr)
    deallocate( w )
  END SUBROUTINE reduce_and_gather

END MODULE density_module
