MODULE density_module

  use wf_module
  use parallel_module, only: comm_grid,comm_band,comm_bzsm,comm_spin,ir_grid,id_grid,ir_spin,id_spin,myrank_g,myrank_s,myrank
  use symmetry_module, only: sym_rho
  USE BasicTypeFactory
#ifdef _USPP_
  use WFDensityG, only: get_rhonks
  use VarSysParameter, only: pp_kind
#endif

  implicit none

  PRIVATE
  PUBLIC :: init_density
  PUBLIC :: normalize_density
  PUBLIC :: calc_density
  PUBLIC :: get_range_density

  real(8),allocatable,PUBLIC :: rho(:,:)
  real(8),PUBLIC :: sum_dspin(2)

  integer :: ML_RHO,ML_0_RHO,ML_1_RHO
  integer :: MS_RHO,MS_0_RHO,MS_1_RHO

  real(8) :: Nelectron_RHO,dV_RHO

  include 'mpif.h'

CONTAINS

!---------------------------------------------------------------------------------------
  SUBROUTINE get_range_density( g_range, s_range )
    implicit none
    type( ArrayRange1D ),intent(OUT) :: g_range,s_range
    g_range%head = ML_0_RHO
    g_range%tail = ML_1_RHO
    g_range%size = ML_1_RHO - ML_0_RHO + 1
    s_range%head = 1
    s_range%tail = MS_RHO
    s_range%size = MS_RHO
    g_range%head_global = 1
    g_range%tail_global = ML_RHO
    g_range%size_global = ML_RHO
    s_range%head_global = 1
    s_range%tail_global = MS_RHO
    s_range%size_global = MS_RHO
  END SUBROUTINE get_range_density

  SUBROUTINE init_density(Nelectron,dV)
    implicit none
    real(8),intent(IN) :: Nelectron,dV
    integer :: i,s

    Nelectron_RHO = Nelectron
    dV_RHO        = dV

    ML_RHO   = sum( ir_grid )
    ML_0_RHO = id_grid(myrank_g) + 1
    ML_1_RHO = id_grid(myrank_g) + ir_grid(myrank_g)
    MS_RHO   = sum( ir_spin )
    MS_0_RHO = id_spin(myrank_s) + 1
    MS_1_RHO = id_spin(myrank_s) + ir_spin(myrank_s)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0_RHO:ML_1_RHO,MS_RHO) )
       call random_number(rho)
       call normalize_density
    end if

#ifdef _SHOWALL_RHO_
    write(410+myrank,*) 'i,s,rho(i,s)'
    do s=MS_0_RHO,MS_1_RHO
       do i=ML_0_RHO,ML_1_RHO
          write(410+myrank,'(2I5,g20.7)') i,s,rho(i,s)
       enddo
    enddo
#endif

  END SUBROUTINE init_density

!---------------------------------------------------------------------------------------
  SUBROUTINE normalize_density
    implicit none
    real(8) :: c,d
    integer :: ierr
    include 'mpif.h'
    c=sum(rho)*dV_RHO
    call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    c=Nelectron_RHO/d
    rho=c*rho
  END SUBROUTINE normalize_density

!---------------------------------------------------------------------------------------
#ifdef _USPP_

  SUBROUTINE calc_density
    ! IN:	unk(:.n.k.s),ML_0,ML_1,MS_0,MS_1
    ! OUT:	rho(n1:n2,s)
    implicit none
    integer :: n,k,s,i
    integer :: n1,n2,n0
    real(8),allocatable :: rhonks(:)

    select case ( pp_kind )
    case ( 'NCPP' )

       rho(:,:)=0.0d0
       do s=MS_0_WF,MS_1_WF
          do k=MK_0_WF,MK_1_WF
             do n=MB_0_WF ,MB_1_WF
                rho(:,s)=rho(:,s)+occ(n,k,s)*abs( unk(:,n,k,s) )**2
             end do
          end do
       end do

    case ( 'USPP' )

       n1 = ML_0_WF
       n2 = ML_1_WF
       n0 = ML_1_WF - ML_0_WF + 1

       allocate( rhonks(n1:n2) ) ; rhonks(:)=0.0d0

       rho(:,:)=0.0d0
       do s=MS_0_WF,MS_1_WF
          do k=MK_0_WF,MK_1_WF
             do n=MB_0_WF,MB_1_WF
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
    rho(:,:)=0.0d0
    do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
          do n=MB_0_WF,MB_1_WF
             rho(:,s)=rho(:,s)+occ(n,k,s)*abs( unk(:,n,k,s) )**2
          end do
       end do
    end do
    call sym_rho( ML_0_RHO, ML_1_RHO, MS_RHO, MS_0_RHO, MS_1_RHO, rho )
    call reduce_and_gather
    call calc_sum_dspin
  END SUBROUTINE calc_density

#endif

!---------------------------------------------------------------------------------------
  SUBROUTINE reduce_and_gather
    implicit none
    integer :: n,k,s,m,ierr
    include 'mpif.h'
    m=ML_1_RHO-ML_0_RHO+1
    do s=MS_0_RHO,MS_1_RHO
       call mpi_allreduce(MPI_IN_PLACE,rho(ML_0_RHO,s),m,mpi_real8,mpi_sum,comm_band,ierr)
       call mpi_allreduce(MPI_IN_PLACE,rho(ML_0_RHO,s),m,mpi_real8,mpi_sum,comm_bzsm,ierr)
    end do
    ! The following assumes all 'MS_1-MS_0+1' are the same
    m=m*(MS_1_RHO-MS_0_RHO+1)
    call mpi_allgather(rho(ML_0_RHO,MS_0_RHO),m,mpi_real8,rho,m,mpi_real8 &
         ,comm_spin,ierr)
  END SUBROUTINE reduce_and_gather
  
#ifdef _USPP_F_TEST_
  SUBROUTINE write_rho(unit,rank)
    implicit none
    integer,intent(IN) :: unit,rank
    integer :: s,i
    if (rank==0) write(200,'(A10,4I7)') 'rho=',MS_0_RHO,MS_1_RHO,ML_0_RHO,ML_1_RHO
    write(unit+rank,'(4I7)') MS_0_RHO,MS_1_RHO,ML_0_RHO,ML_1_RHO
    do s=MS_0_RHO,MS_1_RHO
       do i=ML_0_RHO,ML_1_RHO
          write(unit+rank,'(g20.12)') rho(i,s)
       end do
    end do
  END SUBROUTINE write_rho
#endif


  SUBROUTINE calc_sum_dspin
    implicit none
    integer :: ierr
    real(8) :: tmp(2)
    sum_dspin(:)=0.0d0
    if ( MS_RHO == 2 ) then
       tmp(1)=sum(     rho(:,1)-rho(:,MS_RHO)  )*dV_RHO
       tmp(2)=sum( abs(rho(:,1)-rho(:,MS_RHO)) )*dV_RHO
       call mpi_allreduce( tmp,sum_dspin,2,mpi_real8,mpi_sum,comm_grid,ierr)
!       if ( disp_switch_parallel ) then
!          write(*,*) "sum  dspin(r)  = ",sum_dspin(1)
!          write(*,*) "sum |dspin(r)| = ",sum_dspin(2)
!       end if
    end if
  END SUBROUTINE calc_sum_dspin


END MODULE density_module
