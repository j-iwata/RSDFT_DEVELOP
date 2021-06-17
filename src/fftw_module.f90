module fftw_module

  use,intrinsic :: iso_c_binding
  use grid_module, only: grid, get_range_rgrid, mpi_allgatherv_grid, get_map_3d_to_1d_grid
  use rsdft_mpi_module, only: rsdft_allgatherv

  implicit none

  private
  public :: init_fftw
  public :: finalize_fftw
  public :: plan_forward, plan_backward
  public :: ML1_c, ML2_c, N_ML3_c, ML3_c0
  public :: zwork3_ptr0, zwork3_ptr1
  public :: d1_to_z3_fftw
  public :: z3_to_d1_fftw
  public :: z3_to_z1_fftw
  public :: forward_fftw
  public :: backward_fftw
  public :: forward_2d_fftw

  integer :: comm_fftw
  type(c_ptr) :: plan_forward,plan_backward
  type(c_ptr) :: zwork3_cptr0,zwork3_cptr1
  integer(c_intptr_t) :: ML1_c,ML2_c,ML3_c,ML3_c0,N_ML3_c
  complex(c_double_complex), pointer :: zwork3_ptr0(:,:,:)=>null()
  complex(c_double_complex), pointer :: zwork3_ptr1(:,:,:)=>null()

  integer,allocatable :: ir(:),id(:)

  integer(8) :: plan_2d_forward
  integer :: size_2d(2)=(/0,0/)

  logical :: flag_init_fftw=.false.

contains


  SUBROUTINE init_fftw( Ngrid, Np, comm_grid, myrank_g )

    implicit none
    integer,intent(IN) :: Ngrid(3),Np(3),comm_grid, myrank_g
#ifdef _FFTW_
    integer :: ML1,Ml2,ML3,np1,np2,np3,n,m
    integer :: i1,i2,i3,irank,icolor,ierr
    integer(c_intptr_t) :: alloc_ML3_c
    include "fftw3-mpi.f03"
    include "mpif.h"

    call write_border( 0, " init_fftw(start)" )

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    np1 = Np(1)
    np2 = Np(2)
    np3 = Np(3)

    if ( mod(ML3,np3) /= 0 ) then
      call stop_program( "fftw_init: mod(ML3,np3)/=0 is not supported." )
    end if

    irank=-1
    do i3=1,np3
    do i2=1,np2
    do i1=1,np1
      irank=irank+1
      if ( irank == myrank_g ) icolor=i1+(i2-1)*np1
    end do
    end do
    end do

    call mpi_comm_split(comm_grid,icolor,0,comm_fftw,ierr)

    call fftw_mpi_init()

    ML1_c = ML1
    ML2_c = ML2
    ML3_c = ML3

    alloc_ML3_c = &
         fftw_mpi_local_size_3d(ML3_c,ML2_c,ML1_c,comm_fftw,N_ML3_c,ML3_c0)
    zwork3_cptr0=fftw_alloc_complex(alloc_ML3_c)
    zwork3_cptr1=fftw_alloc_complex(alloc_ML3_c)
    call c_f_pointer(zwork3_cptr0, zwork3_ptr0, [ML1_c,ML2_c,N_ML3_c])
    call c_f_pointer(zwork3_cptr1, zwork3_ptr1, [ML1_c,ML2_c,N_ML3_c])
    zwork3_ptr0=(0.0d0,0.0d0)
    zwork3_ptr1=(0.0d0,0.0d0)

    plan_forward  = fftw_mpi_plan_dft_3d( ML3_c,ML2_c,ML1_c &
         ,zwork3_ptr0,zwork3_ptr1,comm_fftw,fftw_forward,fftw_measure )

    plan_backward = fftw_mpi_plan_dft_3d( ML3_c,ML2_c,ML1_c &
         ,zwork3_ptr0,zwork3_ptr1,comm_fftw,fftw_backward,fftw_measure )

    call mpi_comm_size( comm_fftw, n, ierr )
    call mpi_comm_rank( comm_fftw, m, ierr )
    allocate( ir(0:n-1) ) ; ir=0
    allocate( id(0:n-1) ) ; id=0
!    m=ML1*ML2
!    call mpi_allgather(m*N_ML3_c,1,MPI_INTEGER,ir,1,MPI_INTEGER,comm_fftw,ierr)
!    call mpi_allgather(m*ML3_c0 ,1,MPI_INTEGER,id,1,MPI_INTEGER,comm_fftw,ierr)
    i1 = ML1*ML2*N_ML3_c
    i2 = ML1*ML2*ML3_c0
    call mpi_allgather(i1,1,MPI_INTEGER,ir,1,MPI_INTEGER,comm_fftw,ierr)
    call mpi_allgather(i2,1,MPI_INTEGER,id,1,MPI_INTEGER,comm_fftw,ierr)

    flag_init_fftw=.true.

    call write_border( 0, " init_fftw(end)" )
#endif
  END SUBROUTINE init_fftw


  SUBROUTINE finalize_fftw
#ifdef _FFTW_
    implicit none
    integer :: ierr
    include "fftw3-mpi.f03"
    if ( .not.flag_init_fftw ) return
    call mpi_comm_free(comm_fftw,ierr)
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
    call fftw_free(zwork3_cptr0)
    call fftw_free(zwork3_cptr1)
    call fftw_mpi_cleanup()
#endif
  END SUBROUTINE finalize_fftw


  subroutine d1_to_z3_fftw( d1, z3 )
    implicit none
    real(8),intent(in) :: d1(:)
    complex(8),intent(out) :: z3(:,:,:)
    integer :: i1,i2,i3
    type(grid) :: rgrid
    real(8),allocatable :: work(:)
    integer,allocatable :: LLL(:,:,:)
    call get_range_rgrid( rgrid )
    allocate( work(rgrid%g1%size_global) ); work=0.0d0
    call mpi_allgatherv_grid( d1, work )
    call get_map_3d_to_1d_grid( rgrid, LLL )
    do i3=1,N_ML3_c
    do i2=1,ML2_c
    do i1=1,ML1_c
      z3(i1,i2,i3) = work( LLL(i1-1,i2-1,i3-1+ML3_c0) )
    end do
    end do
    end do
    deallocate( LLL )
    deallocate( work )
  end subroutine d1_to_z3_fftw

  subroutine z3_to_d1_fftw( z3, d1 )
    implicit none
    complex(8),intent(in) :: z3(:,:,:)
    real(8),intent(out) :: d1(:)
    integer :: i1,i2,i3,i
    type(grid) :: rgrid
    call get_range_rgrid( rgrid )
    i=0
    do i3=1,N_ML3_c
    do i2=rgrid%g3%y%head+1,rgrid%g3%y%tail+1
    do i1=rgrid%g3%x%head+1,rgrid%g3%x%tail+1
       i=i+1
       d1(i) = dble( z3(i1,i2,i3) )
    end do
    end do
    end do
  end subroutine z3_to_d1_fftw

  SUBROUTINE z3_to_z1_fftw( z3, z1 )
    implicit none
    complex(8),intent(IN) :: z3(:,:,:)
    complex(8),intent(OUT) :: z1(:)
    integer :: i1,i2,i3,i
    type(grid) :: rgrid
    call get_range_rgrid( rgrid )
    i=0
    do i3=1,N_ML3_c
    do i2=rgrid%g3%y%head+1,rgrid%g3%y%tail+1
    do i1=rgrid%g3%x%head+1,rgrid%g3%x%tail+1
       i=i+1
       z1(i) = z3(i1,i2,i3)
    end do
    end do
    end do
  END SUBROUTINE z3_to_z1_fftw


  SUBROUTINE forward_fftw( z3 )
    implicit none
    complex(8),intent(INOUT) :: z3(:,:,:)
    integer :: i1,i2,i3,j3
    real(8) :: const
#ifdef _FFTW_
    include "fftw3-mpi.f03"
    do i3=1,N_ML3_c
       j3=i3+ML3_c0
       do i2=1,ML2_c
       do i1=1,ML1_c
          zwork3_ptr0(i1,i2,i3) = z3(i1,i2,j3)
       end do
       end do
    end do
    call fftw_mpi_execute_dft( plan_forward, zwork3_ptr0, zwork3_ptr1 )
    call rsdft_allgatherv( zwork3_ptr1, z3, ir, id, comm_fftw )
    const=1.0d0/size(z3)
    z3(:,:,:)=const*z3(:,:,:)
#endif
  END SUBROUTINE forward_fftw


  SUBROUTINE backward_fftw( z3 )
    implicit none
    complex(8),intent(INOUT) :: z3(:,:,:)
    integer :: i1,i2,i3,j3
    real(8) :: const
#ifdef _FFTW_
    include "fftw3-mpi.f03"
    do i3=1,N_ML3_c
       j3=i3+ML3_c0
       do i2=1,ML2_c
       do i1=1,ML1_c
          zwork3_ptr0(i1,i2,i3) = z3(i1,i2,j3)
       end do
       end do
    end do
    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )
    call rsdft_allgatherv( zwork3_ptr1, z3, ir, id, comm_fftw )
    !const=1.0d0/size(z3)
    !z3(:,:,:)=const*z3(:,:,:)
#endif
  END SUBROUTINE backward_fftw


  SUBROUTINE forward_2d_fftw( z2, zw )
    implicit none
    complex(8),intent(INOUT) :: z2(:,:), zw(:,:)
    integer :: m1,m2
    real(8) :: const
#ifdef _FFTW_
    include 'fftw3.f'
    m1 = size( z2, 1 )
    m2 = size( z2, 2 )
    if ( size_2d(1) /= m1 .or. size_2d(2) /= m2 ) then
       call dfftw_plan_dft_2d( plan_2d_forward, m1, m2, z2, z2 &
            , FFTW_FORWARD, FFTW_ESTIMATE )
       size_2d(:) = (/ m1, m2 /)
    end if
    call dfftw_execute_dft( plan_2d_forward, z2, zw )
    const=1.0d0/size(z2)
    z2(:,:)=const*zw(:,:)
#endif
  END SUBROUTINE forward_2d_fftw


END MODULE fftw_module
