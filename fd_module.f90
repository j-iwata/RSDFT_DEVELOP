MODULE fd_module

  implicit none

  PRIVATE
  PUBLIC :: get_coef_laplacian_fd,get_coef_nabla_fd

!  integer :: Md
!  real(8),allocatable :: lap(:),nab(:)

CONTAINS

!  SUBROUTINE read_fd(unit)
!    integer,intent(IN) :: unit
!    read(unit,*) Md
!    write(*,*) "Md =",Md
!  END SUBROUTINE read_fd


!  SUBROUTINE send_fd(rank)
!    integer,intent(IN) :: rank
!    integer :: ierr
!    include 'mpif.h'
!    call mpi_bcast(Md,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
!  END SUBROUTINE send_fd


  SUBROUTINE get_coef_laplacian_fd(Md,lap)
    integer,intent(IN)  :: Md
    real(8),intent(OUT) :: lap(-Md:Md)
    integer :: i,j,k
    real(8) :: t,s,s0
!    allocate( lap(-Md:Md) )
    lap=0.d0
    s=0.d0
    do i=-Md,Md
       if ( i==0 ) cycle
       do j=-Md,Md
          if ( j==i .or. j==0 ) cycle
          s=s+1.d0/dble(i*j)
       end do
    end do
    lap(0)=s
    do j=1,Md
       t=1.d0
       do i=-Md,Md
          if ( i==j ) cycle
          t=t*(j-i)
       end do
       s=0.d0
       do k=-Md,Md
          if ( k==j .or. k==0 ) cycle
          s0=1.d0
          do i=-Md,Md
             if ( i==k .or. i==j .or. i==0 ) cycle
             s0=s0*(-i)
          end do
          s=s+s0
       end do
       lap( j)=2.d0*s/t
       lap(-j)=lap(j)
    end do
  END SUBROUTINE get_coef_laplacian_fd


  SUBROUTINE get_coef_nabla_fd(Md,nab)
    integer,intent(IN)  :: Md
    real(8),intent(OUT) :: nab(-Md:Md)
    integer :: i,j
    real(8) :: t,s
!    allocate( nab(-Md:Md) )
    nab=0.d0
    do j=1,Md
       t=1.d0
       do i=-Md,Md
          if ( i==j ) cycle
          t=t*(j-i)
       end do
       s=1.d0
       do i=-Md,Md
          if ( i==j .or. i==0 ) cycle
          s=s*(-i)
       end do
       nab( j)=s/t
       nab(-j)=nab(j)
    end do
  END SUBROUTINE get_coef_nabla_fd

END MODULE fd_module
