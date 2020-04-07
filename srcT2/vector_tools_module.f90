MODULE vector_tools_module

  use rsdft_mpi_module, only: rsdft_allreduce_sum
  use parallel_module, only: ParaInfo

  implicit none

  PRIVATE
  PUBLIC :: normalize_vector
  PUBLIC :: inner_product_vector
  PUBLIC :: vinfo
  PUBLIC :: ul_check

  INTERFACE normalize_vector
     MODULE PROCEDURE d_normalize_vector, z_normalize_vector
  END INTERFACE

  INTERFACE inner_product_vector
     MODULE PROCEDURE d_inner_product, z_inner_product, &
                      d_inner_product_blas2, z_inner_product_blas2, &
                      d_inner_product_blas3, z_inner_product_blas3
  END INTERFACE

  type vinfo
     real(8) :: factor
     type(ParaInfo) :: pinfo
  end type vinfo

CONTAINS


  SUBROUTINE d_normalize_vector( v, a )
    implicit none
    real(8),intent(INOUT) :: v(:)
    type(vinfo) :: a
    real(8) :: c
    c=sum( abs(v)**2 )*a%factor
    call rsdft_allreduce_sum( c, a%pinfo%comm )
    c=1.d0/sqrt(c)
    v(:)=c*v(:)
  END SUBROUTINE d_normalize_vector

  SUBROUTINE z_normalize_vector( v, a )
    implicit none
    complex(8),intent(INOUT) :: v(:)
    type(vinfo) :: a
    real(8) :: c
    c=sum( abs(v)**2 )*a%factor
    call rsdft_allreduce_sum( c, a%pinfo%comm )
    c=1.d0/sqrt(c)
    v(:)=c*v(:)
  END SUBROUTINE z_normalize_vector


  SUBROUTINE d_inner_product( u, v, uv, a )
    implicit none
    real(8),intent(IN)     :: u(:), v(:)
    real(8),intent(OUT)    :: uv
    type(vinfo),intent(IN) :: a
    uv = sum( u(:)*v(:) )*a%factor
    call rsdft_allreduce_sum( uv, a%pinfo%comm )
  END SUBROUTINE d_inner_product

  SUBROUTINE z_inner_product( u, v, uv, a )
    implicit none
    complex(8),intent(IN)   :: u(:), v(:)
    complex(8),intent(OUT)  :: uv
    type(vinfo),intent(IN)  :: a
    uv = sum( conjg(u(:))*v(:) )*a%factor
    call rsdft_allreduce_sum( uv, a%pinfo%comm )
  END SUBROUTINE z_inner_product


  SUBROUTINE d_inner_product_blas2( a, b, c, v )
    implicit none
    real(8),intent(IN)  :: a(:,:), b(:)
    real(8),intent(OUT) :: c(:)
    type(vinfo),intent(IN) :: v
    call DGEMV( 'T',size(a,1),size(a,2),v%factor,a,size(a,1),b,1,0.0d0,c,1)
    call rsdft_allreduce_sum( c, v%pinfo%comm )
  END SUBROUTINE d_inner_product_blas2

  SUBROUTINE z_inner_product_blas2( a, b, c, v )
    implicit none
    complex(8),intent(IN)  :: a(:,:), b(:)
    complex(8),intent(OUT) :: c(:)
    type(vinfo),intent(IN) :: v
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    call ZGEMV( 'C',size(a,1),size(a,2),v%factor,a,size(a,1),b,1,z0,c,1)
    call rsdft_allreduce_sum( c, v%pinfo%comm )
  END SUBROUTINE z_inner_product_blas2


  SUBROUTINE d_inner_product_blas3( a, b, c, v )
    implicit none
    real(8),intent(IN)  :: a(:,:), b(:,:)
    real(8),intent(OUT) :: c(:,:)
    type(vinfo),intent(IN) :: v
    call DGEMM( 'T','N',size(a,2),size(b,2),size(b,1),v%factor &
                ,a,size(a,1), b,size(b,1), 0.0d0, c,size(c,1) )
    call rsdft_allreduce_sum( c, v%pinfo%comm )
  END SUBROUTINE d_inner_product_blas3

  SUBROUTINE z_inner_product_blas3( a, b, c, v )
    implicit none
    complex(8),intent(IN)  :: a(:,:), b(:,:)
    complex(8),intent(OUT) :: c(:,:)
    type(vinfo),intent(IN) :: v
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    call ZGEMM( 'C','N',size(a,2),size(b,2),size(b,1),v%factor &
                ,a,size(a,1), b,size(b,1), z0, c,size(c,1) )
    call rsdft_allreduce_sum( c, v%pinfo%comm )
  END SUBROUTINE z_inner_product_blas3


  SUBROUTINE ul_check( a )
    implicit none
    real(8),intent(IN) :: a(:,:)
    integer :: n,i,j,u1,u2,l1,l2,myrank
    real(8),parameter :: tol=1.d-10
    include 'mpif.h'
    call mpi_comm_rank(mpi_comm_world,myrank,i)
    n=size(a,1)
    u1=0.0d0
    do j=1,n
       do i=1,j
          if ( abs(a(i,j)) > tol ) u1=u1+1
       end do
    end do
    u2=0.0d0
    do j=1,n
       do i=1,n-j+1
          if ( abs(a(i,j)) > tol ) u2=u2+1
       end do
    end do
    l1=0.0d0
    do j=1,n
       do i=j,n
          if ( abs(a(i,j)) > tol ) l1=l1+1
       end do
    end do
    l2=0.0d0
    do j=1,n
       do i=n-j+1,n
          if ( abs(a(i,j)) > tol ) l2=l2+1
       end do
    end do
    write(*,'(1x,8i8)') u1,l1,u1+l1-n,u2,l2,u2+l2-n,n*(n+1)/2,myrank
  END SUBROUTINE ul_check


END MODULE vector_tools_module
