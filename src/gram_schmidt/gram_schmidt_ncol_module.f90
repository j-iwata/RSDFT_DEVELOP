MODULE gram_schmidt_ncol_module

  use rgrid_module, only: dV
  use parallel_module, only: comm_grid, id_bzsm, myrank_k

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt_ncol

CONTAINS

  SUBROUTINE gram_schmidt_ncol( n0,n1, k0, psi )
    implicit none
    integer,intent(IN) :: n0,n1,k0
    complex(8),intent(INOUT) :: psi(:,:,:,:)
    complex(8) :: uu,vv
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer :: n,m,s,MS,ierr,k
    real(8) :: c,d
    include 'mpif.h'

    call write_border( 1, "gram_schmidt_ncol(start)" )

    MS = size( psi, 4 )
    k  = k0-id_bzsm(myrank_k)

    do n=n0,n1

       do m=n0,n-1
          uu=zero
          do s=1,MS
             uu = uu + sum( conjg(psi(:,m,k,s))*psi(:,n,k,s) )*dV
          end do
          call mpi_allreduce(uu,vv,1,MPI_COMPLEX16,MPI_SUM,comm_grid,ierr)
          do s=1,MS
             psi(:,n,k,s) = psi(:,n,k,s) - psi(:,m,k,s)*vv
          end do
       end do ! m

       c=0.0d0
       do s=1,MS
          c = c + sum( abs(psi(:,n,k,s))**2 )*dV
       end do
       call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       c=1.0d0/sqrt(d)
       do s=1,MS
          psi(:,n,k,s)=c*psi(:,n,k,s)
       end do

    end do ! n

    call write_border( 1, "gram_schmidt_ncol(end)" )

  END SUBROUTINE gram_schmidt_ncol

END MODULE gram_schmidt_ncol_module
