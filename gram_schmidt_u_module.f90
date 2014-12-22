MODULE gram_schmidt_u_module

  use rgrid_module
  use wf_module, only: unk
  use parallel_module
  use localpot2_Smatrix_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt_u

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: uu(:),vv(:)
  real(8),allocatable :: Sunk(:)
  real(8),parameter :: zero=0.0d0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: uu(:),vv(:)
  complex(8),allocatable :: Sunk(:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif

CONTAINS

  SUBROUTINE gram_schmidt_u(n0,n1,k,s)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    integer :: n,m,mm,ierr
    real(8) :: c,d

    mm=size( unk, 1 )
    allocate( Sunk(mm) ) ; Sunk=zero

    allocate( uu(n0:n1) ) ; uu=zero
    allocate( vv(n0:n1) ) ; vv=zero

    do n=n0,n1

       call op_localpot2_Smatrix( unk(:,n,k,s), Sunk )

       if ( n > n0 ) then

          do m=n0,n-1

#ifdef _DRSDFT_
             uu(m)=sum( unk(:,m,k,s)*Sunk(:) )*dV
#else
             uu(m)=sum( conjg(unk(:,m,k,s))*Sunk(:) )*dV
#endif
          end do ! m
          call mpi_allreduce(uu,vv,n-n0,TYPE_MAIN,MPI_SUM,comm_grid,ierr)
          do m=n0,n-1
             unk(:,n,k,s) = unk(:,n,k,s) - unk(:,m,k,s)*vv(m)
          end do ! m

       end if

       call op_localpot2_Smatrix( unk(:,n,k,s), Sunk )
#ifdef _DRSDFT_
       c = sum( unk(:,n,k,s)*Sunk(:) )*dV
#else
       c = sum( conjg(unk(:,n,k,s))*Sunk(:) )*dV
#endif
       call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       c=1.0d0/sqrt(d)
       unk(:,n,k,s)=c*unk(:,n,k,s)

    end do ! n

    deallocate( vv )
    deallocate( uu )
    deallocate( Sunk )

  END SUBROUTINE gram_schmidt_u

END MODULE gram_schmidt_u_module
