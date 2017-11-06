SUBROUTINE allreduce_sub( n, comm, a, num )

  implicit none

  include 'mpif.h'
  integer,intent(IN) :: n, comm, num
  real(8),intent(INOUT) :: a(n)
  real(8),allocatable :: b(:)
  integer :: ierr,i,j,m,nd

  ierr=0

  if ( num <= 1 .or. num > n ) then  

     allocate( b(n) ) ; b=0.0d0

     call mpi_allreduce( a, b, n, MPI_REAL8, MPI_SUM, comm, ierr )

     a=b

  else

     nd = (n+num)/num

     allocate( b(nd) ) ; b=0.0d0

     do i=1,n,nd
        j=min(i+nd-1,n)
        m=j-i+1
        call mpi_allreduce( a(i), b, m, MPI_REAL8, MPI_SUM, comm, ierr )
        a(i:j)=b
     end do

  end if

  deallocate( b )

END SUBROUTINE allreduce_sub
