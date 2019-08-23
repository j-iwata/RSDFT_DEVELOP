module rsdft_allreduce_module

  private
  public :: test_allreduce

  integer :: n_opt, n_opt_h

contains

  SUBROUTINE test_allreduce
    implicit none
    integer :: i,j,k,m,mt,n,i0,ierr,myrank
    real(8),allocatable :: a(:), b(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
    include 'mpif.h'
 
    !return

    call write_border( 0, " test_allreduce(start)" )

    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    n=2**16
    allocate( a(n) ); a(:)=1.0d0
    allocate( b(n) ); b(:)=0.0d0
    ctmin=1.d100
    etmin=1.d100
    do j=16,0,-1
       m = 2**j
       mt=0
       do i=1,2**(16-j)
          mt=mt+m
       end do
       b=0.0d0
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct0) ; et0=mpi_wtime()
       do k=1,16
          do i=1,2**(16-j)
             i0=(i-1)*m+1
             call MPI_Allreduce(a(i0),b(i0),m,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          end do
       end do ! k
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct1) ; et1=mpi_wtime()
       ct=ct1-ct0
       et=et1-et0
       ctmin=min(ct,ctmin)
       if ( et < etmin ) then
          n_opt = m
          etmin = et
       end if
       if ( myrank == 0 ) then
          write(*,'(1x,3i12,4f10.5,2x,g15.8)') m,i-1,mt,ct,ctmin,et,etmin,sum(b(1:mt))
       end if
    end do ! j
    deallocate( b )
    deallocate( a )
    call MPI_BCAST(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h
    call write_border( 0, " test_allreduce(end)" )
  END SUBROUTINE test_allreduce

end module rsdft_allreduce_module
