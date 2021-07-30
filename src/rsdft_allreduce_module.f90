module rsdft_allreduce_module

  use parallel_module, only: comm_grid, comm_band, comm_bzsm, comm_spin

  private
  public :: rsdft_allreduce

  interface rsdft_allreduce
    module procedure d_rsdft_allreduce_0 &
                    ,d_rsdft_allreduce_1, z_rsdft_allreduce_1 &
                    ,d_rsdft_allreduce_2, z_rsdft_allreduce_2
  end interface rsdft_allreduce

  integer :: n_opt, n_opt_h

contains

  subroutine d_rsdft_allreduce_0( a, comm )
    implicit none
    real(8),intent(inout) :: a
    character(*),optional,intent(in) :: comm
    integer :: communicator,ierr
    real(8) :: b
    include 'mpif.h'
    communicator = MPI_COMM_WORLD
    if ( present(comm) ) communicator = set_communicator( comm )
    b=a
    call MPI_Allreduce( b, a, 1, MPI_REAL8, MPI_SUM, communicator, ierr )
  end subroutine d_rsdft_allreduce_0

  subroutine d_rsdft_allreduce_1( a, comm )
    implicit none
    real(8),intent(inout) :: a(:)
    character(*),optional,intent(in) :: comm
    integer :: communicator,n,ierr
    real(8),allocatable :: b(:)
    include 'mpif.h'
    communicator = MPI_COMM_WORLD
    if ( present(comm) ) communicator = set_communicator( comm )
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( b(n) ); b=a
    call MPI_Allreduce( b, a, n, MPI_REAL8, MPI_SUM, communicator, ierr )
#else
    call MPI_Allreduce( MPI_IN_PLACE,a,n,MPI_REAL8,MPI_SUM,communicator,ierr )
#endif
  end subroutine d_rsdft_allreduce_1

  subroutine d_rsdft_allreduce_2( a, comm )
    implicit none
    real(8),intent(inout) :: a(:,:)
    character(*),optional,intent(in) :: comm
    integer :: communicator,n,ierr
    real(8),allocatable :: b(:,:)
    include 'mpif.h'
    communicator = MPI_COMM_WORLD
    if ( present(comm) ) communicator = set_communicator( comm )
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( b(size(a,1),size(a,2)) ); b=a
    call MPI_Allreduce( b, a, n, MPI_REAL8, MPI_SUM, communicator, ierr )
#else
    call MPI_Allreduce( MPI_IN_PLACE,a,n,MPI_REAL8,MPI_SUM,communicator,ierr )
#endif
  end subroutine d_rsdft_allreduce_2

  subroutine z_rsdft_allreduce_1( a, comm )
    implicit none
    complex(8),intent(inout) :: a(:)
    character(*),optional,intent(in) :: comm
    integer :: communicator,n,ierr
    complex(8),allocatable :: b(:)
    include 'mpif.h'
    communicator = MPI_COMM_WORLD
    if ( present(comm) ) communicator = set_communicator( comm )
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( b(n) ); b=a
    call MPI_Allreduce( b,a,n,MPI_COMPLEX16,MPI_SUM,communicator,ierr )
#else
    call MPI_Allreduce( MPI_IN_PLACE,a,n,MPI_COMPLEX16,MPI_SUM,communicator,ierr )
#endif
  end subroutine z_rsdft_allreduce_1

  subroutine z_rsdft_allreduce_2( a, comm )
    implicit none
    complex(8),intent(inout) :: a(:,:)
    character(*),optional,intent(in) :: comm
    integer :: communicator,n,ierr
    complex(8),allocatable :: b(:,:)
    include 'mpif.h'
    communicator = MPI_COMM_WORLD
    if ( present(comm) ) communicator = set_communicator( comm )
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( b(size(a,1),size(a,2)) ); b=a
    call MPI_Allreduce( b,a,n,MPI_COMPLEX16,MPI_SUM,communicator,ierr )
#else
    call MPI_Allreduce( MPI_IN_PLACE,a,n,MPI_COMPLEX16,MPI_SUM,communicator,ierr )
#endif
  end subroutine z_rsdft_allreduce_2


  integer function set_communicator( comm_in )
    implicit none
    character(*),intent(in) :: comm_in
    character(9) :: comm
    comm=comm_in
    call convertToLowercase( comm )
    select case( comm )
    case( 'g', 'grid', 'comm_grid' ); set_communicator = comm_grid
    case( 'b', 'band', 'comm_band' ); set_communicator = comm_band
    case( 'k', 'bzsm', 'comm_bzsm' ); set_communicator = comm_bzsm
    case( 's', 'spin', 'comm_spin' ); set_communicator = comm_spin
    case default
      write(*,*) "comm= ",comm," is not defined."
      call stop_program('set_communicator@rsdft_allreduce_module')
    end select
  end function set_communicator

  subroutine convertToLowercase( cbuf, ckey )
    implicit none
    character(*),intent(inout) :: cbuf
    character(*),optional,intent(out) :: ckey
    integer :: j,k,n
    n=len_trim(cbuf)
    if ( present(ckey) ) ckey=cbuf(1:n)
    do j=1,n
      k=iachar( cbuf(j:j) )
      if ( 65 <= k .and. k <= 90 ) k=k+32
      if ( present(ckey) ) then
        ckey(j:j) = achar(k)
      else
        cbuf(j:j) = achar(k)
      end if
    end do
  end subroutine convertToLowercase

  subroutine test_allreduce
    implicit none
    integer :: i,j,k,m,mt,n,i0,ierr,myrank
    real(8),allocatable :: a(:), b(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
    include 'mpif.h'
 
    !return

    call write_border( 0, " test_allreduce(start)" )

    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
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
    call MPI_Bcast(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h
    call write_border( 0, " test_allreduce(end)" )
  end subroutine test_allreduce

end module rsdft_allreduce_module
