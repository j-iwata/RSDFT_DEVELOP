module rsdft_bcast_module

  use communicator_module, only: set_communicator

  implicit none

  private
  public :: l_rsdft_bcast, l_rsdft_bcast_tmp
  public :: i_rsdft_bcast, i_rsdft_bcast_tmp
  public :: d_rsdft_bcast, d_rsdft_bcast_tmp
  public :: z_rsdft_bcast_tmp
  public :: test_bcast

  integer :: n_opt = 1024
  integer :: n_opt_h = 1024

  interface l_rsdft_bcast
    module procedure l_rsdft_bcast_0, l_rsdft_bcast_1
  end interface l_rsdft_bcast

  interface i_rsdft_bcast
    module procedure i_rsdft_bcast_0, i_rsdft_bcast_1
  end interface i_rsdft_bcast

contains

  subroutine l_rsdft_bcast_0( a, root, comm_chr, comm_int )
    implicit none
    logical,intent(inout) :: a
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: n,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_LOGICAL
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    n = 1
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
#endif
  end subroutine l_rsdft_bcast_0

  subroutine l_rsdft_bcast_1( a, root, comm_chr, comm_int )
    implicit none
    logical,intent(inout) :: a(:)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,n,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_LOGICAL
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine l_rsdft_bcast_1

  subroutine l_rsdft_bcast_tmp( a, n, root, comm_chr, comm_int )
    implicit none
    integer,intent(in) :: n
    logical,intent(inout) :: a(n)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_LOGICAL
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    ! n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine l_rsdft_bcast_tmp

  subroutine i_rsdft_bcast_0( a, root, comm_chr, comm_int )
    implicit none
    integer,intent(inout) :: a
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,n,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_INTEGER
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    n = 1
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
#endif
  end subroutine i_rsdft_bcast_0

  subroutine i_rsdft_bcast_1( a, root, comm_chr, comm_int )
    implicit none
    integer,intent(inout) :: a(:)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,n,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_INTEGER
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine i_rsdft_bcast_1
  
  subroutine i_rsdft_bcast_tmp( a, n, root, comm_chr, comm_int )
    implicit none
    integer,intent(in) :: n
    integer,intent(inout) :: a(n)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_INTEGER
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    ! n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine i_rsdft_bcast_tmp


  subroutine d_rsdft_bcast( a, root, comm_chr, comm_int )
    implicit none
    real(8),intent(inout) :: a(:)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,n,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_REAL8
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine d_rsdft_bcast

  subroutine d_rsdft_bcast_tmp( a, n, root, comm_chr, comm_int )
    implicit none
    integer,intent(in) :: n
    real(8),intent(inout) :: a(n)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_REAL8
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    ! n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine d_rsdft_bcast_tmp


  subroutine z_rsdft_bcast_tmp( a, n, root, comm_chr, comm_int )
    implicit none
    integer,intent(in) :: n
    complex(8),intent(inout) :: a(n)
    integer,optional,intent(in) :: root
    character(*),optional,intent(in) :: comm_chr
    integer,optional,intent(in) :: comm_int
    integer :: m,i,type,comm,ierr,rt
#ifndef _NOMPI_
    include 'mpif.h'
    type = MPI_COMPLEX16
    comm = MPI_COMM_WORLD
    rt = 0
    if ( present(comm_chr) ) comm = set_communicator( comm_chr )
    if ( present(comm_int) ) comm = comm_int
    if ( present(root) ) rt = root
    ! n = size(a)
    call MPI_Bcast(a,n,type,rt,comm,ierr); return
    do i = 1, n, n_opt
      m = min( n-i+1, n_opt )
      call MPI_Bcast(a(i),m,type,rt,comm,ierr)
    end do
#endif
  end subroutine z_rsdft_bcast_tmp


  subroutine test_bcast
    implicit none
    integer :: i,j,k,m,mt,n,ierr,myrank,npow
    real(8),allocatable :: a(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
#ifndef _NOMPI_
    include 'mpif.h'

    return

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    npow=20
    n=2**npow
    allocate( a(n) ); a(:)=1.0d0
    ctmin=1.d100
    etmin=1.d100
    do j=0,npow
      m = 2**j
      mt=0
      do i=1,2**(npow-j)
        mt=mt+m
      end do
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call cpu_time(ct0); et0=mpi_wtime()
      do k=1,10
        do i=1,2**(npow-j)
          call MPI_Bcast(a(i),m,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        end do
      end do ! k
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call cpu_time(ct1); et1=mpi_wtime()
      ct=ct1-ct0
      et=et1-et0
      ctmin=min(ct,ctmin)
      if ( et < etmin ) then
        n_opt = m
        etmin = et
      end if
      if ( myrank == 0 ) then
        write(*,'(1x,3i12,4f10.5)') m,i-1,mt,ct,ctmin,et,etmin
      end if
    end do ! j
    deallocate( a )
    call MPI_Bcast(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h
#endif
  end subroutine test_bcast    

end module rsdft_bcast_module
