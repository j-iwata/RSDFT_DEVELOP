module pzfft3dv_test_module

  use ffte_sub_module, only: comm_fftx, comm_ffty, comm_fftz
  use rsdft_mpi_module, only: rsdft_allgather, rsdft_allgatherv

  implicit none

  private
  public :: init_pzfft3dv_test
  public :: pzfft3dv_test

  complex(8),public,allocatable :: zwork1_ffte(:,:,:)
  complex(8),public,allocatable :: zwork2_ffte(:,:,:)

  integer,allocatable :: ircnx(:), ircny(:), ircnz(:)
  integer,allocatable :: idisx(:), idisy(:), idisz(:)

  integer :: mmx, mmy, mmz
  integer :: nnx, nny, nnz
  integer :: nx, ny, nz
  integer :: nx_0,nx_1,ny_0,ny_1,nz_0,nz_1

  complex(8),allocatable :: a(:)
  complex(8),allocatable :: b(:)
  complex(8),allocatable :: work(:)

  logical :: flag_init_done = .false.

contains

  subroutine init_pzfft3dv_test( nx_in, ny_in, nz_in )
    implicit none
    integer,intent(in) :: nx_in, ny_in, nz_in
    integer :: npx,npy,npz,ierr,i,myrnkx,myrnky,myrnkz
    logical :: disp
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    include 'mpif.h'

    if ( flag_init_done ) return

    call write_border( 0, 'init_pzfft3dv_test(start)' )

    nx = nx_in
    ny = ny_in
    nz = nz_in

    call MPI_Comm_size( comm_fftx, npx, ierr )
    call MPI_Comm_size( comm_ffty, npy, ierr )
    call MPI_Comm_size( comm_fftz, npz, ierr )
    call MPI_Comm_rank( comm_fftx, myrnkx, ierr )
    call MPI_Comm_rank( comm_ffty, myrnky, ierr )
    call MPI_Comm_rank( comm_fftz, myrnkz, ierr )

    allocate( idisx(0:npx-1) ); idisx=0
    allocate( idisy(0:npy-1) ); idisy=0
    allocate( idisz(0:npz-1) ); idisz=0
    allocate( ircnx(0:npx-1) ); ircnx=0
    allocate( ircny(0:npy-1) ); ircny=0
    allocate( ircnz(0:npz-1) ); ircnz=0

    call rsdft_allgather( nx, ircnx, comm_fftx )
    call rsdft_allgather( ny, ircny, comm_ffty )
    call rsdft_allgather( nz, ircnz, comm_fftz )

    do i = 0, npx-1
      idisx(i) = sum( ircnx(0:i) ) - ircnx(i)
    end do
    do i = 0, npy-1
      idisy(i) = sum( ircny(0:i) ) - ircny(i)
    end do
    do i = 0, npz-1
      idisz(i) = sum( ircnz(0:i) ) - ircnz(i)
    end do

    nnx = sum(ircnx)
    nny = sum(ircny)
    nnz = sum(ircnz)

    call get_n235( nnx, mmx )
    call get_n235( nny, mmy )
    call get_n235( nnz, mmz )

    nx_0 = idisx(myrnkx) + 1
    nx_1 = nx_0 + ircnx(myrnkx) - 1
    ny_0 = idisy(myrnky) + 1
    ny_1 = ny_0 + ircny(myrnky) - 1
    nz_0 = idisz(myrnkz) + 1
    nz_1 = nz_0 + ircnz(myrnkz) - 1

    call check_disp_switch( disp, 0 )
    if ( disp ) then
      write(*,'("ircnx",10i4)') ircnx
      write(*,'("ircny",10i4)') ircny
      write(*,'("ircnz",10i4)') ircnz
      write(*,*) "mmx,mmy,mmz=",mmx,mmy,mmz
    end if

    i = max( nnx, nny, nnz )
    allocate( a(i)      ); a=z0
    allocate( b(i)      ); b=z0
    allocate( work(2*i) ); work=z0

    allocate( zwork1_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork1_ffte=z0
    allocate( zwork2_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork2_ffte=z0

    flag_init_done = .true.

    call write_border( 0, 'init_pzfft3dv_test(end)' )

  end subroutine init_pzfft3dv_test

  subroutine pzfft3dv_test( f, iopt )
    implicit none
    complex(8),intent(inout) :: f(:,:,:)
    integer,intent(in) :: iopt
    integer :: ix,iy,iz

    call write_border( 0, 'pzfft3dv_test(start)' )
    if ( .not.flag_init_done ) call stop_program('Call "init_pzfft3dv_test" first.')

    call ZFFT1D( a, nnx, 0, work ) 
    do iz=1,nz
    do iy=1,ny
      a(1:nx) = f(1:nx,iy,iz)
      call rsdft_allgatherv( a(1:nx), b(1:nnx), ircnx, idisx, comm_fftx )
      call ZFFT1D( b, nnx, iopt, work ) 
      f(1:nx,iy,iz) = b(nx_0:nx_1)
    end do
    end do

    call ZFFT1D( a, nny, 0, work ) 
    do iz=1,nz
    do ix=1,nx
      a(1:ny) = f(ix,1:ny,iz)
      call rsdft_allgatherv( a(1:ny), b(1:nny), ircny, idisy, comm_ffty )
      call ZFFT1D( b, nny, iopt, work )
      f(ix,1:ny,iz) = b(ny_0:ny_1)
    end do
    end do

    call ZFFT1D( a, nnz, 0, work ) 
    do iy=1,ny
    do ix=1,nx
      a(1:nz) = f(ix,iy,1:nz)
      call rsdft_allgatherv( a(1:nz), b(1:nnz), ircnz, idisz, comm_fftz )
      call ZFFT1D( b, nnz, iopt, work )
      f(ix,iy,1:nz) = b(nz_0:nz_1)
    end do
    end do

    call write_border( 0, 'pzfft3dv_test(end)' )

  end subroutine pzfft3dv_test

  subroutine get_n235( n, n235 )
    implicit none
    integer,intent(in) :: n
    integer,intent(out) :: n235
    integer :: n2,n3,n5,i2,i3,i5,m
    n2 = nint( log(dble(n))/log(2.0d0) ) + 1
    n3 = nint( log(dble(n))/log(3.0d0) ) + 1
    n5 = nint( log(dble(n))/log(5.0d0) ) + 1
    n235=10000000
    do i5 = 0, n5
    do i3 = 0, n3
    do i2 = 0, n2
      m = 2**i2 * 3**i3 * 5**i5
      if ( m >= n ) n235=min(n235,m)
    end do
    end do
    end do
  end subroutine get_n235

end module pzfft3dv_test_module
