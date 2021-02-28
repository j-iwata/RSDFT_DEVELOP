module pzfft3dv_arb_module

  use ffte_sub_module, only: comm_fftx,comm_ffty,comm_fftz,zwork1_ffte,zwork2_ffte
  use rsdft_mpi_module, only: rsdft_allgather, rsdft_allgatherv

  implicit none

  private
  public :: init_pzfft3dv_arb
  public :: pzfft3dv_arb

  public :: zwork1_ffte, zwork2_ffte

  integer,allocatable :: ircnx(:), ircny(:), ircnz(:)
  integer,allocatable :: idisx(:), idisy(:), idisz(:)

  integer :: mmx, mmy, mmz
  integer :: nnx, nny, nnz
  integer :: nx, ny, nz
  integer :: nx_0,nx_1,ny_0,ny_1,nz_0,nz_1
  real(8) :: inv_ngrd

  complex(8) :: wx, wy, wz

  complex(8),allocatable :: a(:)
  complex(8),allocatable :: b(:)
  complex(8),allocatable :: c(:)
  complex(8),allocatable :: work(:)

  logical :: flag_init_done = .false.

contains

  subroutine init_pzfft3dv_arb( nx_in, ny_in, nz_in )
    implicit none
    integer,intent(in) :: nx_in, ny_in, nz_in
    integer :: npx,npy,npz,ierr,i,myrnkx,myrnky,myrnkz
    logical :: disp
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8) :: pi
    include 'mpif.h'

    if ( flag_init_done ) return

    call write_border( 0, 'init_pzfft3dv_arb(start)' )

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

    inv_ngrd = 1.0d0/(nnx*nny*nnz)

    call get_n235( 2*nnx, mmx )
    call get_n235( 2*nny, mmy )
    call get_n235( 2*nnz, mmz )

    !if ( mod(mmx,nnx) == 0 ) mmx=nnx
    !if ( mod(mmy,nny) == 0 ) mmy=nny
    !if ( mod(mmz,nnz) == 0 ) mmz=nnz

    pi = acos(-1.0d0)
    wx = dcmplx( cos(pi/nnx), -sin(pi/nnx) )
    wy = dcmplx( cos(pi/nny), -sin(pi/nny) )
    wz = dcmplx( cos(pi/nnz), -sin(pi/nnz) )

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
      write(*,*) "nnx,nny,nnz=",nnx,nny,nnz
      write(*,*) "mmx,mmy,mmz=",mmx,mmy,mmz
    end if

    i = max( nnx, nny, nnz, mmx, mmy, mmz )
    allocate( a(i)      ); a=z0
    allocate( b(i)      ); b=z0
    allocate( c(i)      ); c=z0
    allocate( work(2*i) ); work=z0

    if ( allocated(zwork2_ffte) ) deallocate(zwork2_ffte)
    if ( allocated(zwork1_ffte) ) deallocate(zwork1_ffte)
    allocate( zwork1_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork1_ffte=z0
    allocate( zwork2_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork2_ffte=z0

    flag_init_done = .true.

    call write_border( 0, 'init_pzfft3dv_arb(end)' )

  end subroutine init_pzfft3dv_arb

  subroutine pzfft3dv_arb( f, iopt )
    implicit none
    complex(8),intent(inout) :: f(:,:,:)
    integer,intent(in) :: iopt
    integer :: ix,iy,iz,i,j

    complex(8),allocatable :: x(:),y(:),z(:),w(:),x1(:),x2(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    complex(8) :: ww

    call write_border( 0, 'pzfft3dv_arb(start)' )
    if ( .not.flag_init_done ) call stop_program('Call "init_pzfft3dv_arb" first.')

! ---

    ww=wx ; if ( iopt == 1 ) ww=conjg(ww)
    call ZFFT1D( a, mmx, 0, work ) 
    do iz=1,nz
    do iy=1,ny
      a(1:nx) = f(1:nx,iy,iz)
      call rsdft_allgatherv( a(1:nx), b(1:nnx), ircnx, idisx, comm_fftx )
      call ZFFT1D_arb( b(1:mmx),a(1:mmx),c(1:mmx),work(1:2*mmx),ww,mmx,nnx,nx_0-1,nx_1-1 )
      f(1:nx,iy,iz) = b(nx_0:nx_1)
    end do
    end do

    !call ZFFT1D( a, nnx, 0, work ) 
    !do iz=1,nz
    !do iy=1,ny
    !  a(1:nx) = f(1:nx,iy,iz)
    !  call rsdft_allgatherv( a(1:nx), b(1:nnx), ircnx, idisx, comm_fftx )
    !  call ZFFT1D( b, nnx, iopt, work ) 
    !  f(1:nx,iy,iz) = b(nx_0:nx_1)
    !end do
    !end do

! ---

    ww=wy ; if ( iopt == 1 ) ww=conjg(ww)
    call ZFFT1D( a, mmy, 0, work ) 
    do iz=1,nz
    do ix=1,nx
      a(1:ny) = f(ix,1:ny,iz)
      call rsdft_allgatherv( a(1:ny), b(1:nny), ircny, idisy, comm_ffty )
      call ZFFT1D_arb( b(1:mmy),a(1:mmy),c(1:mmy),work(1:2*mmy),ww,mmy,nny,ny_0-1,ny_1-1 )
      f(ix,1:ny,iz) = b(ny_0:ny_1)
    end do
    end do

    !call ZFFT1D( a, nny, 0, work ) 
    !do iz=1,nz
    !do ix=1,nx
    !  a(1:ny) = f(ix,1:ny,iz)
    !  call rsdft_allgatherv( a(1:ny), b(1:nny), ircny, idisy, comm_ffty )
    !  call ZFFT1D( b, nny, iopt, work )
    !  f(ix,1:ny,iz) = b(ny_0:ny_1)
    !end do
    !end do

! ---

    ww=wz ; if ( iopt == 1 ) ww=conjg(ww)
    call ZFFT1D( a, mmz, 0, work ) 
    do iy=1,ny
    do ix=1,nx
      a(1:nz) = f(ix,iy,1:nz)
      call rsdft_allgatherv( a(1:nz), b(1:nnz), ircnz, idisz, comm_fftz )
      call ZFFT1D_arb( b(1:mmz),a(1:mmz),c(1:mmz),work(1:2*mmz),ww,mmz,nnz,nz_0-1,nz_1-1 )
      f(ix,iy,1:nz) = b(nz_0:nz_1)
    end do
    end do

    !call ZFFT1D( a, nnz, 0, work ) 
    !do iy=1,ny
    !do ix=1,nx
    !  a(1:nz) = f(ix,iy,1:nz)
    !  call rsdft_allgatherv( a(1:nz), b(1:nnz), ircnz, idisz, comm_fftz )
    !  call ZFFT1D( b, nnz, iopt, work )
    !  f(ix,iy,1:nz) = b(nz_0:nz_1)
    !end do
    !end do

! ---

    if ( iopt == 1 ) f(:,:,:) = f(:,:,:) * inv_ngrd

    call write_border( 0, 'pzfft3dv_arb(end)' )

  contains

    subroutine DFTJI( b ) !( x --> G )
      implicit none
      complex(8),intent(inout) :: b(0:)
      complex(8) :: wn
      complex(8),allocatable :: bG(:)
      complex(8),parameter :: z0=(0.0d0,0.0d0)
      integer :: n,i,j
      real(8) :: pi2
      n = size(b)
      pi2 = 2.0d0*acos(-1.0d0)
      wn = dcmplx( cos(pi2/n), -sin(pi2/n) )
      allocate( bG(0:n-1) ); bG=z0
      do j = 0, n-1
        do i = 0, n-1
          bG(j) = bG(j) + b(i)*wn**(-i*j)
        end do
      end do
      b = bG
      deallocate( bG )
    end subroutine DFTJI

  end subroutine pzfft3dv_arb

  subroutine ZFFT1D_arb( b, a, c, work, ww, mm, nn, nn_0, nn_1 )
    implicit none
    complex(8),intent(inout) :: b(0:)
    complex(8),intent(inout) :: a(0:), c(0:)
    complex(8),intent(inout) :: work(0:)
    complex(8),intent(in) :: ww
    integer,intent(in) :: mm,nn,nn_0,nn_1
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    integer :: i,j
    a(:) = z0
    do i = -nn+1, nn-1
      j = mod( i+mm, mm )
      a(j) = ww**(-i*i)
    end do
    do i = 0, nn-1
      b(i) = b(i) * ww**(i*i)
    end do
    if ( size(b) > nn ) b(nn:)=z0
    call ZFFT1D( a, mm,-1, work ) 
    call ZFFT1D( b, mm,-1, work ) 
    c = a*b
    call ZFFT1D( c, mm, 1, work ) 
    do i = 0, nn-1
      c(i) = c(i) * ww**(i*i)
    end do
    do i = nn_0, nn_1
      b(i) = c(i)
    end do
  end subroutine ZFFT1D_arb

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

end module pzfft3dv_arb_module
