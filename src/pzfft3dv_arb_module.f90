module pzfft3dv_arb_module

  use ffte_sub_module, only: comm_fftx,comm_ffty,comm_fftz,zwork1_ffte,zwork2_ffte &
                            ,npux,npuy,npuz,me_fftx,me_ffty,me_fftz
  use rsdft_mpi_module, only: rsdft_allgather, rsdft_allgatherv
  !use io_tools_module, only: IOTools_findKeyword

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
  real(8) :: inv_nnx, inv_nny, inv_nnz
  logical :: flag_x235 = .false.
  logical :: flag_y235 = .false.
  logical :: flag_z235 = .false.
  logical :: use_zfft1d_x = .true.
  logical :: use_zfft1d_y = .true.
  logical :: use_zfft1d_z = .true.

  complex(8) :: wx, wy, wz
  complex(8),allocatable :: wxa(:),wxb(:),wya(:),wyb(:),wza(:),wzb(:)

  complex(8),allocatable :: a(:)
  complex(8),allocatable :: b(:)
  complex(8),allocatable :: c(:)
  complex(8),allocatable :: work(:)

  logical :: flag_init_done = .false.

contains

  subroutine init_pzfft3dv_arb( nx_in, ny_in, nz_in )
    implicit none
    integer,intent(in) :: nx_in, ny_in, nz_in
    integer :: npx,npy,npz,ierr,i,j,myrnkx,myrnky,myrnkz
    logical :: disp !, flag_zfft1d
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8) :: pi

    if ( flag_init_done ) return

    call write_border( 0, 'init_pzfft3dv_arb(start)' )

    !call IOTools_findKeyword( 'USE_ZFFT1D', flag_zfft1d, flag_bcast=.true. )

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

    inv_nnx = 1.0d0/nnx
    inv_nny = 1.0d0/nny
    inv_nnz = 1.0d0/nnz

    call get_n235( 2*nnx, mmx )
    call get_n235( 2*nny, mmy )
    call get_n235( 2*nnz, mmz )

    if ( mod(mmx,nnx) == 0 ) then
      mmx=nnx
      flag_x235 = .true.
      !if ( mod(nnx,npx**2) == 0 ) use_zfft1d_x=.false.
    end if
    if ( mod(mmy,nny) == 0 ) then
      mmy=nny
      flag_y235 = .true.
      !if ( mod(nny,npy**2) == 0 ) use_zfft1d_y=.false.
    end if
    if ( mod(mmz,nnz) == 0 ) then
      mmz=nnz
      flag_z235 = .true.
      !if ( mod(nnz,npz**2) == 0 ) use_zfft1d_z=.false.
    end if

    !if ( flag_zfft1d ) then
    !  use_zfft1d_x = .true.
    !  use_zfft1d_y = .true.
    !  use_zfft1d_z = .true.
    !end if

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
      write(*,*) "flag_x235, use_zfft1d_x =",flag_x235, use_zfft1d_x
      write(*,*) "flag_y235, use_zfft1d_y =",flag_y235, use_zfft1d_y
      write(*,*) "flag_z235, use_zfft1d_z =",flag_z235, use_zfft1d_z
    end if

    i = max( nnx, nny, nnz, mmx, mmy, mmz )
    allocate( a(i)      ); a=z0
    allocate( b(i)      ); b=z0
    allocate( c(i)      ); c=z0
    allocate( work(2*i) ); work=z0

    pi = acos(-1.0d0)
    wx = dcmplx( cos(pi/nnx), -sin(pi/nnx) )
    wy = dcmplx( cos(pi/nny), -sin(pi/nny) )
    wz = dcmplx( cos(pi/nnz), -sin(pi/nnz) )

    if ( .not.flag_x235 ) then
      allocate( wxa(0:mmx-1) ); wxa=z0
      allocate( wxb(0:nnx-1) ); wxb=z0
      do i = -nnx+1, nnx-1
        j = mod( i+mmx, mmx )
        wxa(j) = wx**(-i*i)
      end do
      do i = 0, nnx-1
        wxb(i) = wx**(i*i)
      end do
    end if
    if ( .not.flag_y235 ) then
      allocate( wya(0:mmy-1) ); wya=z0
      allocate( wyb(0:nny-1) ); wyb=z0
      do i = -nny+1, nny-1
        j = mod( i+mmy, mmy )
        wya(j) = wy**(-i*i)
      end do
      do i = 0, nny-1
        wyb(i) = wy**(i*i)
      end do
    end if
    if ( .not.flag_z235 ) then
      allocate( wza(0:mmz-1) ); wza=z0
      allocate( wzb(0:nnz-1) ); wzb=z0
      do i = -nnz+1, nnz-1
        j = mod( i+mmz, mmz )
        wza(j) = wz**(-i*i)
      end do
      do i = 0, nnz-1
        wzb(i) = wz**(i*i)
      end do
    end if

    if ( allocated(zwork2_ffte) ) deallocate(zwork2_ffte)
    if ( allocated(zwork1_ffte) ) deallocate(zwork1_ffte)
    allocate( zwork1_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork1_ffte=z0
    allocate( zwork2_ffte(nx_0-1:nx_1-1,ny_0-1:ny_1-1,nz_0-1:nz_1-1) ); zwork2_ffte=z0

    flag_init_done = .true.

    call write_border( 0, 'init_pzfft3dv_arb(end)' )

    contains

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

  end subroutine init_pzfft3dv_arb

  subroutine pzfft3dv_arb( f, iopt )
    implicit none
    complex(8),intent(inout) :: f(:,:,:)
    integer,intent(in) :: iopt
    integer :: ix,iy,iz
    !complex(8) :: ww

    !call write_border( 0, 'pzfft3dv_arb(start)' )
    if ( .not.flag_init_done ) call stop_program('Call "init_pzfft3dv_arb" first.')

! ---

    if ( flag_x235 ) then

      if ( use_zfft1d_x ) then

      call ZFFT1D( a, nnx, 0, work ) 
      do iz=1,nz
      do iy=1,ny
!$omp parallel workshare
        a(1:nx) = f(1:nx,iy,iz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:nx), b(1:nnx), ircnx, idisx, comm_fftx )
        call ZFFT1D( b, nnx, iopt, work ) 
!$omp parallel workshare
        f(1:nx,iy,iz) = b(nx_0:nx_1)
!$omp end parallel workshare
      end do
      end do

      else

      call PZFFT1D( a,b,c,nnx,comm_fftx,me_fftx,npux,0 )
      do iz=1,nz
      do iy=1,ny
        a(1:nx) = f(1:nx,iy,iz)
        call PZFFT1D( a,b,c,nnx,comm_fftx,me_fftx,npux,iopt )
        f(1:nx,iy,iz) = b(1:nx)
      end do
      end do

      end if

    else

      !ww=wx ; if ( iopt == 1 ) ww=conjg(ww)
      if ( iopt == 1 ) then
        wxa = conjg(wxa)
        wxb = conjg(wxb)
      end if
      call ZFFT1D( a, mmx, 0, work ) 
      do iz=1,nz
      do iy=1,ny
!$omp parallel workshare
        a(1:nx) = f(1:nx,iy,iz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:nx), b(1:nnx), ircnx, idisx, comm_fftx )
        !call ZFFT1D_arb( b(1:mmx),a(1:mmx),c(1:mmx),work(1:2*mmx),ww,mmx,nnx,nx_0-1,nx_1-1 )
        call ZFFT1D_arb_test( b(1:mmx),a(1:mmx),work(1:2*mmx),wxa,wxb,nx_0-1,nx_1-1 )
!$omp parallel workshare
        f(1:nx,iy,iz) = b(nx_0:nx_1)
!$omp end parallel workshare
      end do
      end do
      if ( iopt == 1 ) then
        f(:,:,:) = f(:,:,:) * inv_nnx
        wxa = conjg(wxa)
        wxb = conjg(wxb)
      end if

    end if

! ---

    if ( flag_y235 ) then

      if ( use_zfft1d_y ) then

      call ZFFT1D( a, nny, 0, work ) 
      do iz=1,nz
      do ix=1,nx
!$omp parallel workshare
        a(1:ny) = f(ix,1:ny,iz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:ny), b(1:nny), ircny, idisy, comm_ffty )
        call ZFFT1D( b, nny, iopt, work )
!$omp parallel workshare
        f(ix,1:ny,iz) = b(ny_0:ny_1)
!$omp end parallel workshare
      end do
      end do

      else

      call PZFFT1D( a,b,c,nny,comm_ffty,me_ffty,npuy,0 )
      do iz=1,nz
      do ix=1,nx
        a(1:ny) = f(ix,1:ny,iz)
        call PZFFT1D( a,b,c,nny,comm_ffty,me_ffty,npuy,iopt )
        f(ix,1:ny,iz) = b(1:ny)
      end do
      end do

      end if

    else

      !ww=wy ; if ( iopt == 1 ) ww=conjg(ww)
      if ( iopt == 1 ) then
        wya = conjg(wya)
        wyb = conjg(wyb)
      end if
      call ZFFT1D( a, mmy, 0, work ) 
      do iz=1,nz
      do ix=1,nx
!$omp parallel workshare
        a(1:ny) = f(ix,1:ny,iz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:ny), b(1:nny), ircny, idisy, comm_ffty )
        !call ZFFT1D_arb( b(1:mmy),a(1:mmy),c(1:mmy),work(1:2*mmy),ww,mmy,nny,ny_0-1,ny_1-1 )
        call ZFFT1D_arb_test( b(1:mmy),a(1:mmy),work(1:2*mmy),wya,wyb,ny_0-1,ny_1-1 )
!$omp parallel workshare
        f(ix,1:ny,iz) = b(ny_0:ny_1)
!$omp end parallel workshare
      end do
      end do
      if ( iopt == 1 ) then
        f(:,:,:) = f(:,:,:) * inv_nny
        wya = conjg(wya)
        wyb = conjg(wyb)
      end if

    end if

! ---

    if ( flag_z235 ) then

      if ( use_zfft1d_z ) then

      call ZFFT1D( a, nnz, 0, work ) 
      do iy=1,ny
      do ix=1,nx
!$omp parallel workshare
        a(1:nz) = f(ix,iy,1:nz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:nz), b(1:nnz), ircnz, idisz, comm_fftz )
        call ZFFT1D( b, nnz, iopt, work )
!$omp parallel workshare
        f(ix,iy,1:nz) = b(nz_0:nz_1)
!$omp end parallel workshare
      end do
      end do

      else

      call PZFFT1D( a,b,c,nnz,comm_fftz,me_fftz,npuz,0 ) 
      do iy=1,ny
      do ix=1,nx
        a(1:nz) = f(ix,iy,1:nz)
        call PZFFT1D( a,b,c,nnz,comm_fftz,me_fftz,npuz,iopt )
        f(ix,iy,1:nz) = b(1:nz)
      end do
      end do

      end if

    else

      !ww=wz ; if ( iopt == 1 ) ww=conjg(ww)
      if ( iopt == 1 ) then
        wza = conjg(wza)
        wzb = conjg(wzb)
      end if
      call ZFFT1D( a, mmz, 0, work ) 
      do iy=1,ny
      do ix=1,nx
!$omp parallel workshare
        a(1:nz) = f(ix,iy,1:nz)
!$omp end parallel workshare
        call rsdft_allgatherv( a(1:nz), b(1:nnz), ircnz, idisz, comm_fftz )
        !call ZFFT1D_arb( b(1:mmz),a(1:mmz),c(1:mmz),work(1:2*mmz),ww,mmz,nnz,nz_0-1,nz_1-1 )
        call ZFFT1D_arb_test( b(1:mmz),a(1:mmz),work(1:2*mmz),wza,wzb,nz_0-1,nz_1-1 )
!$omp parallel workshare
        f(ix,iy,1:nz) = b(nz_0:nz_1)
!$omp end parallel workshare
      end do
      end do
      if ( iopt == 1 ) then
        f(:,:,:) = f(:,:,:) * inv_nnz
        wza = conjg(wza)
        wzb = conjg(wzb)
      end if

    end if

! ---

    !call write_border( 0, 'pzfft3dv_arb(end)' )

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
!$omp parallel
!$omp workshare
    a(:) = z0
!$omp end workshare
!$omp do private(j)
    do i = -nn+1, nn-1
      j = mod( i+mm, mm )
      a(j) = ww**(-i*i)
    end do
!$omp end do
!$omp do
    do i = 0, nn-1
      b(i) = b(i) * ww**(i*i)
    end do
!$omp end do
!$omp end parallel
    if ( size(b) > nn ) b(nn:)=z0
    call ZFFT1D( a, mm,-1, work )
    call ZFFT1D( b, mm,-1, work )
!$omp parallel workshare
    c = a*b
!$omp end parallel workshare
    call ZFFT1D( c, mm, 1, work )
!$omp parallel
!$omp do 
    do i = 0, nn-1
      c(i) = c(i) * ww**(i*i)
    end do
!$omp end do
!$omp do
    do i = nn_0, nn_1
      b(i) = c(i)
    end do
!$omp end do
!$omp end parallel
  end subroutine ZFFT1D_arb

  subroutine ZFFT1D_arb_test( b, a, work, wwa, wwb, nn_0, nn_1 )
    implicit none
    complex(8),intent(inout) :: b(0:)
    complex(8),intent(inout) :: a(0:)
    complex(8),intent(inout) :: work(0:)
    complex(8),intent(in) :: wwa(0:), wwb(0:)
    integer,intent(in) :: nn_0,nn_1
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    integer :: i,mm,nn
    mm = size(wwa)
    nn = size(wwb)
!$omp parallel
!$omp workshare
    a(:) = wwa(:)
!$omp end workshare
!$omp do
    do i = 0, nn-1
      b(i) = b(i) * wwb(i)
    end do
!$omp end do
!$omp end parallel
    if ( size(b) > nn ) b(nn:)=z0
    call ZFFT1D( a, mm,-1, work ) 
    call ZFFT1D( b, mm,-1, work ) 
!$omp parallel workshare
    a = a*b
!$omp end parallel workshare
    call ZFFT1D( a, mm, 1, work ) 
!$omp parallel do
    do i = nn_0, nn_1
      b(i) = a(i) * wwb(i)
    end do
!$omp end parallel do
  end subroutine ZFFT1D_arb_test


end module pzfft3dv_arb_module
