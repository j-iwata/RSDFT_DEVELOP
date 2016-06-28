MODULE fft_module

  use grid_module
  use rsdft_fft_module
  use fftw_module
  use bb_module
  use hsort_module

  implicit none

  PRIVATE
  PUBLIC :: init_fft
  PUBLIC :: finalize_fft
  PUBLIC :: d1_to_z3_fft, z3_to_d1_fft
  PUBLIC :: z1_to_z3_fft, z3_to_z1_fft
  PUBLIC :: forward_fft, backward_fft
  PUBLIC :: grid_util_fft

  type(grid) :: rgrid
  integer :: ML,ML1,ML2,ML3
  integer,allocatable :: LLL(:,:,:)
!  integer :: ifacx(30),ifacy(30),ifacz(30)
!  integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
!  complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  logical :: keep_LLL=.false.

CONTAINS


  SUBROUTINE init_fft( keep_flag )
    implicit none
    logical,optional,intent(IN) :: keep_flag
    call get_range_rgrid( rgrid )
    if ( .not.keep_LLL ) call get_map_3d_to_1d_grid( rgrid, LLL )
    ML  = rgrid%g1%size_global
    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global
    if ( present(keep_flag) ) then
       keep_LLL = keep_flag
       return
    end if
!    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
!    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )
!    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
!               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
  END SUBROUTINE init_fft


  SUBROUTINE finalize_fft
    implicit none
!    if ( allocated(wsavez) ) deallocate( wsavez,wsavey,wsavex )
!    if ( allocated(lz2)    ) deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    if ( .not.keep_LLL ) then
       if ( allocated(LLL) ) deallocate( LLL )
    end if
  END SUBROUTINE finalize_fft


  SUBROUTINE d1_to_z3_fft( d1, z3 )
    implicit none
    real(8),intent(IN) :: d1(:)
    complex(8),allocatable :: z3(:,:,:) 
    real(8),allocatable :: work(:)
    integer :: i,i1,i2,i3
    allocate( work(ML) ) ; work=0.0d0
    call mpi_allgatherv_grid( d1, work )
    if ( .not.allocated(z3) ) then
       allocate( z3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z3=zero
    end if
    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       i=LLL(i1,i2,i3)
       z3(i1,i2,i3) = work(i)
    end do
    end do
    end do
    deallocate( work )
  END SUBROUTINE d1_to_z3_fft


  SUBROUTINE z3_to_d1_fft( z3, d1 )
    implicit none
    complex(8),intent(IN) :: z3(0:,0:,0:)
    real(8),intent(OUT) :: d1(:)
    integer :: i1,i2,i3,i
    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       d1(i) = dble( z3(i1,i2,i3) )
    end do
    end do
    end do
  END SUBROUTINE z3_to_d1_fft


  SUBROUTINE z1_to_z3_fft( z1, z3 )
    implicit none
    complex(8),intent(IN)  :: z1(:)
    complex(8),allocatable :: z3(:,:,:)
    integer :: i1,i2,i3,i
    complex(8),allocatable :: work(:)
    allocate( work(ML) ) ; work=0.0d0
    call zmpi_allgatherv_grid( z1, work )
    if ( .not.allocated(z3) ) then
       allocate( z3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z3=zero
    end if
    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       i=LLL(i1,i2,i3)
       z3(i1,i2,i3) = work(i)
    end do
    end do
    end do
    deallocate( work )
  END SUBROUTINE z1_to_z3_fft


  SUBROUTINE z3_to_z1_fft( z3, z1 )
    implicit none
    complex(8),intent(IN)  :: z3(0:,0:,0:)
    complex(8),intent(OUT) :: z1(:)
    integer :: i1,i2,i3,i
    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       z1(i) = z3(i1,i2,i3)
    end do
    end do
    end do
  END SUBROUTINE z3_to_z1_fft


  SUBROUTINE forward_fft( z3, w3 )
    implicit none
    complex(8),intent(INOUT) :: z3(0:,0:,0:)
    complex(8),allocatable   :: w3(:,:,:)
#ifdef _FFTW_
    call forward_fftw( z3 )
    return
#endif
    if ( .not.allocated(w3) ) then
       allocate( w3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; w3=zero
    end if
    call rsdft_fft3d( z3,  1 ) ; return
!    call fft3fx(ML1,ML2,ML3,ML,z3,w3,wsavex,wsavey,wsavez &
!               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
  END SUBROUTINE forward_fft


  SUBROUTINE backward_fft( z3, w3 )
    implicit none
    complex(8),intent(INOUT) :: z3(:,:,:)
    complex(8),allocatable   :: w3(:,:,:)
#ifdef _FFTW_
    call backward_fftw( z3 )
    return
#endif
    if ( .not.allocated(w3) ) then
       allocate( w3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; w3=zero
    end if
    call rsdft_fft3d( z3, -1 ) ; return
!    call fft3bx(ML1,ML2,ML3,ML,z3,w3,wsavex,wsavey,wsavez &
!               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
  END SUBROUTINE backward_fft


  SUBROUTINE grid_util_fft( gg, indx )
    implicit none
    real(8),allocatable,intent(INOUT) :: gg(:)
    integer,allocatable,intent(INOUT) :: indx(:)
    integer :: i1,i2,i3,j1,j2,j3,jj
    real(8) :: gx,gy,gz
    if ( .not.allocated(gg) ) allocate( gg(ML) )
    if ( .not.allocated(indx) ) allocate( indx(ML) )
    gg=0.0d0
    do i3=-(ML3-1)/2,ML3/2
    do i2=-(ML2-1)/2,ML2/2
    do i1=-(ML1-1)/2,ML1/2
       j1=mod(i1+ML1,ML1)
       j2=mod(i2+ML2,ML2)
       j3=mod(i3+ML3,ML3)
       jj=1+j1+ML1*j2+ML1*ML2*j3
       gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
       gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
       gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
       gg(jj)=gx*gx+gy*gy+gz*gz
    end do
    end do
    end do
    call indexx( ML, gg, indx )
!    write(*,*) minval(gg),maxval(gg)
!    write(*,*) count(gg/=0.0d0),size(gg)
!    do jj=1,ML
!       write(*,*) jj,gg(jj),gg(indx(jj))
!    end do
!    stop
  END SUBROUTINE grid_util_fft


END MODULE fft_module
