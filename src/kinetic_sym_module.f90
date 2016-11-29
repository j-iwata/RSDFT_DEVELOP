MODULE kinetic_sym_module

!$  use omp_lib
  use rgrid_module
  use parallel_module, only: ir_grid, id_grid, comm_grid
  use rsdft_mpi_module
  use grid_module, only: construct_map_3d_to_1d_grid
  use bz_module, only: kbb
  use enlarge_array_module

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_sym
  PUBLIC :: init2

  logical :: init_flag=.true.
  integer,allocatable :: LLL(:,:,:)
  integer :: n_fd_points
  integer,allocatable :: n_fd_pt(:)
  integer,allocatable :: fd_neighbor_pt(:,:,:)
#ifdef _DRSDFT_
  real(8),allocatable :: fd_coef(:,:)
#else
  complex(8),allocatable :: fd_coef(:,:)
#endif

#ifdef _DRSDFT_
    real(8),allocatable :: work(:)
#else
    complex(8),allocatable :: work(:)
#endif

CONTAINS


  SUBROUTINE init

    implicit none
    integer,parameter :: u=10
    integer :: i,myrank
    include 'mpif.h'

    if ( .not.init_flag ) return
    init_flag=.false.

    call write_border( 0, "init in kinetic_sym_module(start)" )

    call construct_map_3d_to_1d_grid( Ngrid, Igrid, comm_grid, LLL )
    allocate( work(Ngrid(0)) ) ; work=(0.0d0,0.0d0)

    if ( allocated(fd_coef) ) return

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, i )
    if ( myrank == 0 ) then
       open(u,file="kin_coef.dat",status="old")
       read(u,*) n_fd_points
       allocate( fd_coef(n_fd_points,1) ) ; fd_coef=0.0d0
       allocate( fd_neighbor_pt(3,n_fd_points,1) ) ; fd_neighbor_pt=0
       do i=1,n_fd_points
          read(u,*) fd_neighbor_pt(:,i,1), fd_coef(i,1)
       end do
       close(u)
    end if
    call MPI_BCAST( n_fd_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
    if ( myrank /= 0 ) then
       allocate( fd_coef(n_fd_points,1) ) ; fd_coef=0.0d0
       allocate( fd_neighbor_pt(3,n_fd_points,1) ) ; fd_neighbor_pt=0
    end if
    call MPI_BCAST( fd_neighbor_pt, size(fd_neighbor_pt) &
                  , MPI_INTEGER, 0, MPI_COMM_WORLD, i )
#ifdef _DRSDFT_
    call MPI_BCAST( fd_coef, size(fd_coef), MPI_REAL8, 0, MPI_COMM_WORLD, i )
#else
    call MPI_BCAST( fd_coef, size(fd_coef), MPI_COMPLEX16, 0, MPI_COMM_WORLD, i )
#endif
    fd_coef=fd_coef/( Hgrid(1)**2 )

    call write_border( 0, "init in kinetic_sym_module(end)" )

  END SUBROUTINE init


  SUBROUTINE init2( n, nk, k, tmp )
    implicit none
    integer,intent(IN) :: n, nk, k
#ifdef _DRSDFT_
    real(8),intent(IN) :: tmp(-n:n,-n:n,-n:n)
#else
    complex(8),intent(IN) :: tmp(-n:n,-n:n,-n:n)
#endif
    integer :: i1,i2,i3,m
    call write_border( 0, " init2(start)" )
    m = count( abs(tmp) > 1.d-10 )
    if ( n_fd_points == 0 ) then
       n_fd_points=m
       allocate( n_fd_pt(nk)                      ) ; n_fd_pt=0
       allocate( fd_neighbor_pt(3,n_fd_points,nk) ) ; fd_neighbor_pt=0
       allocate( fd_coef(n_fd_points,nk)          ) ; fd_coef=(0.0d0,0.0d0)
    else if ( n_fd_points < m ) then
       n_fd_points=m
       call enlarge_array( fd_neighbor_pt, 3, m, nk )
       call enlarge_array( fd_coef, m, nk )
    end if
    m=0
    do i3=-n,n
    do i2=-n,n
    do i1=-n,n
       if ( abs(tmp(i1,i2,i3)) > 1.d-10 ) then
          m=m+1
          fd_neighbor_pt(1,m,k) = i1
          fd_neighbor_pt(2,m,k) = i2
          fd_neighbor_pt(3,m,k) = i3
          fd_coef(m,k) = tmp(i1,i2,i3)
       end if
    end do
    end do
    end do
    n_fd_pt(k)=m
  END SUBROUTINE init2


  SUBROUTINE op_kinetic_sym( k, tpsi, htpsi )
    implicit none
    integer,intent(IN) :: k
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)    ::  tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8) :: bp
    real(8) :: kr,r(3)
    real(8),parameter :: pi2=6.283185307179586d0
#endif
    integer :: n,i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,n1,n2

    if ( init_flag ) call init

    n1 = Igrid(1,0)
    n2 = Igrid(2,0)

    do n=1,size(tpsi,2)

       work(n1:n2)=tpsi(:,n)
       call rsdft_allgatherv( work(n1:n2), work, ir_grid, id_grid, comm_grid )

!#ifndef _DRSDFT_
!       do k3=0,Ngrid(3)-1
!       do k2=0,Ngrid(2)-1
!       do k1=0,Ngrid(1)-1
!          r(1)=dble(k1)/dble(Ngrid(1))
!          r(2)=dble(k2)/dble(Ngrid(2))
!          r(3)=dble(k3)/dble(Ngrid(3))
!          kr=pi2*( r(1)*kbb(1,k) + r(2)*kbb(2,k) + r(3)*kbb(3,k) )
!          bp=dcmplx( cos(kr), sin(kr) )
!          work( LLL(k1,k2,k3) ) = work( LLL(k1,k2,k3) )*bp
!       end do
!       end do
!       end do
!#endif

       i=0
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1

          do j=1,n_fd_pt(k)

             j1 = i1 + fd_neighbor_pt(1,j,k)
             j2 = i2 + fd_neighbor_pt(2,j,k)
             j3 = i3 + fd_neighbor_pt(3,j,k)

             k1 = mod( j1+Ngrid(1),Ngrid(1) )
             k2 = mod( j2+Ngrid(2),Ngrid(2) )
             k3 = mod( j3+Ngrid(3),Ngrid(3) )

!#ifdef _DRSDFT_
             htpsi(i,n) = htpsi(i,n) + fd_coef(j,k)*work( LLL(k1,k2,k3) )
!#else
!             call calc_bloch_phase( j1,j2,j3, kbb(:,k), bp )
!             htpsi(i,n) = htpsi(i,n) + bp*fd_coef(j,k)*work( LLL(k1,k2,k3) )
!#endif



          end do ! j

       end do ! i1
       end do ! i2
       end do ! i3

!#ifndef _DRSDFT_
!       i=0
!       do i3=Igrid(1,3),Igrid(2,3)
!       do i2=Igrid(1,2),Igrid(2,2)
!       do i1=Igrid(1,1),Igrid(2,1)
!          i=i+1
!          r(1)=dble(i1)/dble(Ngrid(1))
!          r(2)=dble(i2)/dble(Ngrid(2))
!          r(3)=dble(i3)/dble(Ngrid(3))
!          kr=pi2*( r(1)*kbb(1,k) + r(2)*kbb(2,k) + r(3)*kbb(3,k) )
!          bp=dcmplx( cos(kr), -sin(kr) )
!          htpsi(i,n) = htpsi(i,n)*bp
!       end do
!       end do
!       end do
!#endif

    end do ! n

  END SUBROUTINE op_kinetic_sym


  SUBROUTINE calc_bloch_phase( i1,i2,i3, k, bp )
    implicit none
    integer,intent(IN) :: i1,i2,i3
    real(8),intent(IN) :: k(3)
    complex(8),intent(OUT) :: bp
    integer :: n(3)
    real(8) :: kR
    real(8),parameter :: pi2=6.283185307179586d0
    if ( i1 < 0 ) then
       n(1)=-1
    else if ( Ngrid(1) <= i1 ) then
       n(1)= 1
    end if
    if ( i2 < 0 ) then
       n(2)=-1
    else if ( Ngrid(2) <= i2 ) then
       n(2)= 1
    end if
    if ( i3 < 0 ) then
       n(3)=-1
    else if ( Ngrid(3) <= i3 ) then
       n(3)= 1
    end if
    kR=pi2*( k(1)*n(1) + k(2)*n(2) + k(3)*n(3) )
    bp=dcmplx( cos(kR), sin(kR) )
  END SUBROUTINE calc_bloch_phase

END MODULE kinetic_sym_module
