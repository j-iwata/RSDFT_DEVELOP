MODULE kinetic_sym_module

!$  use omp_lib
  use rgrid_module
  use parallel_module, only: ir_grid, id_grid, comm_grid,myrank
  use rsdft_mpi_module
  use grid_module, only: construct_map_3d_to_1d_grid
  use bz_module, only: kbb
  use enlarge_array_module
  use watch_module, only: watchb_omp, time_kine, time_bcfd
  use omp_variables

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
    real(8),allocatable :: work(:),work3(:,:,:)
#else
    complex(8),allocatable :: work(:),work3(:,:,:)
#endif

    integer :: ba1,bb1,ba2,bb2,ba3,bb3

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

    ba1=minval( fd_neighbor_pt(1,:,:) )+Igrid(1,1)
    bb1=maxval( fd_neighbor_pt(1,:,:) )+Igrid(2,1)
    ba2=minval( fd_neighbor_pt(2,:,:) )+Igrid(1,2)
    bb2=maxval( fd_neighbor_pt(2,:,:) )+Igrid(2,2)
    ba3=minval( fd_neighbor_pt(3,:,:) )+Igrid(1,3)
    bb3=maxval( fd_neighbor_pt(3,:,:) )+Igrid(2,3)
    allocate( work3(ba1:bb1,ba2:bb2,ba3:bb3) ) ; work3=(0.0d0,0.0d0)

    if ( allocated(fd_coef) ) then
       call write_border( 0, "init in kinetic_sym_module(return)" )
       return
    end if

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
    call write_border( 0, " init2 in kinetic_sym_module(start)" )
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
       if ( abs(tmp(i1,i2,i3)) > 1.d-8 ) then
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
    if ( init_flag ) call init
    call write_border( 0, " init2 in kinetic_sym_module(end)" )
  END SUBROUTINE init2


  SUBROUTINE op_kinetic_sym_a( tpsi, htpsi, k_in )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),intent(IN)    ::  tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: bp
    real(8) :: kr,r(3)
    real(8),parameter :: pi2=6.283185307179586d0
#endif
    integer,optional,intent(IN) :: k_in
    integer :: n,i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,n1,n2,k
    real(8) :: ttmp(2)

!    if ( init_flag ) call init

    k  = 1 ; if ( present(k_in) ) k=k_in
    n1 = Igrid(1,0)
    n2 = Igrid(2,0)

    do n=1,size(tpsi,2)

       !call watchb_omp( ttmp )

       work(n1:n2)=tpsi(:,n)
       call rsdft_allgatherv( work(n1:n2), work, ir_grid, id_grid, comm_grid )

       !call watchb_omp( ttmp, time_kine(1,3) )

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

       !call watchb_omp( ttmp, time_kine(1,1) )

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

  END SUBROUTINE op_kinetic_sym_a


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


  SUBROUTINE op_kinetic_sym( tpsi, htpsi, k_in )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8),parameter :: zero=0.0d0
    real(8) :: tmp
#else
    complex(8),intent(IN)    ::  tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: tmp
#endif
    integer,optional,intent(IN) :: k_in
    integer :: n,i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,n1,n2,k
    integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp,n1_omp,n2_omp
    integer :: ib1_omp,ib2_omp,nb_omp
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: ttmp(2)

    if ( init_flag ) call init

    k  = 1 ; if ( present(k_in) ) k=k_in
    n1 = Igrid(1,0)
    n2 = Igrid(2,0)

! ---

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    n=0
!$  n=omp_get_thread_num()

    n1_omp = Igrid_omp(1,0,n)
    n2_omp = Igrid_omp(2,0,n)

    a1b_omp = Igrid_omp(1,1,n)
    b1b_omp = Igrid_omp(2,1,n)
    a2b_omp = Igrid_omp(1,2,n)
    b2b_omp = Igrid_omp(2,2,n)
    a3b_omp = Igrid_omp(1,3,n)
    b3b_omp = Igrid_omp(2,3,n)

! ---

    do n=1,size(tpsi,2)

       !call watchb_omp( ttmp )

!$OMP workshare
       work(n1:n2)=tpsi(:,n)
!$OMP end workshare
!$OMP single
       call rsdft_allgatherv( work(n1:n2), work, ir_grid, id_grid, comm_grid )
!$OMP end single

!$OMP do private( i1,i2,i3,k1,k2,k3 )
       do i3=ba3,bb3
       do i2=ba2,bb2
       do i1=ba1,bb1
          k1 = mod( i1+Ngrid(1),Ngrid(1) )
          k2 = mod( i2+Ngrid(2),Ngrid(2) )
          k3 = mod( i3+Ngrid(3),Ngrid(3) )
          work3(i1,i2,i3)=work( LLL(k1,k2,k3) )
       end do
       end do
       end do
!$OMP end do

       !call watchb_omp( ttmp, time_kine(1,3) )

       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          i=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12

          tmp=zero
          do j=1,n_fd_pt(k)
             j1 = i1 + fd_neighbor_pt(1,j,k)
             j2 = i2 + fd_neighbor_pt(2,j,k)
             j3 = i3 + fd_neighbor_pt(3,j,k)
             tmp = tmp + fd_coef(j,k)*work3(j1,j2,j3)
          end do ! j

          htpsi(i,n) = htpsi(i,n) + tmp

       end do ! i1
       end do ! i2
       end do ! i3

       !call watchb_omp( ttmp, time_kine(1,1) )

    end do ! n

  END SUBROUTINE op_kinetic_sym


END MODULE kinetic_sym_module
