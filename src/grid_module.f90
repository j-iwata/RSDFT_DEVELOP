MODULE grid_module

  use rgrid_variables
  use parallel_module
  use basic_type_factory
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: grid, get_range_rgrid
  PUBLIC :: construct_map_3d_to_1d_grid
  PUBLIC :: construct_map_1d_to_3d_grid
  PUBLIC :: get_map_3d_to_1d_grid
  PUBLIC :: mpi_allgatherv_grid, zmpi_allgatherv_grid
  PUBLIC :: inner_product_grid
  PUBLIC :: convert_1d_to_3d_grid
  PUBLIC :: get_info_rs_grid
  PUBLIC :: z_convert_1d_to_3d_grid
  PUBLIC :: z_convert_3d_to_1d_grid
  PUBLIC :: grid_info

  type grid
     type( ArrayRange1D ) :: g1
     type( ArrayRange3D ) :: g3
     real(8) :: spacing(3)
     real(8) :: VolumeElement
     integer :: comm
  end type grid

  type grid_info
     integer :: head(0:3)
     integer :: tail(0:3)
     integer :: size(0:3)
     real(8) :: spacing(3)
     real(8) :: dV
     real(8) :: lattice(3,3)
     real(8) :: Rsize, Zsize
     integer :: comm
     integer :: Norder_FD
     integer,allocatable :: LL(:,:)
     integer,allocatable :: KK(:,:)
  end type grid_info

CONTAINS

  SUBROUTINE get_range_rgrid( rgrid )
    implicit none
    type(grid) :: rgrid
    rgrid%g1%head = Igrid(1,0)
    rgrid%g1%tail = Igrid(2,0)
    rgrid%g1%size = Igrid(2,0)-Igrid(1,0)+1
    rgrid%g1%head_global = 1
    rgrid%g1%tail_global = Ngrid(0)
    rgrid%g1%size_global = Ngrid(0)
    rgrid%g3%x%head = Igrid(1,1)
    rgrid%g3%x%tail = Igrid(2,1)
    rgrid%g3%x%size = Igrid(2,1)-Igrid(1,1)+1
    rgrid%g3%y%head = Igrid(1,2)
    rgrid%g3%y%tail = Igrid(2,2)
    rgrid%g3%y%size = Igrid(2,2)-Igrid(1,2)+1
    rgrid%g3%z%head = Igrid(1,3)
    rgrid%g3%z%tail = Igrid(2,3)
    rgrid%g3%z%size = Igrid(2,3)-Igrid(1,3)+1
    rgrid%g3%x%head_global = 0
    rgrid%g3%x%tail_global = Ngrid(1)-1
    rgrid%g3%x%size_global = Ngrid(1)
    rgrid%g3%y%head_global = 0
    rgrid%g3%y%tail_global = Ngrid(2)-1
    rgrid%g3%y%size_global = Ngrid(2)
    rgrid%g3%z%head_global = 0
    rgrid%g3%z%tail_global = Ngrid(3)-1
    rgrid%g3%z%size_global = Ngrid(3)
    rgrid%spacing(1:3) = Hgrid(1:3)
    rgrid%VolumeElement = dV
    rgrid%comm = comm_grid
  END SUBROUTINE get_range_rgrid


  subroutine construct_map_3d_to_1d_grid( LLL )
    use parallel_module, only: pinfo_grid
    implicit none
    integer,allocatable,intent(inout) :: LLL(:,:,:)
    integer :: i,i1,i2,i3,irank,n1,n2,n3
    if ( .not.allocated(LLL) ) then
      n1 = sum( pinfo_grid(2,:) )
      n2 = sum( pinfo_grid(4,:) )
      n3 = sum( pinfo_grid(6,:) )
      allocate( LLL(0:n1-1,0:n2-1,0:n3-1) )
    end if
    LLL=0
    do irank = 0, size(pinfo_grid,2)-1
      i = pinfo_grid(7,irank)
      do i3 = pinfo_grid(5,irank), pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
      do i2 = pinfo_grid(3,irank), pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
      do i1 = pinfo_grid(1,irank), pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
        i=i+1
        LLL(i1,i2,i3)=i
      end do
      end do
      end do
    end do
  end subroutine construct_map_3d_to_1d_grid


  SUBROUTINE get_map_3d_to_1d_grid( rgrid, LLL )
    implicit none
    type(grid),intent(IN) :: rgrid
    integer,allocatable,intent(INOUT) :: LLL(:,:,:)
    integer :: i,i1,i2,i3,m1,m2,m3,n1,n2,n3
    if ( .not.allocated(LLL) ) then
       m1=rgrid%g3%x%head_global
       n1=rgrid%g3%x%tail_global
       m2=rgrid%g3%y%head_global
       n2=rgrid%g3%y%tail_global
       m3=rgrid%g3%z%head_global
       n3=rgrid%g3%z%tail_global
       allocate( LLL(m1:n1,m2:n2,m3:n3) )
    end if
    LLL=0
    i=rgrid%g1%head-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    call rsdft_allreduce_sum( LLL, rgrid%comm )
  END SUBROUTINE get_map_3d_to_1d_grid


  subroutine construct_map_1d_to_3d_grid( LL )
    use parallel_module, only: pinfo_grid
    implicit none
    integer,allocatable,intent(inout) :: LL(:,:)
    integer :: i,i1,i2,i3,N,irank
    if ( .not.allocated(LL) ) then
      N = sum( pinfo_grid(8,:) )
      allocate( LL(3,N) )
    end if
    LL=0
    do irank = 0, size(pinfo_grid,2)-1
      i = pinfo_grid(7,irank)
      do i3 = pinfo_grid(5,irank), pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
      do i2 = pinfo_grid(3,irank), pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
      do i1 = pinfo_grid(1,irank), pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
        i = i + 1
        LL(1:3,i) = (/ i1, i2, i3 /)
      end do
      end do
      end do
    end do
  end subroutine construct_map_1d_to_3d_grid

  ! SUBROUTINE get_map_3d_to_1d( LLL )
  !    implicit none
  !    integer,allocatable,intent(OUT) :: LLL(:,:,:)
  !    integer :: i,i1,i2,i3,ierr
  !    allocate( LLL(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) ) ; LLL=0
  !    i=Igrid(1,0)-1
  !    do i3=Igrid(1,3),Igrid(2,3)
  !    do i2=Igrid(1,2),Igrid(2,2)
  !    do i1=Igrid(1,1),Igrid(2,1)
  !       i=i+1
  !       LLL(i1,i2,i3)=i
  !    end do
  !    end do
  !    end do
  !    call rsdft_allreduce_sum( LLL, comm_grid )
  ! END SUBROUTINE get_map_3d_to_1d

  SUBROUTINE mpi_allgatherv_grid( f, g )
    implicit none
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: g(:)
    integer :: ierr
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir_grid, id_grid, &
                         MPI_REAL8, comm_grid, ierr )
  END SUBROUTINE mpi_allgatherv_grid

  SUBROUTINE zmpi_allgatherv_grid( zf, zg )
    implicit none
    complex(8),intent(IN)  :: zf(:)
    complex(8),intent(OUT) :: zg(:)
    integer :: ierr
    call MPI_ALLGATHERV( zf, size(zf), RSDFT_MPI_COMPLEX16, zg, ir_grid, &
                         id_grid, RSDFT_MPI_COMPLEX16, comm_grid, ierr )
  END SUBROUTINE zmpi_allgatherv_grid


  SUBROUTINE inner_product_grid( f, g, c, fgc )
    implicit none
    real(8),intent(IN) :: f(:), g(:), c
    real(8),intent(OUT) :: fgc
    real(8) :: fgc_tmp
    integer :: i
    fgc_tmp=sum(f*g)*c
    call MPI_ALLREDUCE( fgc_tmp, fgc, 1, MPI_REAL8, MPI_SUM, comm_grid, i )
  END SUBROUTINE inner_product_grid


  SUBROUTINE convert_1d_to_3d_grid( i, iii )
    implicit none
    integer,intent(IN)  :: i
    integer,intent(OUT) :: iii(3)
    integer :: i1,i2,i3,j,jjj(3)
    jjj=0
    if ( Igrid(1,0) <= i .and. i <= Igrid(2,0) ) then
       j=Igrid(1,0)-1
       loop_i3 : do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          j=j+1
          if ( j == i ) then
             jjj(:)=(/ i1,i2,i3 /)
             exit loop_i3
          end if
       end do
       end do
       end do loop_i3
    end if
    call MPI_ALLREDUCE( jjj, iii, 3, MPI_INTEGER, MPI_SUM, comm_grid, j ) 
  END SUBROUTINE convert_1d_to_3d_grid


  SUBROUTINE get_info_rs_grid( Ngrid_rs, Igrid_rs, dV_rs, comm_rs )
    implicit none
    integer,optional,intent(OUT) :: Ngrid_rs(0:3), Igrid_rs(2,0:3)
    real(8),optional,intent(OUT) :: dV_rs
    integer,optional,intent(OUT) :: comm_rs
    if ( present(Ngrid_rs) ) Ngrid_rs = Ngrid
    if ( present(Igrid_rs) ) Igrid_rs = Igrid
    if ( present(dV_rs)    ) dV_rs    = dV
    if ( present(comm_rs)  ) comm_rs  = comm_grid
  END SUBROUTINE get_info_rs_grid


  SUBROUTINE z_convert_1d_to_3d_grid( lat, f1, f3 )
    implicit none
    integer,intent(IN) :: lat(:,:)
    complex(8),intent(IN) :: f1(:)
    complex(8),intent(OUT) :: f3(0:,0:,0:)
    integer :: i
    do i=1,size(f1)
       f3( lat(1,i), lat(2,i), lat(3,i) ) = f1(i)
    end do
  END SUBROUTINE z_convert_1d_to_3d_grid


  SUBROUTINE z_convert_3d_to_1d_grid( lat, f3, f1 )
    implicit none
    integer,intent(IN) :: lat(:,:)
    complex(8),intent(IN) :: f3(0:,0:,0:)
    complex(8),intent(OUT) :: f1(:)
    integer :: i,j
    do i=1,size(f1)
       j=i-1+Igrid(1,0)
       f1(i) = f3( lat(1,j), lat(2,j), lat(3,j) )
    end do
  END SUBROUTINE z_convert_3d_to_1d_grid


END MODULE grid_module
