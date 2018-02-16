MODULE noncollinear_module

  use rgrid_module, only: dV, Ngrid
  use xc_module, only: calc_xc, DCxc, Exc
  use parallel_module, only: comm_grid,comm_band,comm_bzsm,ir_grid,id_grid &
                            ,RSDFT_MPI_COMPLEX16,myrank_g
  use mixing_dm_module
  use rsdft_mpi_module
  use watch_module
  use io_ctrl_parameters
  use grid_module, only: z_convert_1d_to_3d_grid,z_convert_3d_to_1d_grid

  implicit none

  PRIVATE
  PUBLIC :: init_noncollinear
  PUBLIC :: calc_xc_noncollinear
  PUBLIC :: op_xc_noncollinear
  PUBLIC :: io_write_noncollinear, io_read_noncollinear

  logical,PUBLIC :: flag_noncollinear=.false.

  complex(8),allocatable :: den_mat(:,:,:)
  complex(8),allocatable :: vxc_mat(:,:,:)
  complex(8),allocatable :: old_mat(:,:,:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)

CONTAINS


  SUBROUTINE init_noncollinear( rho, vxc )
    implicit none
    real(8),intent(IN) :: rho(:,:), vxc(:,:)
    integer :: ml,i
    ml=size(rho,1)
    if ( .not.allocated(den_mat) ) call allocate_array( ml )
    do i=1,2
       den_mat(:,i,i)=rho(:,i)
    end do
    do i=1,2
       vxc_mat(:,i,i)=vxc(:,i)
    end do
    call calc_dc_xc( vxc_mat, den_mat, DCxc )
  END SUBROUTINE init_noncollinear


  SUBROUTINE calc_xc_noncollinear( psi, occ, rho_out, vxc_out )

    implicit none
    complex(8),intent(IN) :: psi(:,:,:,:)
    real(8),intent(IN) :: occ(:,:)
    real(8),optional,intent(OUT) :: rho_out(:,:)
    real(8),optional,intent(OUT) :: vxc_out(:,:)
    integer :: ml,i,a,b
    real(8) :: phi,theta,vxc_0,vxc_1
    real(8),allocatable :: rot_ang(:,:),rho(:,:),vxc(:,:)

    if ( .not.flag_noncollinear ) return

    call write_border( 1, "calc_xc_noncolinear(start)" )

    ml=size(psi,1)
    if ( .not.allocated(den_mat) ) call allocate_array( ml )

    allocate( rot_ang(ml,2) ) ; rot_ang=0.0d0
    allocate( rho(ml,2) ) ; rho=0.0d0
    allocate( vxc(ml,2) ) ; vxc=0.0d0

    call calc_dm_noncollinear( psi, occ, den_mat )

! ---

    call init_mixing_dm( ml, 4, comm_grid, den_mat )
    call pulay_mixing_dm( ml, 4, den_mat )

! ---

    a=1
    b=2
    do i=1,ml

       phi = -atan( aimag(den_mat(i,a,b))/real(den_mat(i,a,b)) )
       theta = atan( 2.0d0*( real(den_mat(i,a,b))*cos(phi) &
                           -aimag(den_mat(i,a,b))*sin(phi) ) &
                     /real( den_mat(i,a,a)-den_mat(i,b,b) ) )

       rho(i,1) = 0.5d0*real( den_mat(i,a,a)+den_mat(i,b,b) ) &
                + 0.5d0*real( den_mat(i,a,a)-den_mat(i,b,b) )*cos(theta) &
                + (  real(den_mat(i,a,b))*cos(phi) &
                   -aimag(den_mat(i,a,b))*sin(phi) )*sin(theta)

       rho(i,2) = 0.5d0*real( den_mat(i,a,a)+den_mat(i,b,b) ) &
                - 0.5d0*real( den_mat(i,a,a)-den_mat(i,b,b) )*cos(theta) &
                - (  real(den_mat(i,a,b))*cos(phi) &
                   -aimag(den_mat(i,a,b))*sin(phi) )*sin(theta)

       rot_ang(i,1) = phi
       rot_ang(i,2) = theta

    end do

    call calc_xc( rho, vxc, Exc )

    do i=1,ml

       phi = rot_ang(i,1)
       theta = rot_ang(i,2)

       vxc_0 = 0.5d0*( vxc(i,1) + vxc(i,2) )
       vxc_1 = 0.5d0*( vxc(i,1) - vxc(i,2) )

       vxc_mat(i,1,1) = vxc_0 + vxc_1*cos(theta)
       vxc_mat(i,2,1) = vxc_1*dcmplx( cos(phi), sin(phi) )*sin(theta)
       vxc_mat(i,1,2) = vxc_1*dcmplx( cos(phi),-sin(phi) )*sin(theta)
       vxc_mat(i,2,2) = vxc_0 - vxc_1*cos(theta)

    end do

    call calc_dc_xc( vxc_mat, den_mat, DCxc )

    if ( present(rho_out) ) rho_out=rho
    if ( present(vxc_out) ) then
       vxc_out=0.0d0
!       vxc_out(:,1)   = vxc_mat(:,1,1)
!       vxc_mat(:,1,1) = zero
!       vxc_out(:,2)   = vxc_mat(:,2,2)
!       vxc_mat(:,2,2) = zero
    end if

    deallocate( vxc )
    deallocate( rho )
    deallocate( rot_ang )

    call write_border( 1, "calc_xc_noncolinear(end)" )

  END SUBROUTINE calc_xc_noncollinear


  SUBROUTINE calc_dm_noncollinear( psi, occ, dm )
    implicit none
    complex(8),intent(IN) :: psi(:,:,:,:)
    real(8),intent(IN) :: occ(:,:)
    complex(8),intent(OUT) :: dm(:,:,:)
    integer :: n,k,i,j,mb,mk

    mb = size(psi,2)
    mk = size(psi,3)

    dm(:,:,:)=0.0d0

    do j=1,size(dm,3)
    do i=1,size(dm,2)

       do k=1,mk
       do n=1,mb

          dm(:,i,j) = dm(:,i,j) + occ(n,k)*conjg( psi(:,n,k,i) )*psi(:,n,k,j)

       end do
       end do

    end do
    end do

    do j=1,size(dm,3)
    do i=1,size(dm,2)
       call rsdft_allreduce_sum( dm(:,i,j), comm_bzsm )
    end do
    end do

  END SUBROUTINE calc_dm_noncollinear


  SUBROUTINE calc_dc_xc( vxc, den, dc_xc )
    implicit none
    complex(8),intent(IN) :: vxc(:,:,:), den(:,:,:)
    real(8),intent(OUT) :: dc_xc
    complex(8) :: ztmp,ztmp0,zmat(2,2)
    integer :: i,j
    include 'mpif.h'
    do j=1,2
    do i=1,2
       zmat(i,j) = sum( vxc(:,i,1)*den(:,1,j)+vxc(:,i,2)*den(:,2,j) )*dV
    end do
    end do
    ztmp0=zmat(1,1)+zmat(2,2)
    call MPI_ALLREDUCE( ztmp0,ztmp,1,MPI_COMPLEX16, MPI_SUM, comm_grid, i )
    dc_xc = ztmp
  END SUBROUTINE calc_dc_xc


  SUBROUTINE op_xc_noncollinear( tpsi, hpsi )
    implicit none
    complex(8),intent(IN) :: tpsi(:,:,:)
    complex(8),intent(INOUT) :: hpsi(:,:,:)
    if ( .not.flag_noncollinear ) return
    if ( .not.allocated(vxc_mat) ) return
    !call write_border( 1, "op_xc_noncollinear(start)" )
    hpsi(:,1,1) = hpsi(:,1,1) + vxc_mat(:,1,1)*tpsi(:,1,1) + vxc_mat(:,1,2)*tpsi(:,1,2)
    hpsi(:,1,2) = hpsi(:,1,2) + vxc_mat(:,2,1)*tpsi(:,1,1) + vxc_mat(:,2,2)*tpsi(:,1,2)
    !call write_border( 1, "op_xc_noncollinear(end)" )
  END SUBROUTINE op_xc_noncollinear


  SUBROUTINE io_write_noncollinear( rank, flag )

    implicit none
    integer,intent(IN) :: rank
    logical,intent(IN) :: flag
    character(21) :: fln
    character( 2) :: nn
    logical :: flag_exist
    integer :: unit=1, i,j, ml
    complex(8),allocatable :: zw(:)
    type(time) :: tt
!    character( 8) :: yyyymmdd
!    character(10) :: hhmmsspsss
!    character(12) :: label

    if ( OC2 < 1 .or. OC < 1 .or. OC > 15 ) return
    if ( .not.( flag .or. icount==OC2 ) ) return

    call write_border( 0, " io_write_noncollinear(start)" )

    call start_timer( tt )

!    call date_and_time( DATE=yyyymmdd, TIME=hhmmsspsss )
!    label=yyyymmdd(3:8)//hhmmsspsss(1:6)
!    fln="ncol.dat."//label

    if ( rank == 0 ) then

       do i=1,99
          write(nn,'(i2.2)') i
          fln="ncol.dat"//nn
          inquire( FILE=fln, EXIST=flag_exist )
          if ( .not.flag_exist ) exit
       end do

       open(unit,file=fln,form='unformatted')

    end if

    ml=sum(ir_grid)
    allocate( zw(ml) ) ; zw=(0.0d0,0.0d0)

    do j=1,size(den_mat,3)
    do i=1,size(den_mat,2)
       call rsdft_allgatherv( den_mat(:,i,j), zw, ir_grid, id_grid, comm_grid )
       if ( rank == 0 ) write(unit) zw
    end do
    end do
    do j=1,size(vxc_mat,3)
    do i=1,size(vxc_mat,2)
       call rsdft_allgatherv( vxc_mat(:,i,j), zw, ir_grid, id_grid, comm_grid )
       if ( rank == 0 ) write(unit) zw
    end do
    end do

    deallocate( zw )

    if ( rank == 0 ) close(unit)

    call result_timer( "io_write_noncollinear", tt )

    call write_border( 0, " io_write_noncollinear(end)" )

  END SUBROUTINE io_write_noncollinear


  SUBROUTINE io_read_noncollinear( rank, flag_read )

    implicit none
    integer,intent(IN)  :: rank
    logical,intent(OUT) :: flag_read
    logical :: flag_exist
    character(11) :: fln
    character( 2) :: nn
    integer,parameter :: unit=2
    integer :: ml,i,j,ierr
    complex(8),allocatable :: zw1(:),zw3(:,:,:)
    include 'mpif.h'

    flag_read = .false.

    if ( IC <= 0 ) return
    if ( .not.flag_noncollinear ) return

    call write_border( 0, " io_read_noncollinear(start)" )

    if ( rank == 0 ) then
       fln="ncol.dat"
       inquire( FILE=fln, EXIST=flag_exist )
       if ( .not. flag_exist ) then
          do i=99,1,-1
             write(nn,'(i2.2)') i
             fln="ncol.dat"//nn
             inquire( FILE=fln, EXIST=flag_exist )
             if ( flag_exist ) exit
          end do
       end if
    end if
    call MPI_BCAST( flag_exist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr )

    if ( .not. flag_exist ) goto 900

    if ( rank == 0 ) then
       write(*,'(1x,"Read from ",a)') fln
       open(unit,file=fln,status="old",form="unformatted")
    end if

    if ( .not.allocated(den_mat) ) call allocate_array( ir_grid(myrank_g) )

    ml=sum(ir_grid)
    allocate( zw1(ml) ) ; zw1=zero
    allocate( zw3(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) ) ; zw3=zero

    do j=1,size(den_mat,3)
    do i=1,size(den_mat,2)
       if ( rank == 0 ) read(unit) zw1
       call MPI_BCAST( zw1, ml, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr )
       call z_convert_1d_to_3d_grid( lat_old, zw1, zw3 )
       call z_convert_3d_to_1d_grid( lat_new, zw3, den_mat(:,i,j) )
    end do
    end do
    do j=1,size(vxc_mat,3)
    do i=1,size(vxc_mat,2)
       if ( rank == 0 ) read(unit) zw1
       call MPI_BCAST( zw1, ml, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr )
       call z_convert_1d_to_3d_grid( lat_old, zw1, zw3 )
       call z_convert_3d_to_1d_grid( lat_new, zw3, vxc_mat(:,i,j) )
    end do
    end do

    deallocate( zw3 )
    deallocate( zw1 )
    call deallocate_io_ctrl_parameters

    if ( rank == 0 ) close(unit)

    flag_read = .true.

    call init_mixing_dm( size(den_mat,1), 4, comm_grid, den_mat )

900 call write_border( 0, " io_read_noncollinear(end)" )

  END SUBROUTINE io_read_noncollinear


  SUBROUTINE allocate_array( ml )
    implicit none
    integer,intent(IN) :: ml
    allocate( den_mat(ml,2,2) ) ; den_mat=zero
    allocate( vxc_mat(ml,2,2) ) ; vxc_mat=zero
    allocate( old_mat(ml,2,2) ) ; old_mat=zero
  END SUBROUTINE allocate_array


END MODULE noncollinear_module
