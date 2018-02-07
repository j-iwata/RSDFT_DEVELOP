MODULE localpot_module

!  use rgrid_module, only: Ngrid, Igrid
  use xc_module, only: Vxc
  use hartree_variables, only: Vh
  use ps_local_module, only: Vion
 !use localpot2_module

  implicit none

  PRIVATE
  PUBLIC :: Vloc, init_localpot, construct_matrix_localpot
  PUBLIC :: op_localpot
  PUBLIC :: update_localpot

  real(8),allocatable :: Vloc(:,:)
  integer :: ML_0,ML_1,MS_0,MS_1
  integer :: MS
  integer,allocatable :: Ngrid(:),Igrid(:,:)

CONTAINS


  SUBROUTINE init_localpot( ML_0_in,ML_1_in, MS_0_in,MS_1_in )
    implicit none
    integer,intent(IN) :: ML_0_in,ML_1_in,MS_0_in,MS_1_in
    call write_border( 0, " init_localpot(start)" )
    ML_0 = ML_0_in
    ML_1 = ML_1_in
    MS_0 = MS_0_in
    MS_1 = MS_1_in
    allocate( Vloc(ML_0:ML_1,MS_0:MS_1) )
    Vloc=0.0d0
    call write_border( 0, " init_localpot(end)" )
  END SUBROUTINE init_localpot


  SUBROUTINE op_localpot( f, Vf, s_in )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN) :: f(:,:)
    real(8),intent(INOUT) :: Vf(:,:)
#else
    complex(8),intent(IN) :: f(:,:)
    complex(8),intent(INOUT) :: Vf(:,:)
#endif
    integer,optional,intent(IN) :: s_in
    integer :: s,m,n,i,j

    m = size( f, 1 )
    n = size( f, 2 )
    s = 1 ; if ( present(s_in) ) s=s_in

!    if ( flag_localpot2 ) then
!
!       allocate( ft(Ngrid(0)) )
!       ft(:)=zero
!
!!!$OMP PARALLEL
!       do n=1,nn
!!$OMP SINGLE
!       call mpi_allgatherv(f(1,n),mm,TYP,ft,ir_grid,id_grid,TYP,comm_grid,ierr)
!!$OMP END SINGLE
!!$OMP DO
!          do i=1,mm
!          do j=1,MLpot
!             vf(i,n)=vf(i,n)+vloc_nl(j,i+ML_0-1)*ft(Lpot(j,i+ML_0-1))
!          end do
!          end do
!!$OMP END DO
!       end do
!!!$OMP END PARALLEL
!
!       deallocate( ft )
!
!    else

!!$OMP PARALLEL
    do j=1,n
!$OMP DO
       do i=1,m
          Vf(i,j) = Vf(i,j) + Vloc(ML_0-1+i,s)*f(i,j)
       end do
!$OMP END DO
    end do
!!$OMP END PARALLEL

!    end if

  END SUBROUTINE op_localpot


  SUBROUTINE read_localpot(file_name,rank)
    implicit none
    character(*),intent(IN) :: file_name
    integer,intent(IN) :: rank
    integer,parameter :: unit=80
    integer :: ML_tmp(0:3),ierr,i1,i2,i3,i,s
    integer,allocatable :: LL_tmp(:,:)
    real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
    include 'mpif.h'
    if ( rank == 0 ) then
       open(unit,file=file_name,form='unformatted')
       read(unit) ML_tmp(0:3)
       write(*,*) "ML_TMP=",ML_tmp
    end if
    call mpi_bcast(ML_tmp,4,mpi_integer,0,mpi_comm_world,ierr)
    if ( ML_tmp(0) /= Ngrid(0) ) stop "at read_localpot"
    allocate( LL_tmp(3,ML_tmp(0)) )
    if ( rank == 0 ) then
       read(unit) LL_tmp(:,:)
    end if
    call mpi_bcast(LL_tmp,3*ML_tmp(0),mpi_integer,0,mpi_comm_world,ierr)
    allocate( rtmp3(0:ML_tmp(1)-1,0:ML_tmp(2)-1,0:ML_tmp(3)-1) )
    allocate( rtmp(ML_tmp(0)) )
    do s=1,MS
       if ( rank == 0 ) then
          read(unit) rtmp(:) ! rho
          read(unit) rtmp(:) ! Vloc
       end if
       call mpi_bcast(rtmp,ML_tmp(0),mpi_real8,0,mpi_comm_world,ierr)
       do i=1,ML_tmp(0)
          rtmp3( LL_tmp(1,i),LL_tmp(2,i),LL_tmp(3,i) ) = rtmp(i)
       end do
       if ( MS_0 <= s .and. s <= MS_1 ) then
          i=ML_0-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             Vloc(i,s) = rtmp3(i1,i2,i3)
          end do
          end do
          end do
       end if
    end do ! s
    if ( rank == 0 ) then
       close(unit)
    end if
    deallocate( rtmp   )
    deallocate( rtmp3  )
    deallocate( LL_tmp )
  END SUBROUTINE read_localpot


  SUBROUTINE construct_matrix_localpot( s, ML, Hmat )
    implicit none
    integer,intent(IN) :: s, ML
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Hmat(ML,ML)
#else
    complex(8),intent(INOUT) :: Hmat(ML,ML)
#endif
    integer :: i

    do i=1,ML
       Hmat(i,i) = Hmat(i,i) + Vloc(i,s)
    end do

  END SUBROUTINE construct_matrix_localpot


  SUBROUTINE update_localpot
    implicit none
    integer :: s
    do s=MS_0,MS_1
       Vloc(:,s) = Vion(:) + Vh(:) + Vxc(:,s)
    end do
  END SUBROUTINE update_localpot


END MODULE localpot_module
