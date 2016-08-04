MODULE hpc_module

  use io_tools_module
  use aa_module, only: aa

  implicit none

  PRIVATE
  PUBLIC :: force_projection_hpc

  real(8), allocatable :: vnormal(:,:)
  integer :: hpc_switch=0
  logical :: flag_init=.true.

CONTAINS


  SUBROUTINE read_hpc( Natom )
    implicit none
    integer,intent(IN) :: Natom
    integer,parameter :: unit_vnormal=973
    integer :: rank,i
    include 'mpif.h'
    call IOTools_readIntegerKeyword( "HPC", hpc_switch )
    if ( hpc_switch == 1 ) then
       allocate( vnormal(1:3,1:Natom) ) ; vnormal=0.0d0
       call MPI_COMM_RANK( MPI_COMM_WORLD, rank, i )
       if ( rank == 0 ) then
          write(*,*) repeat("-",10),' hpc normal vector: vnormal ',repeat("-",10)
          do i=1,Natom
             read(unit_vnormal,*) vnormal(1:3,i)
          end do
       end if
       call MPI_BCAST( vnormal, size(vnormal), MPI_REAL8, 0, MPI_COMM_WORLD, i )
    end if
    flag_init = .false.
  END SUBROUTINE read_hpc


  SUBROUTINE force_projection_hpc( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(INOUT) :: force(1:3,1:MI)
    real(8) :: vnormal_xyz(1:3,1:MI)        
    integer :: i,it
    real(8) :: sq_vnormal, projection       

    if ( flag_init ) call read_hpc( MI )
    if ( hpc_switch /= 1 ) return

    call write_border( 1, " force_projection_hpc(start)" )

    vnormal_xyz(:,:) = matmul( aa(:,:), vnormal(:,:) )

    sq_vnormal = sum( vnormal_xyz**2 )
    projection = sum( vnormal_xyz*force )
    projection = projection/sq_vnormal

    do it=1,MI
       do i=1,3
          force(i,it) = force(i,it) - projection*vnormal_xyz(i,it)
       end do
    end do

    call write_border( 1, " force_projection_hpc(end)" )

  END SUBROUTINE force_projection_hpc


END MODULE hpc_module
