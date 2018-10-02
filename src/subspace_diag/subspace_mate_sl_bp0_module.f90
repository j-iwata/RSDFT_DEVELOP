module subspace_mate_sl_bp0_module

  use vector_tools_module, only: vinfo
  use rsdft_mpi_module
  use subspace_diag_variables, only: Hsub
  use sl_tools_module, only: distribute_matrix
  use scalapack_module, only: UPLO, sl
  use hamiltonian_module

  implicit none

  private
  public :: subspace_mate_sl_bp0

contains

  subroutine subspace_mate_sl_bp0( k, s, u, v )
    implicit none
    integer,intent(in) :: k,s
    real(8),intent(inout) :: u(:,:)
    type(vinfo),intent(in) :: v(2)
    integer :: ib,ml,nb,m1,m2,np_band,comm_band,comm_grid,i
    integer,allocatable :: ir(:),id(:)
    real(8),allocatable :: psiu(:,:), psit(:,:), H(:,:)
    real(8) :: dv,zz
    include 'mpif.h'

    call write_border( 1, " subspace_mate_sl_bp0(start)" )

    ml = size( u, 1 )
    nb = sum( v(2)%pinfo%ir )
    m1 = v(1)%pinfo%id( v(1)%pinfo%me ) + 1
    m2 = m1 + v(1)%pinfo%ir( v(1)%pinfo%me ) - 1
    dv = v(1)%factor
    zz = 0.5d0*dv
    np_band   = v(2)%pinfo%np
    comm_band = v(2)%pinfo%comm
    comm_grid = v(1)%pinfo%comm

    allocate( psiu(ml,nb) ); psiu=0.0d0
    allocate( psit(ml,nb) ); psit=0.0d0
    allocate( H(nb,nb) ); H=0.0d0

    allocate( ir(0:np_band-1) ); ir=ml*v(2)%pinfo%ir
    allocate( id(0:np_band-1) ); id=ml*v(2)%pinfo%id
    call rsdft_allgatherv( u, psiu, ir, id, comm_band )
    !call MPI_Allgatherv( u,size(u),MPI_REAL8,psiu,ir,id,MPI_REAL8,comm_band,i)
    deallocate( id )
    deallocate( ir )

    do ib=1,nb
       call hamiltonian(k,s,psiu(1,ib),psit(1,ib),m1,m2,ib,ib)
    end do

    UPLO="L"
    call DSYR2K(UPLO,'T',nb,ml,zz,psiu,ml,psit,ml,0.0d0,H,nb)
    !call DGEMM('T','N',nb,nb,ml,dv,psiu,ml,psit,ml,0.0d0,H,nb)
    call rsdft_allreduce_sum( H, comm_grid )

    call distribute_matrix( sl, H, Hsub )

    deallocate( H )
    deallocate( psit )
    deallocate( psiu )

    call write_border( 1, " subspace_mate_sl_bp0(end)" )

  end subroutine subspace_mate_sl_bp0

end module subspace_mate_sl_bp0_module
