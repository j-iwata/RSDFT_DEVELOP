module subspace_rotv_sl_bp0_module

  use vector_tools_module, only: vinfo
  use scalapack_module, only: sl
  use subspace_diag_variables, only: Vsub
  use rsdft_mpi_module
  use sl_tools_module, only: gather_matrix

  implicit none

  private
  public :: subspace_rotv_sl_bp0

contains

  subroutine subspace_rotv_sl_bp0( u, v )
    implicit none
    real(8),intent(inout) :: u(:,:)
    type(vinfo),intent(in) :: v(2)
    integer :: ml,nb,nb0,m1,m2,np_band,comm_band,n1,n2,n0
    integer,allocatable :: ir(:),id(:)
    real(8) :: zz,dv
    real(8),allocatable :: psiu(:,:), R(:,:)

    call write_border( 1, " subspace_mate_sl_bp0(start)" )

    ml = size( u, 1 )
    nb = sum( v(2)%pinfo%ir )
    n0 = v(2)%pinfo%ir( v(2)%pinfo%me )
    n1 = v(2)%pinfo%id( v(2)%pinfo%me ) + 1
    n2 = n1 + n0 - 1
    m1 = v(1)%pinfo%id( v(1)%pinfo%me ) + 1
    m2 = m1 + v(1)%pinfo%ir( v(1)%pinfo%me ) - 1
    dv = v(1)%factor
    zz = 0.5d0*dv
    np_band   = v(2)%pinfo%np
    comm_band = v(2)%pinfo%comm

    allocate( psiu(ml,nb) ); psiu=0.0d0
    allocate( R(nb,nb) ); R=0.0d0

    allocate( ir(0:np_band-1) ); ir=ml*v(2)%pinfo%ir
    allocate( id(0:np_band-1) ); id=ml*v(2)%pinfo%id
    call rsdft_allgatherv( u, psiu, ir, id, comm_band )
    !call MPI_Allgatherv( u,size(u),MPI_REAL8,psiu,ir,id,MPI_REAL8,comm_band,i)
    deallocate( id )
    deallocate( ir )

    call gather_matrix( sl, Vsub, R, sl%icontxt_a )

    call DGEMM('N','N',ml,n0,nb,1.0d0,psiu,ml,R(1,n1),nb,0.0d0,u,ml)

    deallocate( R )
    deallocate( psiu )

    call write_border( 1, "subspace_rotv_sl_bp0(end)" )

  end subroutine subspace_rotv_sl_bp0

end module subspace_rotv_sl_bp0_module
