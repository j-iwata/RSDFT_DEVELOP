module gram_schmidt_m_bp_module

  use vector_tools_module, only: vinfo, inner_product_vector, normalize_vector
  use rsdft_mpi_module

  implicit none

  private
  public :: gram_schmidt_m_bp

contains

  subroutine gram_schmidt_m_bp( u, v )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(inout) :: u(:,:)
    real(8),allocatable :: w(:,:)
    real(8) :: uu,vv
#else
    complex(8),intent(inout) :: u(:,:)
    complex(8),allocatable :: w(:,:)
    complex(8) :: uu,vv
#endif
    type(vinfo),intent(in) :: v(2)
    integer :: m, nb, ib, jb, n1, n2
    integer :: np_band, comm_band
    integer,allocatable :: ir(:),id(:)

    call write_border( 1, "gram_schmidt_m_bp(start)" )

    m = size( u, 1 )
    nb = sum( v(2)%pinfo%ir )

    allocate( w(m,nb) ); w=0.0d0

    np_band   = v(2)%pinfo%np
    comm_band = v(2)%pinfo%comm
    allocate( ir(0:np_band-1) ); ir=m*v(2)%pinfo%ir
    allocate( id(0:np_band-1) ); id=m*v(2)%pinfo%id
    call rsdft_allgatherv( u, w, ir, id, comm_band )
    deallocate( id )
    deallocate( ir )

    do jb=1,nb
       do ib=1,jb-1
          call inner_product_vector( w(:,ib), w(:,jb), uu, v(1) )
          w(:,jb) = w(:,jb) - w(:,ib)*uu
       end do ! ib
       call normalize_vector( w(:,jb), v(1) )
    end do ! jb

    n1 = v(2)%pinfo%id( v(2)%pinfo%me ) + 1
    n2 = n1 + v(2)%pinfo%ir( v(2)%pinfo%me ) - 1

    u(:,:) = w(:,n1:n2)

    deallocate( w )

    call write_border( 1, "gram_schmidt_m_bp(end)" )

  end subroutine gram_schmidt_m_bp

end module gram_schmidt_m_bp_module
