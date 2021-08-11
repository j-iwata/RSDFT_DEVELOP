module subspace_rotv_sl2_module

  implicit none

  private
  public :: subspace_rotv_sl2

contains

  subroutine subspace_rotv_sl2( psi )
    use sl_variables, only: sl2, Hsub
    use parallel_module, only: get_range_parallel
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    complex(8),intent(inout) :: psi(:,:)
    complex(8),allocatable :: Vsub(:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0), z1=(1.0d0,0.0d0)
    complex(8),allocatable :: psiRot(:,:)
    integer :: nb_0, nb_1, nb, ng, i0, i1, ii, ib

    call write_border( 1, ' z_subspace_rotv_sl2(start)' )

    ng = size( psi, 1 )
    call get_range_parallel( nb_0, nb_1, 'b' )
    nb = nb_1 - nb_0 + 1

    allocate( Vsub(sl2%nband,nb) ); Vsub=z0

    if ( sl2%myrow >= 0 ) then
      i0 = sl2%myrow * sl2%mbsize + 1
      i1 = min( i0 + sl2%mbsize - 1, sl2%nband )
      ii = i1 - i0 + 1
      Vsub(i0:i1,1:nb) = Hsub(1:ii,1:nb) 
    end if
    do ib = 1, nb
      call rsdft_allreduce( Vsub(:,ib), 'g' )
    end do

    allocate( psiRot(ng,nb) ); psiRot=z0

    call ZGEMM( 'N','N',ng,nb,sl2%nband, z1, psi, ng, Vsub, sl2%nband, z0, psiRot, ng )

    psi(:,nb_0:nb_1) = psiRot(:,:)

    deallocate( psiRot )
    deallocate( Vsub )

    call write_border( 1, ' z_subspace_rotv_sl2(end)' )

  end subroutine subspace_rotv_sl2

end module subspace_rotv_sl2_module