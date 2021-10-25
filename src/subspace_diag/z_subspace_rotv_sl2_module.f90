module z_subspace_rotv_sl2_module

  implicit none

  private
  public :: z_subspace_rotv_sl2

contains

  subroutine z_subspace_rotv_sl2( psi )
    use sl_variables, only: sl2, Hsub
    use parallel_module, only: get_range_parallel
    use rsdft_allreduce_module, only: rsdft_allreduce
    use watch_module, only: watchb
    implicit none
    complex(8),intent(inout) :: psi(:,:)
    complex(8),allocatable :: Vsub(:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0), z1=(1.0d0,0.0d0)
    complex(8),allocatable :: psiRot(:,:)
    integer :: nb_0, nb_1, nb, ng, i0, i1, ii, ib
    real(8) :: ttmp(2),tt(2,3),ttmin(2,3),ttmax(2,3)
    logical :: disp_on

    call write_border( 1, ' z_subspace_rotv_sl2(start)' )

    call watchb( ttmp, barrier='on' )

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

    call watchb( ttmp, tt(:,1) )

    call ZGEMM( 'N','N',ng,nb,sl2%nband, z1, psi, ng, Vsub, sl2%nband, z0, psiRot, ng )

    call watchb( ttmp, tt(:,2) )

    psi(:,nb_0:nb_1) = psiRot(:,:)

    deallocate( psiRot )
    deallocate( Vsub )

    call watchb( ttmp, tt(:,3) )

    ! ttmin(:,1:3)=tt(:,1:3); call rsdft_ar( ttmin, op_in='min' )
    ! ttmax(:,1:3)=tt(:,1:3); call rsdft_ar( ttmax, op_in='max' )
    ! call check_disp_switch( disp_on, 0 )
    ! if ( disp_on ) then
    !   write(*,'(1x,"time_rot(1)=",2f12.5,2x,2f12.5)') ttmin(:,1),ttmax(:,1)
    !   write(*,'(1x,"time_rot(2)=",2f12.5,2x,2f12.5)') ttmin(:,2),ttmax(:,2)
    !   write(*,'(1x,"time_rot(3)=",2f12.5,2x,2f12.5)') ttmin(:,3),ttmax(:,3)
    ! end if

    call write_border( 1, ' z_subspace_rotv_sl2(end)' )

  end subroutine z_subspace_rotv_sl2

end module z_subspace_rotv_sl2_module
