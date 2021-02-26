module ps_local_fft0_module

  use ggrid_module, only: LLG,construct_Ggrid,destruct_Ggrid,allgatherV_Ggrid,MGL
  use watch_module
  use fft_module, only: backward_fft, forward_fft,z3_to_d1_fft,init_fft,finalize_fft

  implicit none

  private
  public :: construct_ps_local_fft0

contains

  subroutine construct_ps_local_fft0( vqlg, SGK, Vion )
    implicit none
    real(8),intent(in) :: vqlg(:,:)
    complex(8),intent(in) :: SGK(:,:)
    real(8),intent(out) :: Vion(:)
    integer :: i,j,ik,MG,MG_0,MG_1,nelem
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:),vg(:)
    real(8) :: ctt(0:3),ett(0:3)
    logical :: disp_sw

    call write_border( 0, " construct_ps_local_fft0(start)" )

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    call construct_Ggrid(2)

    MG = size( LLG, 2 )
    MG_0 = lbound( SGK, 1 )
    MG_1 = ubound( SGK, 1 )
    nelem = size( SGK, 2 )

    allocate( vg(MG) )

    do i=MG_0,MG_1
      j=MGL(i)
      vg(i)=vqlg(j,1)*SGK(i,1)
    end do
    do ik=2,nelem
      do i=MG_0,MG_1
        j=MGL(i)
        vg(i)=vg(i)+vqlg(j,ik)*SGK(i,ik)
      end do
    end do
    call allgatherv_Ggrid(vg)

    call init_fft( allocarray=zwork0 )

    do i=1,MG
      zwork0(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    call watch(ctt(1),ett(1))

    call backward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

    call z3_to_d1_fft( zwork0, Vion )

    call finalize_fft

    if ( allocated(zwork1) ) deallocate( zwork1 )
    if ( allocated(zwork0) ) deallocate( zwork0 )

    call watch(ctt(3),ett(3))

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!      write(*,*) "time(construct_ps_local_fft0_1)",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(construct_ps_local_fft0_2)",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(construct_ps_local_fft0_3)",ctt(3)-ctt(2),ett(3)-ett(2)
!    end if

    call write_border( 0, " construct_ps_local_fft0(end)" )

  end subroutine construct_ps_local_fft0

end module ps_local_fft0_module
