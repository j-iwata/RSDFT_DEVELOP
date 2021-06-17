module hartree_sol_fft0_module

  use hartree_variables, only: E_hartree, Vh
  use rgrid_module, only: Igrid,dV
  use ggrid_module, only: GG,LLG,MGL,construct_ggrid,destruct_ggrid
  use parallel_module, only: comm_grid
  use fft_module, only: d1_to_z3_fft, z3_to_d1_fft, forward_fft, backward_fft, init_fft, finalize_fft
  use rsdft_mpi_module, only: rsdft_allreduce
  use watch_module

  implicit none

  private
  public :: calc_hartree_sol_fft0

contains

  subroutine calc_hartree_sol_fft0( rho )
    implicit none
    real(8),intent(in) :: rho(:,:)
    integer :: i,i1,i2,i3,ng,ns
    real(8) :: Eh0,pi4,g2,ctt(0:5),ett(0:5)
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    logical :: disp_sw

    call write_border( 1, " calc_hartree_sol_fft0(start)" )

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    ng = size( rho, 1 )
    ns = size( rho, 2 )

    allocate( work(ng) ); work=0.0d0

    do i=1,ns
      work(:) = work(:) + rho(:,i)
    end do

    call init_fft
    call d1_to_z3_fft( work, zwork0 )

    call watch(ctt(1),ett(1))

    call forward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

! ---

    call construct_Ggrid(2)

    pi4 = 4.0d0*acos(-1.0d0)

    zwork1(:,:,:)=z0
    do i=1,size(MGL)
      g2=GG(MGL(i))
      if ( g2 == 0.0d0 ) cycle
      i1=LLG(1,i)
      i2=LLG(2,i)
      i3=LLG(3,i)
      zwork1(i1,i2,i3)=zwork0(i1,i2,i3)*pi4/g2
    end do

    call destruct_Ggrid

! ---

    call watch(ctt(3),ett(3))

    call backward_fft( zwork1, zwork0 )

    call watch(ctt(4),ett(4))

    call finalize_fft

    if ( allocated(zwork0) ) deallocate( zwork0 )

    call z3_to_d1_fft( zwork1, Vh )

    Eh0=0.0d0
    i=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Eh0 = Eh0 + dble( zwork1(i1,i2,i3) )*work(i)
    end do
    end do
    end do
    E_hartree = 0.5d0*Eh0*dV
    call rsdft_allreduce( E_hartree, comm_grid )

    if ( allocated(zwork1) ) deallocate( zwork1 )
    deallocate( work )

    call watch(ctt(5),ett(5))

!    call check_disp_switch( disp_sw, 0 )
!    if ( disp_sw ) then
!      write(*,*) "time(calc_hatree_sol_fft0_1)=",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(calc_hatree_sol_fft0_2)=",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(calc_hatree_sol_fft0_3)=",ctt(3)-ctt(2),ett(3)-ett(2)
!      write(*,*) "time(calc_hatree_sol_fft0_4)=",ctt(4)-ctt(3),ett(4)-ett(3)
!      write(*,*) "time(calc_hatree_sol_fft0_5)=",ctt(5)-ctt(4),ett(5)-ett(4)
!    end if

    call write_border( 1, " calc_hartree_sol_fft0(end)" )

  end subroutine calc_hartree_sol_fft0

end module hartree_sol_fft0_module
