MODULE force_correction_module

  use ps_initrho_module
  use rgrid_variables, only: Hgrid, Igrid, dV, Ngrid
  use hartree_module, only: calc_hartree
  use xc_module, only: calc_xc
  use density_module, only: rho, rho_in
  use rsdft_mpi_module
  use parallel_module, only: comm_grid

  implicit none

  PRIVATE
  PUBLIC :: force_correction

CONTAINS


  SUBROUTINE force_correction( force )

    implicit none
    real(8),intent(INOUT) :: force(:,:)
    integer :: natm,iatm,i,Ngr,Nsp
    real(8) :: d1(3,3),d2(3,3),fc(3),Edummy
    real(8),allocatable :: vht1(:),vht2(:),vxc1(:,:),vxc2(:,:)
    real(8),allocatable :: rho_ini(:,:),rho_tmp(:,:),rhot(:,:)
    logical :: disp

    call check_disp_switch( disp, 0 )

    d1=0.0d0
    d2=0.0d0
    do i=1,3
       d1(i,i)= 1.0d0/Ngrid(i)
       d2(i,i)=-1.0d0/Ngrid(i)
    end do

    natm = size(force,2)

    Ngr=size(rho,1)
    Nsp=size(rho,2)

    allocate( rho_ini(Ngr,Nsp) ) ; rho_ini=0.0d0
    allocate( rho_tmp(Ngr,Nsp) ) ; rho_tmp=0.0d0
    allocate( rhot(Ngr,Nsp) ) ; rhot=0.0d0
    call get_ps_initrho( rho_ini )

    allocate( vht1(Ngr) ) ; vht1=0.0d0
    allocate( vht2(Ngr) ) ; vht2=0.0d0
    allocate( vxc1(Ngr,1) ) ; vxc1=0.0d0
    allocate( vxc2(Ngr,1) ) ; vxc2=0.0d0

    do iatm=1,natm

       if ( disp ) write(*,'(1x,"iatm/natm=",i3,"/",i3)') iatm,natm

       do i=1,3
          rho_tmp=rho_ini
          call derivative_ps_initrho_r( rho_tmp, iatm, d1(:,i) )
          call calc_hartree( Igrid(1,0),Igrid(2,0),Nsp, rho_tmp, vht1 )
          call calc_xc( rho_tmp, vxc1, Edummy )
          rho_tmp=rho_ini
          call derivative_ps_initrho_r( rho_tmp, iatm, d2(:,i) )
          call calc_hartree( Igrid(1,0),Igrid(2,0),Nsp, rho_tmp, vht2 )
          call calc_xc( rho_tmp, vxc2, Edummy )
          vht1(:) = ( vht1(:) - vht2(:) )/( 2.0d0*Hgrid(i) )
          vxc1(:,1) = ( vxc1(:,1) - vxc2(:,1) )/( 2.0d0*Hgrid(i) )
          fc(i) = sum( vht1(:)*(rho(:,1)-rho_in(:,1)) )*dV &
                + sum( vxc1(:,1)*(rho(:,1)-rho_in(:,1)) )*dV
       end do

       call rsdft_allreduce_sum( fc, comm_grid )

       force(1:3,iatm) = force(1:3,iatm) - fc(1:3)

    end do ! iatm

    deallocate( vxc2 )
    deallocate( vxc1 )
    deallocate( vht2 )
    deallocate( vht1 )
    deallocate( rho_tmp )
    deallocate( rho_ini )

  END SUBROUTINE force_correction

END MODULE force_correction_module
