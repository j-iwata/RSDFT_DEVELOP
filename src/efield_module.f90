module efield_module

  use aa_module, only: aa
  use rgrid_variables, only: Igrid, Ngrid, Hgrid
  use io_tools_module
  use parallel_module, only: comm_grid

  implicit none

  private
  public :: sawtooth_efield
  public :: force_sawtooth_efield

  real(8) :: Evec(3)
  integer :: icontrol = 0
  real(8) :: zero_pot_pos(3)
  real(8),allocatable :: Vpot(:,:,:)

contains

  subroutine init_efield
    implicit none
    real(8) :: tmp(6)
    integer :: a1,a2,a3,b1,b2,b3,ierr,i1,i2,i3
    real(8) :: ee(3,3),x,y,z,vx,vy,vz,x0,y0,z0
    real(8) :: aax,aay,aaz
    Evec(:)=0.0d0
    tmp(:)=0.0d0
    call IOTools_readReal8Keyword( "EFIELD", tmp )
    Evec(1:3)=tmp(1:3)
    icontrol = -1
    if ( any( Evec /= 0.0d0 ) ) icontrol = 1
    zero_pot_pos = tmp(4:6)
    if ( icontrol == -1 ) return
    !
    a1=Igrid(1,1); b1=Igrid(2,1)
    a2=Igrid(1,2); b2=Igrid(2,2)
    a3=Igrid(1,3); b3=Igrid(2,3)
    allocate( Vpot(a1:b1,a2:b2,a3:b3) ); Vpot=0.0d0
    ee(:,1)=aa(:,1)/Ngrid(1)
    ee(:,2)=aa(:,2)/Ngrid(2)
    ee(:,3)=aa(:,3)/Ngrid(3)
    aax=sqrt(sum(aa(:,1)**2))
    aay=sqrt(sum(aa(:,2)**2))
    aaz=sqrt(sum(aa(:,3)**2))
    !
    x0 = zero_pot_pos(1)
    y0 = zero_pot_pos(2)
    z0 = zero_pot_pos(3)
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
      x = ee(1,1)*i1+ee(1,2)*i2+ee(1,3)*i3 - x0
      y = ee(2,1)*i1+ee(2,2)*i2+ee(2,3)*i3 - y0
      z = ee(3,1)*i1+ee(3,2)*i2+ee(3,3)*i3 - z0
      vx = x*Evec(1)
      if ( x > aax*0.5d0 ) then
        vx = x*Evec(1) - Evec(1)*aax
      else if ( x < -aax*0.5d0 ) then
        vx = x*Evec(1) + Evec(1)*aax
      end if
      vy = y*Evec(2)
      if ( y > aay*0.5d0 ) then
        vy = y*Evec(2) - Evec(2)*aay
      else if ( y < -aay*0.5d0 ) then
        vy = y*Evec(2) + Evec(2)*aay
      end if
      vz = z*Evec(3)
      if ( z > aaz*0.5d0 ) then
        vz = z*Evec(3) - Evec(3)*aaz
      else if ( z < -aaz*0.5d0 ) then
        vz = z*Evec(3) + Evec(3)*aaz
      end if
      Vpot(i1,i2,i3) = vx + vy + vz
    end do
    end do
    end do
    !rewind 10
    !do i3=a3,b3
    !  write(10,*) i3*Hgrid(3),Vpot(a1,a2,i3)
    !end do
    !call stop_program('init_efield')
  end subroutine init_efield


  subroutine sawtooth_efield( v )
    implicit none
    real(8),intent(inout) :: v(:)
    real(8) :: ee(3,3),x,y,z
    integer :: i1,i2,i3,i
    if ( icontrol ==  0 ) call init_efield
    if ( icontrol == -1 ) return
    !ee(:,1)=aa(:,1)/Ngrid(1)
    !ee(:,2)=aa(:,2)/Ngrid(2)
    !ee(:,3)=aa(:,3)/Ngrid(3)
    i=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
      i=i+1
      !x=ee(1,1)*i1+ee(1,2)*i2+ee(1,3)*i3 - zero_pot_pos(1)
      !y=ee(2,1)*i1+ee(2,2)*i2+ee(2,3)*i3 - zero_pot_pos(2)
      !z=ee(3,1)*i1+ee(3,2)*i2+ee(3,3)*i3 - zero_pot_pos(3)
      !v(i) = v(i) + x*Evec(1) + y*Evec(2) + z*Evec(3)
      v(i) = v(i) + Vpot(i1,i2,i3)
    end do
    end do
    end do
  end subroutine sawtooth_efield


  subroutine force_sawtooth_efield( force, ki, zion )
    implicit none
    real(8),intent(inout) :: force(:,:)
    integer,intent(in) :: ki(:)
    real(8),intent(in) :: zion(:)
    integer :: a
    if ( icontrol == -1 ) return
    do a = 1, size(ki)
      force(:,a) = force(:,a) + zion(ki(a))*Evec(:)
    end do
  end subroutine force_sawtooth_efield

end module efield_module
