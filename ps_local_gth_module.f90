MODULE ps_local_gth_module

  use pseudopot_module
  use ps_gth_module
  use aa_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_local_gth

CONTAINS


  SUBROUTINE init_ps_local_gth(NMGL,MKI,GG,vqlg)
    implicit none
    integer,intent(IN)  :: NMGL,MKI
    real(8),intent(IN)  :: GG(NMGL)
    real(8),intent(OUT) :: vqlg(NMGL,MKI)
    real(8) :: pi,rloc,const,C1,C2,C3,C4,G,v
    integer :: ig,ielm

    pi = acos(-1.d0)
    const = 1.d0/abs(Va)

    vqlg(:,:) = 0.d0

    do ielm=1,MKI

       rloc = Rcloc(ielm)
       C1   = parloc(1,ielm)
       C2   = parloc(2,ielm)
       C3   = parloc(3,ielm)
       C4   = parloc(4,ielm)

       vqlg(1,ielm) = const*( 2.d0*pi*Zps(ielm)*rloc**2 &
            + sqrt((2.d0*pi)**3)*rloc**3*( C1+C2*3.d0+C3*15.d0+C4*105.d0 ) )

       do ig=2,NMGL

          G = sqrt( GG(ig) )

!          if ( G == 0.d0 ) then
!             v = 2.d0*pi*Zps(ielm)*rloc**2
!          else
             v = -4.d0*pi*Zps(ielm)*exp(-0.5d0*(G*rloc)**2)/G**2
!          end if

          v = v + sqrt((2.d0*pi)**3)*rloc**3*exp(-0.5d0*(rloc*G)**2) &
            *( C1 &
              +C2*(3.d0-(rloc*G)**2) &
              +C3*(15.d0-10.d0*(rloc*G)**2+(rloc*G)**4) &
              +C4*(105.d0-105.d0*(rloc*G)**2+21.d0*(rloc*G)**4-(rloc*G)**6) )

          vqlg(ig,ielm) = const*v

       end do ! ig

    end do ! ielm

    return

  END SUBROUTINE init_ps_local_gth


END MODULE ps_local_gth_module
