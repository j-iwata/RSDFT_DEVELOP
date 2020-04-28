module ps_nloc3_gth_module

  implicit none
  private
  public :: init_ps3_gth

contains

  subroutine init_ps3_gth( norb, no, lo, ro, hnl, gk, viodgk )
    implicit none
    integer,intent(in) :: norb(:)
    integer,intent(in) :: no(:,:)
    integer,intent(in) :: lo(:,:)
    real(8),intent(in) :: ro(:,:)
    real(8),intent(in) :: hnl(:,0:,:)
    real(8),intent(in) :: gk(:)
    real(8),intent(inout) :: viodgk(:,:,:)
    integer :: ielm,Nelm,iorb,l,i,mpgk,n
    real(8) :: q,pi,c

    pi = acos(-1.0d0)

    Nelm = size(norb)
    mpgk = size(viodgk,1)

    do ielm = 1, Nelm
      do iorb = 1, norb(ielm)
      
        n = no(iorb,ielm)
        l = lo(iorb,ielm)
        c = sqrt(hnl(n,l,ielm))

        select case( l )
        case( 0 )
        
          select case( n )
          case( 1 )
            do i = 1, mpgk
              q = gk(i)
              viodgk(i,iorb,ielm) = c*4.0d0*ro(iorb,ielm) &
              *sqrt(2.0d0*ro(iorb,ielm))*pi**1.25d0 &
              *exp(-0.5d0*(q*ro(iorb,ielm))**2)
            end do
          case( 2 )
            do i = 1, mpgk
              q = gk(i)
              viodgk(i,iorb,ielm) = c*8.0d0*ro(iorb,ielm) &
              *sqrt(2.0d0*ro(iorb,ielm)/15.0d0)*pi**1.25d0 &
              *exp(-0.5d0*(q*ro(iorb,ielm))**2) &
              *(3.0d0-(q*ro(iorb,ielm))**2)
            end do
          end select
        
        case( 1 )
        
          do i = 1, mpgk
            q = gk(i)
            viodgk(i,iorb,ielm) = c*8.0d0*ro(iorb,ielm)**2 &
            *sqrt(ro(iorb,ielm)/3.0d0)*pi**1.25d0 &
            *exp(-0.5d0*(q*ro(iorb,ielm))**2) &
            *q
          end do
        
        end select

      end do !iorb
    end do !ielm

  end subroutine init_ps3_gth

end module ps_nloc3_gth_module
