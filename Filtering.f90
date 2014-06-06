MODULE Filtering
  implicit none

CONTAINS

!------------------------------------
  SUBROUTINE opFiltering( qc,L,NRc,NRps,rad,rad1,vrad,tarIN )
    implicit none

    real(8),intent(IN) :: qc
    integer,intent(IN) :: L,NRc,NRps
    real(8),allocatable,intent(IN) :: rad(:),rad1(:),vrad(:)
    real(8),allocatable,intent(INOUT) :: tarIN(:)
    integer :: i,j
    real(8) :: r,r1
    real(8) :: sb0x,sb1x,sb0y,sb1y,sum0
    real(8),allocatable :: tmp(:)
    real(8),parameter :: const=0.d0 ! HERE

    allocate( tmp(NRc) ) ; tmp(:)=0.d0
    do i=1,NRps
      r=rad1(i)

      select case(L)
      case(0)

      end select
    end do


  END SUBROUTINE opFiltering

END MODULE Filtering
