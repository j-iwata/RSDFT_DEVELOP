MODULE phase_control_module

  use parallel_module, only: comm_grid

  implicit none

  PRIVATE
  PUBLIC :: phase_control

CONTAINS

  SUBROUTINE phase_control( m0, m1, unk )
    implicit none
    integer,intent(IN) :: m0,m1
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: unk(:,:,:,:)
#else
    complex(8),intent(INOUT) :: unk(:,:,:,:)
#endif
    integer :: s,k,n,ierr
    real(8) :: a,b,r,co,si
    complex(8) :: z,p
    include 'mpif.h'

    do s=1,size(unk,4)
    do k=1,size(unk,3)
    do n=1,size(unk,2)
       if ( m0 == 1 ) then
          z=unk(1,n,k,s)
          r=abs(z)
          co=real(z)/r
          si=aimag(z)/r
          p=dcmplx( co, -si )
!          z=unk(1,n,k,s)
!          r=abs(z)
!          co=real(z)/r
!          si=aimag(z)/r
!          write(*,'(1x,3i4,2f20.15)') n,k,s,acos(co),asin(si)
       end if
       call MPI_BCAST( p,1,MPI_COMPLEX16,0,comm_grid,ierr )
       unk(:,n,k,s)=unk(:,n,k,s)*p
    end do
    end do
    end do

  END SUBROUTINE phase_control

END MODULE phase_control_module
