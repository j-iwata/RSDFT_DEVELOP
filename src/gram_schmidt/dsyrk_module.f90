module dsyrk_module

  use rgrid_variables, only: dV

  implicit none

  private
  public :: calc_dsyrk3

  integer,public :: ialgo_dsyrk=1
  integer,public :: nblk_dsyrk=1

contains

  subroutine calc_dsyrk3( u, S )
    implicit none
    real(8),intent(in) :: u(:,:)
    real(8),intent(inout) :: S(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: utmp(:,:)
    integer :: m,n,nblk,iblk,jblk
    integer :: i0,i1,j0,j1,ni,nj

    m=size(u,1)
    n=size(u,2)
    nblk = max( n/nblk_dsyrk, 1 )

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk, n )
       nj = j1 - j0 + 1

       do iblk = 1, jblk, nblk

          i0 = iblk
          i1 = min( iblk+nblk, n )
          ni = i1 - i0 + 1

          if ( iblk == jblk ) then
             if ( ialgo_dsyrk == 1 ) then
                call DGEMM( 'T','N',ni,nj,m,dV,u(1,i0),m,u(1,j0),m,zero,S(i0,j0),n )
             else
                call DSYRK( 'U', 'C', ni, m, dV, u(1,i0), m, zero, S(i0,i0), n )
             end if
          else
             call DGEMM( 'T','N',ni,nj,m,dV,u(1,i0),m,u(1,j0),m,zero,S(i0,j0),n )
          end if

       end do !iblk

    end do !jblk

!    allocate( utmp(n,m) ) !; utmp=transpose(u)
!    call tr
!    call DGEMM( 'N', 'N', n, n, m, dV, utmp, n, u, m, zero, S, n )
!    call DGEMM( 'T', 'N', n, n, m, dV, u, m, u, m, zero, S, n )
!    deallocate( utmp )
!  contains
!    subroutine tr
!       utmp=transpose(u)
!    end subroutine tr
  end subroutine calc_dsyrk3

end module dsyrk_module

