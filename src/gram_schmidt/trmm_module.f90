module trmm_module

  implicit none

  private
  public :: calc_dtrmm3, calc_ztrmm3

  integer :: ialgo_dtrmm=1
  integer :: nblk_dtrmm=1

contains

  subroutine calc_dtrmm( u, S )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: S(:,:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n
    real(8),allocatable :: utmp(:,:)
    m=size(u,1)
    n=size(u,2)
    call dtrmm( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )
!    call alloc_utmp
!    utmp=u
!    call DGEMM( 'N', 'N', m, n, n, one, utmp, m, S, n, zero, u, m )
!    call dealloc_utmp
  contains
    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      deallocate( utmp )
    end subroutine dealloc_utmp
  end subroutine calc_dtrmm


  subroutine calc_dtrmm3( u, S, b0 )
    implicit none
    real(8),intent(inout) :: u(:,:)
    integer,intent(in) :: b0
    real(8),intent(in) :: S(:,b0:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n,j0,j1,i0,i1,b1,iblk,jblk
    integer :: nblk,nj,ni
    real(8),allocatable :: utmp(:,:)

    m=size(u,1)
    n=size(u,2)

    b1 = b0 + size(S,2) - 1

    nblk = max( n/nblk_dtrmm, 1 )

    call alloc_utmp

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk-1, n )

       if ( b1 < j0 .or. j1 < b0 ) cycle
       j0 = max( j0, b0 )
       j1 = min( j1, b1 )

       nj = j1 - j0 + 1

       do iblk = 1, jblk, nblk

          i0 = iblk
          i1 = min( iblk+nblk-1, n )
          ni = i1 - i0 + 1

          if ( iblk == jblk ) then
             if ( ialgo_dtrmm == 1 ) then
                call DGEMM( 'N','N',m,ni,ni,one,u(1,i0),m,S(i0,i0),n,one,utmp(1,i0),m )
             else
                call DTRMM( 'R', 'U', 'N', 'N', m, ni, one, S(i0,i0), n, u(1,i0), m )
             end if
          else
             call DGEMM( 'N','N',m,nj,ni,one,u(1,i0),m,S(i0,j0),n,one,utmp(1,j0),m )
          end if

       end do !iblk

    end do !jblk

    !u=utmp
    call dealloc_utmp

  contains
    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,b0:b1) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      u(:,b0:b1)=utmp
      deallocate( utmp )
    end subroutine dealloc_utmp
  end subroutine calc_dtrmm3

  subroutine calc_ztrmm3( u, S, nblk_in )
    implicit none
    complex(8),intent(inout) :: u(:,:)
    complex(8),intent(in) :: S(:,:)
    integer,optional,intent(in) :: nblk_in
    complex(8),parameter :: zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0)
    complex(8),allocatable :: utmp(:,:)
    integer :: m,n,j0,j1,i0,i1,iblk,jblk
    integer :: nblk,nj,ni

    m=size(u,1)
    n=size(u,2)

    nblk = n
    if ( present(nblk_in) ) nblk = max( n/nblk_in, 1 )

    call alloc_utmp

    do jblk = 1, n, nblk

      j0 = jblk
      j1 = min( jblk+nblk-1, n )
      nj = j1 - j0 + 1

      do iblk = 1, jblk, nblk

        i0 = iblk
        i1 = min( iblk+nblk-1, n )
        ni = i1 - i0 + 1

        if ( iblk == jblk ) then
          call ZGEMM( 'N','N',m,ni,ni,one,u(1,i0),m,S(i0,i0),n,one,utmp(1,i0),m )
         !call ZTRMM( 'R', 'U', 'N', 'N', m, ni, one, S(i0,i0), n, u(1,i0), m )
        else
          call ZGEMM( 'N','N',m,nj,ni,one,u(1,i0),m,S(i0,j0),n,one,utmp(1,j0),m )
        end if

      end do !iblk

    end do !jblk

    !u=utmp
    call dealloc_utmp

  contains
    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      u=utmp
      deallocate( utmp )
    end subroutine dealloc_utmp
  end subroutine calc_ztrmm3


  subroutine calc_dtrmm4( u, S )
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: S(:,:)
    real(8),parameter :: zero=0.0d0,one=1.0d0
    integer :: m,n,j0,j1,i0,i1,iblk,jblk
    integer :: nblk,nj,ni
    real(8),allocatable :: utmp(:,:)
    real(8),allocatable :: stmp(:,:)

    m=size(u,1)
    n=size(u,2)

    nblk = max( n/2, 1 )

    call alloc_utmp

    do iblk = 1, n, nblk
       i0 = iblk
       i1 = min( iblk+nblk-1, n )
       ni = i1 - i0 + 1
       utmp(:,i0:i1) = u(:,i0:i1)
       call DTRMM( 'R', 'U', 'N', 'N', m, ni, one, S(i0,i0), n, utmp(1,i0), m )
    end do !iblk

 ! ---

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk-1, n )
       nj = j1 - j0 + 1

       do iblk = 1, jblk-1, nblk

          i0 = iblk
          i1 = min( iblk+nblk-1, n )
          ni = i1 - i0 + 1

          call DGEMM( 'N','N',m,nj,ni,one,u(1,i0),m,S(i0,j0),n,one,utmp(1,j0),m )

       end do !iblk

    end do !jblk

    !u=utmp
    call dealloc_utmp

  contains

    subroutine alloc_utmp
      implicit none
      allocate( utmp(m,n) ); utmp=zero
    end subroutine alloc_utmp
    subroutine dealloc_utmp
      u=utmp
      deallocate( utmp )
    end subroutine dealloc_utmp

!    subroutine alloc_stmp
!      implicit none
!      allocate( stmp(i0:i1,i0:i1) ); stmp=S(i0:i1,i0:i1)
!    end subroutine alloc_stmp
!    subroutine dealloc_stmp
!      deallocate( stmp )
!    end subroutine dealloc_stmp

  end subroutine calc_dtrmm4

end module trmm_module
