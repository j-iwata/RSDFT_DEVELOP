MODULE calc_overlap_module

  use parallel_module
  use watch_module
  use rsdft_mpi_module
  use transpose_module

  implicit none

  PRIVATE
  PUBLIC :: calc_overlap
  PUBLIC :: calc_overlap_no_mpi

  integer :: nblk0
  integer :: nblk1
  real(8),parameter :: alpha=1.0d0, beta=0.0d0

CONTAINS

  SUBROUTINE calc_overlap(m,n,a,b,dv,c)
    implicit none
    integer,intent(IN)  :: m,n
    real(8),intent(IN)  :: a(m,n),b(m,n),dv
    real(8),intent(INOUT) :: c(n,n)
    integer :: nme,i,j,k,ierr
    real(8),allocatable :: s(:)
    real(8),save,allocatable :: c_tmp(:,:), u(:,:)
    real(8) :: ttmp(2),tttt(2,5)

    !call write_border( 1, "calc_overlap(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    nblk0 = n
    nblk1 = 4

    call calc_overlap_sub(m,n,nblk0,a,b,c)

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    nme = n*(n+1)/2
    allocate( s(nme) )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       s(k)=c(i,j)*dv
    end do
    end do

    !call watchb( ttmp, tttt(:,2), barrier="on" )

    call rsdft_allreduce_sum( s, comm_grid )

    !call watchb( ttmp, tttt(:,3), barrier="on" )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       c(i,j)=s(k)
    end do
    end do

    deallocate( s )

    !call watchb( ttmp, tttt(:,4), barrier="on" )

! ---(1)
!
!!$omp parallel do
!    do j=1,n-1
!       do i=j+1,n
!          c(j,i) = c(i,j)
!       end do
!    end do
!!$omp end parallel do
!
! ---(2)
!
!    if ( .not.allocated(c_tmp) ) then
!       allocate( c_tmp(n,n) ); c_tmp=0.0d0
!    end if
!    c_tmp=transpose(c)
!!$omp parallel do
!    do j=1,n
!       do i=j+1,n
!          c(i,j)=c_tmp(i,j)
!       end do
!    end do
!!$omp end parallel do
!
! ---(3)
!
!    if ( .not.allocated(u) ) then
!       allocate( u(n,n) ); u=0.0d0
!       do i=1,n
!          u(i,i)=1.0d0
!       end do
!    end if
!    do j=1,n-1
!       do i=j+1,n
!          c_tmp(i,j)=c(i,j)
!       end do
!    end do
!    call DGEMM('T','N',n,n,n,1.0d0,c_tmp,n,u,n,1.0d0,c,n)
!
! ---(4)
!
    call rsdft_transpose( c, 64 ) 
!
! ------

    !call watchb( ttmp, tttt(:,5), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,5
!          write(*,'(4x,"time_calc_overlap(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "calc_overlap(end)" )

  END SUBROUTINE calc_overlap

  RECURSIVE SUBROUTINE calc_overlap_sub(m,n,nblk,a,b,c)
    implicit none
    integer,intent(IN)  :: m,n,nblk
    real(8),intent(IN)  :: a(m,n),b(m,n)
    real(8),intent(OUT) :: c(n,n)
    integer :: i,j,i0,i1,j0,j1,ni,nj,nblkh
    real(8),allocatable :: ctmp(:,:)

    do i0=1,n,nblk
       i1=min(i0+nblk-1,n)
       ni=i1-i0+1

       do j0=1,n,nblk
          j1=min(j0+nblk-1,n)
          nj=j1-j0+1

          if ( j0 > i1 ) then

             cycle

          else if ( j1 <= i0 ) then

             call dgemm &
                  ('T','N',ni,nj,m,alpha,a(1,i0),m,b(1,j0),m,beta,c(i0,j0),n)

          else

             if ( ni > nblk1 ) then
                allocate( ctmp(ni,ni) ) ; ctmp=0.d0
                nblkh=nblk/2
                call calc_overlap_sub(m,ni,nblkh,a(1,i0),b(1,j0),ctmp)
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
                do j=j0,j1
                do i=j ,i1
                    c(i,j)=sum( a(:,i)*b(:,j) )
                end do
                end do
             end if

          end if

       end do ! j0

    end do ! i0

  END SUBROUTINE calc_overlap_sub


  SUBROUTINE calc_overlap_no_mpi( a, b, dv, c )
    implicit none
    real(8),intent(IN)  :: a(:,:),b(:,:),dv
    real(8),intent(OUT) :: c(:,:)
    real(8) :: ttmp(2),tttt(2,1)
    integer :: m,n

    !call watchb( ttmp ); tttt=0.0d0

    m = size( a, 1 )
    n = size( a, 2 )

    nblk0 = n
    nblk1 = 4

    call calc_overlap_sub( m, n, nblk0, a, b, c )

    c=c*dv

    !call watchb( ttmp, tttt(:,1) )

!    if ( myrank == 0 ) then
!       do i=1,1
!          write(*,'(4x,"time_calc_overlap_no_mpi(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

  END SUBROUTINE calc_overlap_no_mpi


END MODULE calc_overlap_module
