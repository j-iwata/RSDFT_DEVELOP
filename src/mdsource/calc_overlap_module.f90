MODULE calc_overlap_module

  use parallel_module
  use watch_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: calc_overlap
  PUBLIC :: calc_overlap_no_mpi

  integer :: nblk0
  integer :: nblk1
  real(8),parameter :: alpha=1.d0,beta=0.d0

CONTAINS

  SUBROUTINE calc_overlap(m,n,a,b,dv,c)
    implicit none
    integer,intent(IN)  :: m,n
    real(8),intent(IN)  :: a(m,n),b(m,n),dv
    real(8),intent(OUT) :: c(n,n)
    integer :: nme,i,j,k,ierr
    real(8),allocatable :: s(:),r(:)
    real(8) :: ttmp(2),tttt(2,9)

    tttt=0.0d0
    !call watchb( ttmp )

    nblk0 = n
    nblk1 = 4

    call calc_overlap_sub(m,n,nblk0,a,b,c)

    !call watchb( ttmp, tttt(:,1) )

    nme = n*(n+1)/2
    allocate( s(nme),r(nme) )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       s(k)=c(i,j)*dv
    end do
    end do

    !call watchb( ttmp, tttt(:,2) )

    !call mpi_allreduce(MPI_IN_PLACE,s,nme,mpi_real8,mpi_sum,comm_grid,ierr)
    call rsdft_allreduce_sum( s, comm_grid )

    !call watchb( ttmp, tttt(:,3) )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       c(i,j)=s(k)
    end do
    end do

    deallocate( r,s )

    !call watchb( ttmp, tttt(:,4) )

!    if ( myrank == 0 ) then
!       do i=1,4
!          write(*,'(4x,"time_calc_overlap(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

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
    real(8) :: ttmp(2),tttt(2,9)
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
