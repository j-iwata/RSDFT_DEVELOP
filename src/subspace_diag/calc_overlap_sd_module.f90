module calc_overlap_sd_module

  use parallel_module, only: comm_grid
  use rsdft_mpi_module, only: rsdft_allreduce
!  use transpose_module, only: rsdft_transpose
!  use watch_module

  implicit none

  private
  public :: calc_overlap_sd
  !public :: calc_overlap_sd_blk
  !public :: calc_overlap_sd_dsyr2k

  integer :: nblk0
  integer :: nblk1
  !real(8),parameter :: alpha=1.0d0, beta=0.0d0

  integer :: nblk_ovlp=2
  !integer,public :: nblk_ovlp=2

  interface calc_overlap_sd
    module procedure d_calc_overlap_sd, z_calc_overlap_sd
  end interface

contains

  subroutine d_calc_overlap_sd( a,b,dv,c )
    implicit none
    real(8),intent(in)    :: a(:,:),b(:,:)
    real(8),intent(in)    :: dv
    real(8),intent(inout) :: c(:,:)
    integer :: nme,i,j,k,m,n,ierr
    real(8),allocatable :: s(:)
    real(8),save,allocatable :: c_tmp(:,:), u(:,:)
    real(8) :: ttmp(2),tttt(2,5)

    call write_border( 1, "d_calc_overlap_sd(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m = size( a, 1 )
    n = size( a, 2 )

    nblk0 = n
    nblk1 = nblk_ovlp; if ( nblk1 < 0 ) nblk1=n

    call d_calc_overlap_sd_sub(nblk0,a,b,c)

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

    call rsdft_allreduce( s, comm_grid )

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
!    call rsdft_transpose( c, 64 ) 
!
! ------

    !call watchb( ttmp, tttt(:,5), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,5
!          write(*,'(4x,"time_calc_overlap(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "calc_overlap_sd(end)" )

  end subroutine d_calc_overlap_sd

  recursive subroutine d_calc_overlap_sd_sub(nblk,a,b,c)
    implicit none
    integer,intent(in)  :: nblk
    real(8),intent(in)  :: a(:,:),b(:,:)
    real(8),intent(out) :: c(:,:)
    integer :: i,j,i0,i1,j0,j1,ni,nj,m,n,nblkh
    real(8),allocatable :: ctmp(:,:)
    real(8),parameter :: one=1.0d0, zero=0.0d0

    m = size( a, 1 )
    n = size( a, 2 )

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
                  ('T','N',ni,nj,m,one,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)

          else

             if ( ni > nblk1 ) then
                allocate( ctmp(ni,nj) ); ctmp=zero
                nblkh=nblk/2
                call d_calc_overlap_sd_sub(nblkh,a(:,i0:i1),b(:,j0:j1),ctmp)
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
!                do j=j0,j1
!                do i=j ,i1
!                    c(i,j)=sum( a(:,i)*b(:,j) )
!                end do
!                end do
                call dgemm &
                     ('T','N',ni,nj,m,one,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)
             end if

          end if

       end do ! j0

    end do ! i0

  end subroutine d_calc_overlap_sd_sub


  subroutine calc_overlap_sd_blk( a, b, dv, c )
    implicit none
    real(8),intent(in) :: a(:,:), b(:,:)
    real(8),intent(in) :: dv
    real(8),intent(inout) :: c(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: s(:)
    integer :: m,n,nblk,iblk,jblk
    integer :: i0,i1,j0,j1,ni,nj,nme,i,j,k

    m=size(a,1)
    n=size(a,2)
    nblk = n/nblk_ovlp

    do jblk = 1, n, nblk

       j0 = jblk
       j1 = min( jblk+nblk, n )
       nj = j1 - j0 + 1

       do iblk = jblk, n, nblk

          i0 = iblk
          i1 = min( iblk+nblk, n )
          ni = i1 - i0 + 1

          if ( iblk == jblk ) then
             call DGEMM( 'T','N',ni,nj,m,dv,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n )
          else
             call DGEMM( 'T','N',ni,nj,m,dv,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n )
          end if

       end do !iblk

    end do !jblk

    nme = n*(n+1)/2
    allocate( s(nme) )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       s(k)=c(i,j)*dv
    end do
    end do

    !call watchb( ttmp, tttt(:,2), barrier="on" )

    call rsdft_allreduce( s, comm_grid )

    !call watchb( ttmp, tttt(:,3), barrier="on" )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       c(i,j)=s(k)
    end do
    end do

    deallocate( s )

  end subroutine calc_overlap_sd_blk

  subroutine calc_overlap_sd_dsyr2k( a,b,dv,c )
    implicit none
    real(8),intent(in) :: a(:,:),b(:,:)
    real(8),intent(inout) :: c(:,:)
    real(8),intent(in) :: dv
    real(8),parameter :: zero=0.0d0
    real(8) :: z
    integer :: m,n
    m=size(a,1)
    n=size(a,2)
    z=dv*0.5d0
    call dsyr2k('L','T',n,m,z,a,m,b,m,zero,c,n)
  end subroutine calc_overlap_sd_dsyr2k


  subroutine z_calc_overlap_sd( a,b,dv,c )
    implicit none
    complex(8),intent(in)    :: a(:,:),b(:,:)
    real(8),intent(in)       :: dv
    complex(8),intent(inout) :: c(:,:)
    integer :: nme,i,j,k,m,n,ierr
    complex(8),allocatable :: s(:)
    complex(8),save,allocatable :: c_tmp(:,:), u(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    real(8) :: ttmp(2),tttt(2,5)

    call write_border( 1, "z_calc_overlap_sd(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m = size( a, 1 )
    n = size( a, 2 )

    nblk0 = n
    nblk1 = nblk_ovlp; if ( nblk1 < 0 ) nblk1=n

    call z_calc_overlap_sd_sub(nblk0,a,b,c)

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

    call rsdft_allreduce( s, comm_grid )

    !call watchb( ttmp, tttt(:,3), barrier="on" )

    do j=1,n
    do i=j,n
       k=(j-1)*n-(j*(j-1))/2+i
       c(i,j)=s(k)
    end do
    end do

    deallocate( s )

    !call watchb( ttmp, tttt(:,4), barrier="on" )

    !call write_border( 1, "calc_overlap_sd(end)" )

  end subroutine z_calc_overlap_sd

  recursive subroutine z_calc_overlap_sd_sub(nblk,a,b,c)
    implicit none
    integer,intent(in)  :: nblk
    complex(8),intent(in)  :: a(:,:),b(:,:)
    complex(8),intent(out) :: c(:,:)
    complex(8),allocatable :: ctmp(:,:)
    complex(8),parameter :: one=(1.0d0,0.0d0), zero=(0.0d0,0.0d0)
    integer :: i,j,i0,i1,j0,j1,m,n,ni,nj,nblkh

    m = size( a, 1 )
    n = size( a, 2 )

    do i0=1,n,nblk
       i1=min(i0+nblk-1,n)
       ni=i1-i0+1

       do j0=1,n,nblk
          j1=min(j0+nblk-1,n)
          nj=j1-j0+1

          if ( j0 > i1 ) then

             cycle

          else if ( j1 <= i0 ) then

             call zgemm &
                  ('C','N',ni,nj,m,one,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)

          else

             if ( ni > nblk1 ) then
                allocate( ctmp(ni,nj) ); ctmp=zero
                nblkh=nblk/2
                call z_calc_overlap_sd_sub(nblkh,a(:,i0:i1),b(:,j0:j1),ctmp)
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
!                do j=j0,j1
!                do i=j ,i1
!                    c(i,j)=sum( conjg(a(:,i))*b(:,j) )
!                end do
!                end do
                call zgemm &
                     ('C','N',ni,nj,m,one,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)
             end if

          end if

       end do ! j0

    end do ! i0

  end subroutine z_calc_overlap_sd_sub


end module calc_overlap_sd_module
