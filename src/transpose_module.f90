module transpose_module

  implicit none

  private
  public :: rsdft_transpose
  public :: estimate_transpose

contains

  subroutine rsdft_transpose( a, nblk )

    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,intent(in) :: nblk
    integer :: n,i,j,i0,i1,j0,j1,ii,jj
    real(8),allocatable :: blk(:,:),trp(:,:)

    n=size(a,1)
    allocate( blk(nblk,nblk) ); blk=0.0d0
    allocate( trp(nblk,nblk) ); trp=0.0d0

    do j0=1,n,nblk
       j1=min(j0+nblk-1,n)
       jj=j1-j0+1

       do i0=j0,n,nblk
          i1=min(i0+nblk-1,n)
          ii=i1-i0+1

          if ( i0 == j0 ) then

             do j=j0,j1
                do i=j+1,i1
                   blk(i-i0+1,j-j0+1)=a(i,j)
                end do
             end do
             trp(1:jj,1:ii)=transpose( blk(1:ii,1:jj) )
             do i=i0,i1
                do j=j0,i-1
                   a(j,i)=trp(j-j0+1,i-i0+1)
                end do
             end do

          else

             do j=j0,j1
                do i=i0,i1
                   blk(i-i0+1,j-j0+1)=a(i,j)
                end do
             end do
             trp(1:jj,1:ii)=transpose( blk(1:ii,1:jj) )
             do i=i0,i1
                do j=j0,j1
                   a(j,i)=trp(j-j0+1,i-i0+1)
                end do
             end do

          end if

       end do ! i0

    end do ! j0

    deallocate( trp )
    deallocate( blk )

  end subroutine rsdft_transpose


  subroutine estimate_transpose( n, nblk )
    implicit none
    integer,intent(in)  :: n
    integer,intent(out) :: nblk
    real(8),allocatable :: b(:,:)
    integer :: itry, it_best(3), i, j
    real(8) :: ct0,ct1,ct,ct_best(3)

    allocate( b(n,n) ); b=0.0d0

    nblk = 1
    ct_best(:) = 1.d100
    it_best(:) = 0

    do itry=1,10

       nblk = nblk*2
       if ( nblk >= n ) exit

       call cpu_time( ct0 )
       call rsdft_transpose( b, nblk )
       call cpu_time( ct1 )

       ct=ct1-ct0
       do i=1,3
          if ( ct < ct_best(i) ) then
             do j=3,i+1,-1
                ct_best(j) = ct_best(j-1)
                it_best(j) = it_best(j-1)
             end do
             ct_best(i) = ct
             it_best(i) = itry
          end if
       end do ! i

    end do ! itry

!    write(*,'(1x,3(f10.5,i10))') (ct_best(i), it_best(i), i=1,3)

    nblk = 2**it_best(1)

    deallocate( b )

  end subroutine estimate_transpose


end module transpose_module
