MODULE sl_tools_module

  implicit none

  PRIVATE
  PUBLIC :: fill_matrix
  PUBLIC :: gather_matrix
  PUBLIC :: gatherA_matrix
  PUBLIC :: distribute_matrix
  PUBLIC :: slinfo
  PUBLIC :: sl_block_map

  complex(8),parameter :: zero=(0.0d0,0.0d0)

  type slinfo
     integer :: myrow,  mycol
     integer :: nprow,  npcol
     integer :: mbsize, nbsize
     integer,allocatable :: map_1to2(:,:)
     integer,allocatable :: map_2to1(:,:)
     integer :: icontxt_a, icontxt_b
     integer :: desca(9), descb(9)
  end type slinfo

  type blkinfo
     integer :: iband(2), jband(2)
     integer :: istor(2), jstor(2)
     integer :: iproc, jproc
  end type blkinfo

  type(blkinfo),allocatable :: blk_map(:,:)

  INTERFACE fill_matrix
     MODULE PROCEDURE d_fill_matrix, z_fill_matrix
  END INTERFACE

  INTERFACE gather_matrix
     MODULE PROCEDURE d_gather_matrix, z_gather_matrix
  END INTERFACE

  INTERFACE gatherA_matrix
     MODULE PROCEDURE d_gatherA_matrix
  END INTERFACE

  INTERFACE distribute_matrix
     MODULE PROCEDURE d_distribute_matrix, z_distribute_matrix
  END INTERFACE

CONTAINS


  SUBROUTINE d_fill_matrix( h )
    implicit none
    real(8),intent(INOUT) :: h(:,:)
    integer :: n,i,j
    n=size(h,1)
    do j=1,n
       do i=j+1,n
          h(i,j) = h(j,i)
       end do
    end do
  END SUBROUTINE d_fill_matrix

  SUBROUTINE z_fill_matrix( h )
    implicit none
    complex(8),intent(INOUT) :: h(:,:)
    integer :: n,i,j
    n=size(h,1)
    do j=1,n
       do i=j+1,n
          h(i,j) = conjg( h(j,i) )
       end do
    end do
  END SUBROUTINE z_fill_matrix


  SUBROUTINE d_distribute_matrix( sl, H, Hsub )
    implicit none
    type(slinfo),intent(IN) :: sl
    real(8),intent(IN) :: H(:,:)
    real(8),intent(OUT) :: Hsub(:,:)
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    m = size( H, 1 )
    n = size( H, 2 )
    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
             Hsub(is+i,js+j) = H(i0+i-1,j0+j-1)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do
  END SUBROUTINE d_distribute_matrix

  SUBROUTINE z_distribute_matrix( sl, H, Hsub )
    implicit none
    type(slinfo),intent(IN) :: sl
    complex(8),intent(IN) :: H(:,:)
    complex(8),intent(OUT) :: Hsub(:,:)
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    m = size( H, 1 )
    n = size( H, 2 )
    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
             Hsub(is+i,js+j) = H(i0+i-1,j0+j-1)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do
  END SUBROUTINE z_distribute_matrix


  SUBROUTINE d_gather_matrix( sl, Hsub, H, icontxt )
    implicit none
    type(slinfo) :: sl
    real(8),intent(IN) :: Hsub(:,:)
    real(8),intent(OUT) :: H(:,:)
    integer,optional,intent(IN) :: icontxt
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    include 'mpif.h'

    H(:,:) = 0.0d0

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
             H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do

    if ( present(icontxt) ) then
       call dgsum2d( icontxt, 'All', ' ', m, n, H, m, -1, -1 )
    else
       call mpi_allreduce( mpi_in_place, H, size(H), mpi_real8 &
            ,mpi_sum, mpi_comm_world, i )
    end if

  END SUBROUTINE d_gather_matrix

  SUBROUTINE z_gather_matrix( sl, Hsub, H )
    implicit none
    type(slinfo) :: sl
    complex(8),intent(IN) :: Hsub(:,:)
    complex(8),intent(OUT) :: H(:,:)
    integer :: m,n,m0,n0,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    include 'mpif.h'

    H(:,:) = zero

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
             H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do

    call mpi_allreduce( mpi_in_place, H, size(H), mpi_complex16 &
         ,mpi_sum, mpi_comm_world, i )

  END SUBROUTINE z_gather_matrix

  SUBROUTINE d_gatherA_matrix( sl, Hsub, H, icontxt )
    implicit none
    type(slinfo) :: sl
    real(8),intent(IN) :: Hsub(:,:)
    real(8),intent(OUT) :: H(:,:)
    integer,optional,intent(IN) :: icontxt
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
    include 'mpif.h'

    H(:,:) = 0.0d0

    m  = size( H, 1 )
    n  = size( H, 2 )

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
!            if ( i0+i-1 < j0+j-1 ) cycle
             H(i0+i-1,j0+j-1) = Hsub(is+i,js+j)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do

    if ( present(icontxt) ) then
       call dgsum2d( icontxt, 'All', ' ', m, n, H, m, -1, -1 )
    else
       call mpi_allreduce( mpi_in_place, H, size(H), mpi_real8 &
            ,mpi_sum, mpi_comm_world, i )
    end if

  END SUBROUTINE d_gatherA_matrix


  SUBROUTINE sl_block_map( m,n,sl )
    implicit none
    integer,intent(IN) :: m,n
    type(slinfo),intent(IN) :: sl
    integer :: i0,j0,i1,j1,ip,jp,ii,jj
    integer :: ib,jb,mblk,nblk
    integer,allocatable :: is(:,:),js(:,:)

    mblk = (m+sl%mbsize-1)/sl%mbsize
    nblk = (n+sl%nbsize-1)/sl%nbsize

    if ( allocated(blk_map) ) deallocate(blk_map)
    allocate( blk_map(mblk,nblk) )

    allocate( is(0:sl%nprow-1,0:sl%npcol-1) ) ; is=0
    allocate( js(0:sl%nprow-1,0:sl%npcol-1) ) ; js=0

    jp=-1
    jb= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jb=jb+1
       jp=mod(jp+1,sl%npcol)
       ip=-1
       is(:,jp)=0
       ib= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ib=ib+1
          ip=mod(ip+1,sl%nprow)
          blk_map(ib,jb)%iband(1)=i0
          blk_map(ib,jb)%iband(2)=i1
          blk_map(ib,jb)%jband(1)=j0
          blk_map(ib,jb)%jband(2)=j1
          blk_map(ib,jb)%istor(1)=is(ip,jp)+1
          blk_map(ib,jb)%istor(2)=is(ip,jp)+ii
          blk_map(ib,jb)%jstor(1)=js(ip,jp)+1
          blk_map(ib,jb)%jstor(2)=js(ip,jp)+jj
          blk_map(ib,jb)%iproc=ip
          blk_map(ib,jb)%jproc=jp
          is(ip,jp)=is(ip,jp)+ii
       end do
       js(:,jp)=js(:,jp)+jj
    end do
    deallocate( js, is )
  END SUBROUTINE sl_block_map


END MODULE sl_tools_module
