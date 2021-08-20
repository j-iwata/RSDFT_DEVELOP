module ps_nloc2_comm_module

  implicit none
  private
  public :: d_ps_nloc2_comm
  public :: z_ps_nloc2_comm

contains

  subroutine d_ps_nloc2_comm( w )
    use ps_nloc2_variables, only: nrlma_xyz, num_2_rank, lma_nsend, &
    sendmap, recvmap, d_sbufnl3, d_rbufnl3, z_sbufnl3, z_rbufnl3
    use parallel_module, only: comm_grid
    implicit none
    real(8),intent(inout) :: w(:,:,:)
    real(8),allocatable :: v(:,:,:)
    integer :: i,j,m,n1,n2,n3,i0,i1,i2,i3
    integer :: irank,jrank,nreq,ierr
    integer,allocatable :: istatus(:,:), ireq(:)
    include 'mpif.h'

    n1 = size( w, 1 )
    n2 = size( w, 2 )
    n3 = size( w, 3 )

    allocate( v(n1,n2,n3) ); v=0.0d0

    m = max( 6, 6*product(nrlma_xyz) ) * 2
    allocate( istatus(MPI_STATUS_SIZE,m) ); istatus=0
    allocate( ireq(m) ); ireq=0
    
    do i = 1, 6

      select case(i)
      case(1,3,5)
        j=i+1
        v(:,:,:)=w(:,:,:)
      case(2,4,6)
        j=i-1
      end select
      
      do m = 1, nrlma_xyz(i)
        nreq=0
        irank=num_2_rank(m,i)
        jrank=num_2_rank(m,j)
        if( irank >= 0 )then
          i0=0
          do i3 = 1, n3
          do i2 = 1, lma_nsend(irank)
          do i1 = 1, n1
            i0=i0+1
            d_sbufnl3(i0,m,i) = v( i1, sendmap(i2,irank), i3 )
          end do
          end do
          end do
          nreq=nreq+1
          call MPI_Isend(d_sbufnl3(1,m,i),n1*lma_nsend(irank)*n3, &
          MPI_REAL8,irank,1,comm_grid,ireq(nreq),ierr)
        end if
        if( jrank >= 0 )then
          nreq=nreq+1
          call MPI_Irecv(d_rbufnl3(1,m,i),n1*lma_nsend(jrank)*n3, &
          MPI_REAL8,jrank,1,comm_grid,ireq(nreq),ierr)
        end if
        call MPI_Waitall(nreq,ireq,istatus,ierr)
        if( jrank >= 0 )then
          i0=0
          do i3 = 1, n3
          do i2 = 1, lma_nsend(jrank)
          do i1 = 1, n1
            i0=i0+1
            w(i1,recvmap(i2,jrank),i3) = w(i1,recvmap(i2,jrank),i3) + d_rbufnl3(i0,m,i)
          end do
          end do
          end do
        end if
      end do ! m
    end do ! i

    deallocate( ireq )
    deallocate( istatus )
    deallocate( v )

  end subroutine d_ps_nloc2_comm

  subroutine z_ps_nloc2_comm( w )
    use ps_nloc2_variables, only: nrlma_xyz, num_2_rank, lma_nsend, &
    sendmap, recvmap, z_sbufnl3, z_rbufnl3
    use parallel_module, only: comm_grid
    implicit none
    complex(8),intent(inout) :: w(:,:,:)
    complex(8),allocatable :: v(:,:,:)
    integer :: i,j,m,n1,n2,n3,i0,i1,i2,i3
    integer :: irank,jrank,nreq,ierr
    integer,allocatable :: istatus(:,:), ireq(:)
    include 'mpif.h'

    n1 = size( w, 1 )
    n2 = size( w, 2 )
    n3 = size( w, 3 )

    allocate( v(n1,n2,n3) ); v=(0.0d0,0.0d0)

    m = max( 6, 6*product(nrlma_xyz) ) * 2
    allocate( istatus(MPI_STATUS_SIZE,m) ); istatus=0
    allocate( ireq(m) ); ireq=0
    
    do i = 1, 6

      select case(i)
      case(1,3,5)
        j=i+1
        v(:,:,:)=w(:,:,:)
      case(2,4,6)
        j=i-1
      end select
      
      do m = 1, nrlma_xyz(i)
        nreq=0
        irank=num_2_rank(m,i)
        jrank=num_2_rank(m,j)
        if( irank >= 0 )then
          i0=0
          do i3 = 1, n3
          do i2 = 1, lma_nsend(irank)
          do i1 = 1, n1
            i0=i0+1
            z_sbufnl3(i0,m,i) = v( i1, sendmap(i2,irank), i3 )
          end do
          end do
          end do
          nreq=nreq+1
          call MPI_Isend(z_sbufnl3(1,m,i),n1*lma_nsend(irank)*n3, &
          MPI_COMPLEX16,irank,1,comm_grid,ireq(nreq),ierr)
        end if
        if( jrank >= 0 )then
          nreq=nreq+1
          call MPI_Irecv(z_rbufnl3(1,m,i),n1*lma_nsend(jrank)*n3, &
          MPI_COMPLEX16,jrank,1,comm_grid,ireq(nreq),ierr)
        end if
        call MPI_Waitall(nreq,ireq,istatus,ierr)
        if( jrank >= 0 )then
          i0=0
          do i3 = 1, n3
          do i2 = 1, lma_nsend(jrank)
          do i1 = 1, n1
            i0=i0+1
            w(i1,recvmap(i2,jrank),i3) = w(i1,recvmap(i2,jrank),i3) + z_rbufnl3(i0,m,i)
          end do
          end do
          end do
        end if
      end do ! m
    end do ! i

    deallocate( ireq )
    deallocate( istatus )
    deallocate( v )

  end subroutine z_ps_nloc2_comm

end module ps_nloc2_comm_module
