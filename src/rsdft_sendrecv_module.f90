module rsdft_sendrecv_module

  implicit none

  private
  public :: test_sendrecv

  integer :: n_opt, n_opt_h

contains

#ifdef test
  subroutine test_sendrecv( node_partition )
    implicit none
    integer,intent(in) :: node_partition(:)
    integer :: i,j,k,m,mt,n,i0,ierr,myrank,nprocs,npow,p
    integer :: irank,jrank,itag,ireq(100),nreq
    real(8),allocatable :: a(:), b(:), c(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,100)
    integer,allocatable :: proc_map(:,:,:)
    integer :: n1,n2,n3,i1,i2,i3,idir
    integer :: pos_sendto(3), pos_recvfrom(3), pos_myrank(3)

    !return

    call write_border( 0, " test_sendrecv_a(start)" )

    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)

    n1=node_partition(1)
    n2=node_partition(2)
    n3=node_partition(3)
    allocate( proc_map(0:n1-1,0:n2-1,0:n3-1) ); proc_map=0
    n=-1
    do i3=0,n3-1
    do i2=0,n2-1
    do i1=0,n1-1
       n=n+1
       proc_map(i1,i2,i3)=n
       if ( n == myrank ) pos_myrank(1:3)=(/ i1,i2,i3 /)
    end do
    end do
    end do

    npow=16

    n=2**npow
    allocate( a(n) ); a(:)=myrank+1
    allocate( b(n) ); b(:)=0.0d0

    ctmin=1.d100
    etmin=1.d100

    nreq=0

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call cpu_time(ct0); et0=mpi_wtime()

    do idir=1,6

       pos_sendto(1:3) = pos_myrank(1:3)
       pos_recvfrom(1:3) = pos_myrank(1:3)
       select case( idir )
       case( 1 )
          pos_sendto(1) = mod( pos_myrank(1)+1+n1, n1 )
          pos_recvfrom(1) = mod( pos_myrank(1)-1+n1, n1 )
          m=9588
       case( 2 )
          pos_sendto(1) = mod( pos_myrank(1)-1+n1, n1 )
          pos_recvfrom(1) = mod( pos_myrank(1)+1+n1, n1 )
          m=9588
       case( 3 )
          pos_sendto(2) = mod( pos_myrank(2)+1+n2, n2 )
          pos_recvfrom(2) = mod( pos_myrank(2)-1+n2, n2 )
          m=15792
       case( 4 )
          pos_sendto(2) = mod( pos_myrank(2)-1+n2, n2 )
          pos_recvfrom(2) = mod( pos_myrank(2)+1+n2, n2 )
          m=15792
       case( 5 )
          pos_sendto(3) = mod( pos_myrank(3)+1+n3, n3 )
          pos_recvfrom(3) = mod( pos_myrank(3)-1+n3, n3 )
          m=11424
       case( 6 )
          pos_sendto(3) = mod( pos_myrank(3)-1+n3, n3 )
          pos_recvfrom(3) = mod( pos_myrank(3)+1+n3, n3 )
          m=11424
       end select

       irank = proc_map( pos_sendto(1)  , pos_sendto(2)  , pos_sendto(3)   )
       jrank = proc_map( pos_recvfrom(1), pos_recvfrom(2), pos_recvfrom(3) )
       itag = idir*10

       nreq = nreq + 1
       call MPI_Isend(a,2*m,MPI_REAL8,irank,itag,MPI_COMM_WORLD,ireq(nreq),ierr)

       !write(100+myrank,'(1x,20i6)') idir,1,m,irank,jrank,myrank,nreq

       nreq = nreq + 1
       call MPI_Irecv(b,2*m,MPI_REAL8,jrank,itag,MPI_COMM_WORLD,ireq(nreq),ierr)

       !write(200+myrank,'(1x,20i6)') idir,1,m,irank,jrank,myrank,nreq

    end do ! idir

    call MPI_Waitall( nreq, ireq, istatus, ierr )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call cpu_time(ct1); et1=mpi_wtime()

    write(*,*) "time=", ct1-ct0, et1-et0

    deallocate( b )
    deallocate( a )

    call write_border( 0, " test_sendrecv_a(end)" )

  end subroutine test_sendrecv

#else

  subroutine test_sendrecv
    implicit none
    integer :: i,j,k,m,mt,n,i0,ierr,myrank,nprocs,npow,p
    integer :: irank,jrank,itag,ireq(100),nreq
    real(8),allocatable :: a(:), b(:), c(:)
    real(8) :: ct,ct0,ct1,ctmin,et,et0,et1,etmin
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,100)

    !return

    call write_border( 0, " test_sendrecv(start)" )

    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)

    npow=18

    n=2**npow
    allocate( a(n) ); a(:)=myrank+1
    allocate( b(n) ); b(:)=0.0d0
    ctmin=1.d100
    etmin=1.d100
    do j=npow,0,-1
       m = 2**j
       mt=0
       do i=1,2**(npow-j)
          mt=mt+m
       end do
       b=0.0d0
       itag=0
       !irank=0
       !jrank=1
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct0) ; et0=mpi_wtime()
       do k=1,10
          do i=1,2**(npow-j)
             i0=(i-1)*m+1
             nreq=0
             do irank=0,nprocs-1,12
                jrank=mod(irank+12+nprocs,nprocs)
                !jrank=mod(myrank+12+nprocs,nprocs)
                if ( irank == myrank ) then
                   nreq=nreq+1
                   call MPI_Isend(a(i0),m,MPI_REAL8,jrank,itag,MPI_COMM_WORLD,ireq(nreq),ierr)
                end if
                !irank=mod(myrank-12+nprocs,nprocs)
                if ( jrank == myrank ) then
                   nreq=nreq+1
                   call MPI_Irecv(b(i0),m,MPI_REAL8,irank,itag,MPI_COMM_WORLD,ireq(nreq),ierr)
                end if
                !if ( nreq > 0 ) then
                !   call MPI_Waitall( nreq, ireq, istatus, ierr )
                !end if
             end do ! irank
             if ( nreq > 0 ) then
                call MPI_Waitall( nreq, ireq, istatus, ierr )
             end if
          end do
       end do ! k
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       call cpu_time(ct1) ; et1=mpi_wtime()
       ct=ct1-ct0
       et=et1-et0
       ctmin=min(ct,ctmin)
       if ( et < etmin ) then
          n_opt = m
          etmin = et
       end if
       if ( myrank == 0 ) then
          write(*,'(1x,3i12,4f10.5,2x,i8,g15.8)') m,i-1,mt,ct,ctmin,et,etmin,count(b/=0.0d0),sum(b)
       end if
    end do ! j
    deallocate( b )
    deallocate( a )
    call MPI_BCAST(n_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_opt_h=n_opt/2
    if ( myrank == 0 ) write(*,*) "n_opt, n_opt_h=",n_opt,n_opt_h

    call write_border( 0, " test_sendrecv(end)" )

  end subroutine test_sendrecv
#endif

end module rsdft_sendrecv_module
