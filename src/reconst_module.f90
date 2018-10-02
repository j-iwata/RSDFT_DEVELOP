MODULE reconst_module

  use vector_tools_module

  implicit none

  private
  public :: reconstruct_init, reconstruct_scatter, reconstruct_gather
  public :: reconstruct_allgatherv

  integer :: comm            ! MPI
  integer :: nprocs, myrank  ! MPI
  integer :: ierr            ! MPI

CONTAINS


  SUBROUTINE reconstruct_init(commIN)
     implicit none 
     integer, intent(IN) :: commIN
     comm = commIN
     call mpi_comm_size(comm,nprocs,ierr)
     call mpi_comm_rank(comm,myrank,ierr)
  END SUBROUTINE reconstruct_init
  

  SUBROUTINE reconstruct_scatter( nsize1, nsize2, nsize3, nblock, u )
! nsize1: size of 1st-dimension
! nsize2: size of 2nd-dimension
! nsize3: size of 2nd-dimension
! nblock: block size
!
    implicit none
    integer,intent(IN) :: nsize1,nsize2,nsize3
    integer,intent(INOUT) :: nblock
    real(8),allocatable,intent(INOUT) :: u(:,:)
    real(8),allocatable :: utmp(:,:)
    integer :: i,j,k,ncycle,nsize3_min,nsize3_max

    call write_border( 1, "reconstruct_scatter(start)" )

    allocate( utmp(nsize1,nsize3) ) ; utmp=0.0d0

    if ( nblock <= 0 .or. nblock >= nsize2 .or. nblock == nsize3 ) then

       ncycle = nprocs
       nblock = nsize3

    else

       ncycle = nsize2/nblock

       if ( mod(ncycle,nprocs) /= 0 ) then
          write(*,*) "load imbalance is occurred"
          nsize3_min = ncycle/nprocs * nblock
          nsize3_max = (ncycle/nprocs+1) * nblock
          write(*,*) "nsize3_min,nsize3_max,nsize3=",nsize3_min,nsize3_max,nsize3
          if ( nsize3_max > nsize3 ) call stop_program( "stop@reconstruct_scatter" )
       end if

    end if

    k=0
    do j=1,ncycle
       if ( mod(j-1,nprocs) == myrank ) then
          do i=1,nblock
             k=k+1
             utmp(:,k)=u(:,nblock*(j-1)+i)
!             write(*,*) myrank,nblock*(j-1)+i,k,nblock,nprocs
          end do
       end if
    end do

    deallocate( u )
    allocate( u(nsize1,nsize3) ) ; u=0.0d0

    u = utmp
    deallocate( utmp )

    call write_border( 1, "reconstruct_scatterv(end)" )

  END SUBROUTINE reconstruct_scatter


  SUBROUTINE reconstruct_gather( nsize1, nsize2, nsize3, nblock, u )
! nsize1: size of 1st-dimension
! nsize2: size of 2nd-dimension
! nsize3: size of 2nd-dimension  !MB
! nblock: block size

    implicit none
    integer,intent(IN) :: nsize1, nsize2, nsize3, nblock
    real(8),allocatable,intent(INOUT) :: u(:,:)
    real(8),allocatable :: utmp1(:,:),utmp2(:,:)
    integer :: ntmp_size, iroot, ierr
    integer :: i, j, k, ncycle
    include 'mpif.h'

    call write_border( 1, "reconstruct_gather(start)" )

    allocate( utmp1(nsize1,nsize2) ) ; utmp1=0.0d0  ! for save
    allocate( utmp2(nsize1,nsize2) ) ; utmp2=0.0d0  ! for bcast
    ntmp_size=nsize1*nsize2

    utmp1=u

    deallocate( u )
    allocate( u(nsize1,nsize3) ) ; u=0.0d0

    ncycle=nsize2/nblock

    do k=1,nprocs
       iroot=k-1
       utmp2=0.0d0
       if ( iroot == myrank ) utmp2=utmp1
       call MPI_BCAST( utmp2, ntmp_size, MPI_REAL8, iroot, comm, ierr )
       do j=1,ncycle
          do i=1,nblock
             u(:,nprocs*nblock*(j-1)+nblock*(k-1)+i)=utmp2(:,nblock*(j-1)+i)
          end do
       end do
    end do

    deallocate( utmp1, utmp2 )

    call write_border( 1, "reconstruct_gather(end)" )

  END SUBROUTINE reconstruct_gather


  SUBROUTINE reconstruct_allgatherv( u, v )
    implicit none
    real(8),allocatable,intent(INOUT) :: u(:,:)
    type(vinfo),intent(IN) :: v(2)
    real(8),allocatable :: utmp(:,:)
    integer :: m0,n0,nn,np,i
    integer,allocatable :: ir(:),id(:)
    include 'mpif.h'
    m0=size( u, 1 )
    n0=size( u, 2 )
    nn=sum( v(2)%pinfo%ir )
    np=v(2)%pinfo%np
    allocate( utmp(m0,n0) ) ; utmp=0.0d0
    utmp=u
    deallocate( u )
    allocate( u(m0,nn)   ) ; u=0.0d0
    allocate( ir(0:np-1) ) ; ir=0
    allocate( id(0:np-1) ) ; id=0
    ir(:)=m0*v(2)%pinfo%ir(:)
    id(:)=m0*v(2)%pinfo%id(:)
    call MPI_ALLGATHERV( utmp, size(utmp), MPI_REAL8, u, ir, id, MPI_REAL8, v(2)%pinfo%comm, i )
    deallocate( id   )
    deallocate( ir   )
    deallocate( utmp )
  END SUBROUTINE reconstruct_allgatherv


END MODULE reconst_module
