module rsdft_mpi_module

  implicit none

  private
  public :: init_rsdft_mpi
  public :: rsdft_allreduce
  public :: rsdft_allreduce_sum
  public :: d_rsdft_allreduce_sum_5, d_rsdft_allreduce_sum_6
  ! public :: rsdft_bcast
  public :: rsdft_allgather
  public :: rsdft_allgatherv

  interface rsdft_allreduce  
    module procedure d_rsdft_allreduce_0, &
                     i_rsdft_allreduce_0, &
                     d_rsdft_allreduce_1, z_rsdft_allreduce_1, &
                     d_rsdft_allreduce_2, z_rsdft_allreduce_2, &
                     i_rsdft_allreduce_2, &
                     z_rsdft_allreduce_3
  end interface

  INTERFACE rsdft_allreduce_sum
     MODULE PROCEDURE d_rsdft_allreduce_sum_0, z_rsdft_allreduce_sum_0, &
                      i_rsdft_allreduce_sum_1,                          &
                      d_rsdft_allreduce_sum_1, z_rsdft_allreduce_sum_1, &
                      i_rsdft_allreduce_sum_2,                          &
                      d_rsdft_allreduce_sum_2, z_rsdft_allreduce_sum_2, &
                      i_rsdft_allreduce_sum_3,                          &
                      d_rsdft_allreduce_sum_3, z_rsdft_allreduce_sum_3
  END INTERFACE

  ! INTERFACE rsdft_bcast
  !    MODULE PROCEDURE d_rsdft_bcast2, z_rsdft_bcast2
  ! END INTERFACE

  INTERFACE rsdft_allgather
     MODULE PROCEDURE l_rsdft_allgather12, &
                      i_rsdft_allgather01, &
                      i_rsdft_allgather12, &
                      d_rsdft_allgather12, &
                      d_rsdft_allgather23, &
                      d_rsdft_allgather34, &
                      d_rsdft_allgather11, &
                      d_rsdft_allgather22, &
                      z_rsdft_allgather23
  END INTERFACE

  INTERFACE rsdft_allgatherv
     MODULE PROCEDURE i_rsdft_allgatherv22, &
                      d_rsdft_allgatherv11, &
                      d_rsdft_allgatherv22, &
                      d_rsdft_allgatherv33, &
                      z_rsdft_allgatherv11, &
                      z_rsdft_allgatherv22, &
                      z_rsdft_allgatherv33
  END INTERFACE

  integer :: RSDFT_MPI_COMPLEX16
  integer :: ndiv_allreduce = 100

CONTAINS

!-------------------------------------------------- rsdft_allreduce_sum

  subroutine init_rsdft_mpi
    implicit none
#ifndef _NOMPI_
    include 'mpif.h'
#ifdef _NO_MPI_COMPLEX16_
    RSDFT_MPI_COMPLEX16=MPI_DOUBLE_COMPLEX
#else
    RSDFT_MPI_COMPLEX16=MPI_COMPLEX16
#endif
#endif
  end subroutine init_rsdft_mpi

  subroutine i_rsdft_allreduce_0( a, comm_in, op_in )
    implicit none
    integer,intent(inout) :: a
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, comm, op, ierr
    integer :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    a0=a
    m=1
    call MPI_Allreduce(a0,a,m,MPI_INTEGER,op,comm,ierr)
#endif
  end subroutine i_rsdft_allreduce_0

  subroutine d_rsdft_allreduce_0( a, comm_in, op_in )
    implicit none
    real(8),intent(inout) :: a
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, comm, op, ierr
    real(8) :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    a0=a
    m=1
    call MPI_Allreduce(a0,a,m,MPI_REAL8,op,comm,ierr)
#endif
  end subroutine d_rsdft_allreduce_0

  subroutine d_rsdft_allreduce_1( a, comm_in, op_in )
    implicit none
    real(8),intent(inout) :: a(:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, comm, op, ierr
    real(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m) ); a0=a
    call MPI_Allreduce( a0, a, m, MPI_REAL8, op, comm, ierr )
    deallocate( a0 )
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m,MPI_REAL8,op,comm,ierr)
#endif
#endif
  end subroutine d_rsdft_allreduce_1

  subroutine z_rsdft_allreduce_1( a, comm_in, op_in )
    implicit none
    complex(8),intent(inout) :: a(:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, comm, op, ierr
    complex(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m) ); a0=a
    call MPI_Allreduce( a0, a, m, RSDFT_MPI_COMPLEX16, op, comm, ierr )
    deallocate( a0 )
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m,RSDFT_MPI_COMPLEX16,op,comm,ierr)
#endif
#endif
  end subroutine z_rsdft_allreduce_1

  subroutine i_rsdft_allreduce_2( a, comm_in, op_in )
    implicit none
    integer,intent(inout) :: a(:,:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, n, comm, op, ierr
    integer,allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ); a0=a
    call MPI_Allreduce( a0, a, m*n, MPI_INTEGER, op, comm, ierr )
    deallocate( a0 )
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m*n,MPI_INTEGER,op,comm,ierr)
#endif
#endif
  end subroutine i_rsdft_allreduce_2

  subroutine d_rsdft_allreduce_2( a, comm_in, op_in )
    implicit none
    real(8),intent(inout) :: a(:,:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, n, comm, op, ierr
    real(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ); a0=a
    call MPI_Allreduce( a0, a, m*n, MPI_REAL8, op, comm, ierr )
    deallocate( a0 )
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m*n,MPI_REAL8,op,comm,ierr)
#endif
#endif
  end subroutine d_rsdft_allreduce_2

  subroutine z_rsdft_allreduce_2( a, comm_in, op_in )
    implicit none
    complex(8),intent(inout) :: a(:,:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    integer :: m, n, comm, op, ierr
    complex(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ); a0=a
    call MPI_Allreduce(a0,a,m*n,RSDFT_MPI_COMPLEX16,op,comm,ierr)
    deallocate( a0 )
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m*n,RSDFT_MPI_COMPLEX16,op,comm,ierr)
#endif
#endif
  end subroutine z_rsdft_allreduce_2

  subroutine z_rsdft_allreduce_3( a, comm_in, op_in, out )
    implicit none
    complex(8),intent(inout) :: a(:,:,:)
    integer,optional,intent(in) :: comm_in
    character(*),optional,intent(in) :: op_in
    complex(8),optional,intent(out) :: out(:,:,:)
    integer :: m1,m2,m3, comm, op, ierr
    complex(8),allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    if ( present(out) ) out=a
    return
#else
    include 'mpif.h'
    comm=MPI_COMM_WORLD; if ( present(comm_in) ) comm=comm_in
    op=MPI_SUM; if ( present(op_in) ) op=get_op_id(op_in)
    m1=size(a,1)
    m2=size(a,2)
    m3=size(a,3)
#ifdef _NO_MPI_INPLACE_
    if ( present(out) ) then
      call MPI_Allreduce(a,out,m1*m2*m3,RSDFT_MPI_COMPLEX16,op,comm,ierr)
    else
      allocate( a0(m1,m2,m3) ); a0=a
      call MPI_Allreduce(a0,a,m1*m2*m3,RSDFT_MPI_COMPLEX16,op,comm,ierr)
      deallocate( a0 )
    end if
#else
    call MPI_Allreduce(MPI_IN_PLACE,a,m1*m2*m3,RSDFT_MPI_COMPLEX16,op,comm,ierr)
#endif
#endif
  end subroutine z_rsdft_allreduce_3

  function get_op_id( op_keyword )
    implicit none
    integer :: get_op_id
    character(*),intent(in) :: op_keyword
#ifdef _NOMPI_
    get_op_id = 0
    return
#else
    include 'mpif.h'
    select case( op_keyword )
    case default ; get_op_id=MPI_SUM
    case( 'min' ); get_op_id=MPI_MIN
    case( 'max' ); get_op_id=MPI_MAX
    end select
#endif
  end function get_op_id


  SUBROUTINE d_rsdft_allreduce_sum_0( a, comm )
    implicit none
    real(8),intent(INOUT) :: a
    integer,intent(IN) :: comm
    integer :: ierr
    real(8) :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    a0=a
    call MPI_ALLREDUCE( a0, a, 1, MPI_REAL8, MPI_SUM, comm, ierr )
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_0

  SUBROUTINE z_rsdft_allreduce_sum_0( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a
    integer,intent(IN) :: comm
    integer :: ierr
    complex(8) :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    a0=a
    call MPI_ALLREDUCE( a0, a, 1, RSDFT_MPI_COMPLEX16, MPI_SUM, comm, ierr )
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_0


  SUBROUTINE i_rsdft_allreduce_sum_1( a, comm )
    implicit none
    integer,intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    integer,allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, MPI_INTEGER, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,MPI_INTEGER,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE i_rsdft_allreduce_sum_1

  SUBROUTINE d_rsdft_allreduce_sum_1( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    real(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_1

  SUBROUTINE z_rsdft_allreduce_sum_1( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    complex(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, RSDFT_MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,RSDFT_MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_1


  SUBROUTINE i_rsdft_allreduce_sum_2( a, comm )
    implicit none
    integer,intent(INOUT) :: a(:,:)
    integer,intent(IN) :: comm
    integer :: m, n, ierr
    integer,allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, m*n, MPI_INTEGER, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,m*n,MPI_INTEGER,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE i_rsdft_allreduce_sum_2

  SUBROUTINE d_rsdft_allreduce_sum_2( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: comm
    integer :: m, n, ierr
    real(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    m=size(a,1)
    n=size(a,2)
#ifdef _DIV_ALLREDUCE_
    call allreduce_sub( m*n, comm, a, ndiv_allreduce )
#else
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, m*n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,m*n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_2

  SUBROUTINE z_rsdft_allreduce_sum_2( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: comm
    integer :: m, n, ierr
    complex(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    m=size(a,1)
    n=size(a,2)
#ifdef _DIV_ALLREDUCE_
    call allreduce_sub( 2*m*n, comm, a, ndiv_allreduce )
#else
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, m*n, RSDFT_MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,m*n,RSDFT_MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_2


  SUBROUTINE i_rsdft_allreduce_sum_3( a, comm )
    implicit none
    integer,intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    integer,allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, MPI_INTEGER, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,MPI_INTEGER,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE i_rsdft_allreduce_sum_3

  SUBROUTINE d_rsdft_allreduce_sum_3( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    real(8),allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _DIV_ALLREDUCE_
    call allreduce_sub( l*m*n, comm, a, ndiv_allreduce )
#else
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_3

  SUBROUTINE z_rsdft_allreduce_sum_3( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    complex(8),allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _DIV_ALLREDUCE_
    call allreduce_sub( 2*l*m*n, comm, a, ndiv_allreduce )
#else
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, RSDFT_MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,RSDFT_MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_3


  SUBROUTINE d_rsdft_allreduce_sum_5( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:,:,:)
    integer,intent(IN) :: comm
    integer :: i, m(6), mm, ierr
    real(8),allocatable :: a0(:,:,:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    mm=1
    do i=1,5
       m(i)=size(a,i)
       mm=mm*m(i)
    end do
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m(1),m(2),m(3),m(4),m(5)) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, mm, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,mm,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_5


  SUBROUTINE d_rsdft_allreduce_sum_6( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:,:,:,:)
    integer,intent(IN) :: comm
    integer :: i, m(6), mm, ierr
    real(8),allocatable :: a0(:,:,:,:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    mm=1
    do i=1,6
       m(i)=size(a,i)
       mm=mm*m(i)
    end do
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m(1),m(2),m(3),m(4),m(5),m(6)) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, mm, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,mm,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_6

! !-------------------------------------------------- rsdft_bcast
!
!   SUBROUTINE d_rsdft_bcast2( a, src, comm )
!     implicit none
!     real(8),intent(INOUT) :: a(:,:)
!     integer,intent(IN) :: src, comm
!     integer :: ierr
! #ifdef _NOMPI_
!     return
! #else
!     include 'mpif.h'
!     call MPI_BCAST( a, size(a), MPI_REAL8, src, comm, ierr )
! #endif
!   END SUBROUTINE d_rsdft_bcast2
!
!   SUBROUTINE z_rsdft_bcast2( a, src, comm )
!     implicit none
!     complex(8),intent(INOUT) :: a(:,:)
!     integer,intent(IN) :: src, comm
!     integer :: ierr
! #ifdef _NOMPI_
!     return
! #else
!     include 'mpif.h'
!     call MPI_BCAST( a, size(a), RSDFT_MPI_COMPLEX16, src, comm, ierr )
! #endif
!   END SUBROUTINE z_rsdft_bcast2

!-------------------------------------------------- rsdft_allgather

  SUBROUTINE l_rsdft_allgather12( f, g, comm )
    implicit none
    logical,intent(INOUT) :: f(:), g(:,0:)
    integer,intent(IN) :: comm
    integer :: i
    logical,allocatable :: t(:)
#ifdef _NO_MPI_
    g(:,0) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f)) ) ; t=f
    call MPI_ALLGATHER( t,size(f),MPI_LOGICAL, g,size(f),MPI_LOGICAL, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f,size(f),MPI_LOGICAL, g,size(f),MPI_LOGICAL, comm, i )
#endif
#endif
  END SUBROUTINE l_rsdft_allgather12

  subroutine i_rsdft_allgather12( f, g, comm )
    implicit none
    integer,intent(inout) :: f(:), g(:,0:)
    integer,intent(in) :: comm
    integer :: i
    integer,allocatable :: t(:)
#ifdef _NO_MPI_
    g(:,0) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f)) ) ; t=f
    call MPI_Allgather( t,size(f),MPI_INTEGER,g,size(f),MPI_INTEGER,comm,i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f,size(f),MPI_INTEGER,g,size(f),MPI_INTEGER,comm,i )
#endif
#endif
  end subroutine i_rsdft_allgather12

  subroutine i_rsdft_allgather01( f, g, comm )
    implicit none
    integer,intent(inout) :: f, g(0:)
    integer,intent(in) :: comm
    integer :: i
#ifdef _NO_MPI_
    g(0) = f
    return
#else
    include 'mpif.h'
    call MPI_Allgather( f,1,MPI_INTEGER, g,1,MPI_INTEGER, comm, i )
#endif
  end subroutine i_rsdft_allgather01

  SUBROUTINE d_rsdft_allgather12( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:,0) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather12

  SUBROUTINE d_rsdft_allgather23( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:,0) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather23

  SUBROUTINE d_rsdft_allgather34( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:,:), g(:,:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:,0) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather34

  SUBROUTINE d_rsdft_allgather11( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather11

  SUBROUTINE d_rsdft_allgather22( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather22

  SUBROUTINE z_rsdft_allgather23( f, g, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:), g(:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    complex(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:,0) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER(t,size(f),RSDFT_MPI_COMPLEX16,g,size(f),RSDFT_MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHER(f,size(f),RSDFT_MPI_COMPLEX16,g,size(f),RSDFT_MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgather23

!-------------------------------------------------- rsdft_allgatherv

  SUBROUTINE i_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    integer,intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    integer,allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV( t,size(f),MPI_INTEGER, g,ir,id,MPI_INTEGER, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f,size(f),MPI_INTEGER, g,ir,id,MPI_INTEGER, comm, i )
#endif
#endif
  END SUBROUTINE i_rsdft_allgatherv22

  SUBROUTINE d_rsdft_allgatherv11( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv11

  SUBROUTINE d_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv22

  SUBROUTINE d_rsdft_allgatherv33( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:,:), g(:,:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv33

  SUBROUTINE z_rsdft_allgatherv11( f, g, ir, id, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:), g(:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    complex(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_COMPLEX16, g, ir, id, MPI_COMPLEX16, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_COMPLEX16, g, ir, id, MPI_COMPLEX16, comm, i )
#endif
#endif
  END SUBROUTINE z_rsdft_allgatherv11

  SUBROUTINE z_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    complex(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV(t,size(f),RSDFT_MPI_COMPLEX16,g,ir,id,RSDFT_MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHERV(f,size(f),RSDFT_MPI_COMPLEX16,g,ir,id,RSDFT_MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgatherv22

  SUBROUTINE z_rsdft_allgatherv33( f, g, ir, id, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:,:), g(:,:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    complex(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHERV(t,size(f),RSDFT_MPI_COMPLEX16,g,ir,id,RSDFT_MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHERV(f,size(f),RSDFT_MPI_COMPLEX16,g,ir,id,RSDFT_MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgatherv33


end module rsdft_mpi_module
