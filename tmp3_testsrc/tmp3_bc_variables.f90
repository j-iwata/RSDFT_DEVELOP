MODULE bc_variables

  implicit none

  PRIVATE
  PUBLIC :: Md,n_neighbor,fdinfo_send,fdinfo_recv,www,sbuf,rbuf,zero &
           ,TYPE_MAIN,init_bc_test,allocate_bc_test

  integer :: Md
  integer :: n_neighbor(6)
  integer,allocatable :: fdinfo_send(:,:,:)
  integer,allocatable :: fdinfo_recv(:,:,:)
  include 'mpif.h'

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: www(:,:,:,:)
  real(8),allocatable :: sbuf(:,:,:)
  real(8),allocatable :: rbuf(:,:,:)
  real(8),parameter :: zero=0.d0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: www(:,:,:,:)
  complex(8),allocatable :: sbuf(:,:,:)
  complex(8),allocatable :: rbuf(:,:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
#endif

CONTAINS


  SUBROUTINE init_bc_test(n,fdinfo_send_in,fdinfo_recv_in,n_neighbor_in)
    implicit none
    integer,intent(IN) :: n
    integer,intent(IN) :: fdinfo_send_in(10,n,6),fdinfo_recv_in(10,n,6)
    integer,intent(IN) :: n_neighbor_in(6)
    n_neighbor(:)=n_neighbor_in(:)
    allocate( fdinfo_send(10,n,6) ) ; fdinfo_send=0
    allocate( fdinfo_recv(10,n,6) ) ; fdinfo_recv=0
    fdinfo_send(:,:,:)=fdinfo_send_in(:,:,:)
    fdinfo_recv(:,:,:)=fdinfo_recv_in(:,:,:)
  END SUBROUTINE init_bc_test


  SUBROUTINE allocate_bc_test(MB_d,Igrid)
    implicit none
    integer,intent(IN) :: MB_d,Igrid(2,0:3)
    integer :: a1,a2,a3,b1,b2,b3
    a1=Igrid(1,1)-Md ; a2=Igrid(1,2)-Md ; a3=Igrid(1,3)-Md
    b1=Igrid(2,1)+Md ; b2=Igrid(2,2)+Md ; b3=Igrid(2,3)+Md
    if ( allocated(www) ) deallocate(www)
    allocate( www(a1:b1,a2:b2,a3:b3,MB_d) )
    www(:,:,:,:)=zero
    a1=maxval( n_neighbor(1:6) )
    a2=maxval( fdinfo_send(9,1:a1,1:6) )*MB_d
    a3=maxval( fdinfo_recv(9,1:a1,1:6) )*MB_d
    if ( allocated(sbuf) ) deallocate(sbuf)
    if ( allocated(rbuf) ) deallocate(rbuf)
    allocate( sbuf(a2,a1,6) )
    allocate( rbuf(a3,a1,6) )
    sbuf(:,:,:)=zero
    rbuf(:,:,:)=zero
  END SUBROUTINE allocate_bc_test

END MODULE bc_variables
