module kinetic_allgatherv_module

  use kinetic_allgatherv_sub_module
  use watch_module, only: watchb, time_kine

  implicit none

  private
  public :: init_kinetic_allgatherv
  public :: op_kinetic_allgatherv

#ifdef _DRSDFT_
  real(8),allocatable :: tpsi(:)
  real(8),parameter :: zero=0.0d0
  real(8),allocatable :: kin_coef(:,:,:)
#else
  complex(8),allocatable :: tpsi(:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),allocatable :: kin_coef(:,:,:)
#endif

  integer :: comm
  integer :: Igrid(2,0:3)
  integer,allocatable :: ir_grid(:), id_grid(:)
  integer,allocatable :: list_vec(:,:)
  logical :: init_flag=.false.

contains


  subroutine init_kinetic_allgatherv( Igrid_in, comm_in )

    implicit none
    integer,intent(in) :: Igrid_in(1:2,0:3),comm_in
    integer :: n, i, ierr
    include 'mpif.h'

    call write_border( 0, "init_kinetic_allgatherv(start)" )

    comm = comm_in 
    call MPI_Comm_size( comm, n, ierr )
    allocate( ir_grid(0:n-1) ); ir_grid=0
    allocate( id_grid(0:n-1) ); id_grid=0

    Igrid = Igrid_in
    n=Igrid(2,0)-Igrid(1,0)+1
    call MPI_Allgather( n,1,MPI_INTEGER, ir_grid,1,MPI_INTEGER, comm,ierr )

    do i=0,size(ir_grid)-1
       id_grid(i) = sum( ir_grid(0:i) ) - ir_grid(i)
    end do

    n=sum(ir_grid)
    allocate( tpsi(n) ); tpsi=zero

    call kinetic_allgatherv_sub( Igrid, comm, list_vec, kin_coef )

    init_flag=.true.

    call write_border( 0, "init_kinetic_allgatherv(end)" )

  end subroutine init_kinetic_allgatherv


  subroutine op_kinetic_allgatherv( psi, Kpsi, k )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(in) ::  psi(:,:)
    real(8),intent(inout) :: Kpsi(:,:)
    real(8) :: tmp
#else
    complex(8),intent(in) ::  psi(:,:)
    complex(8),intent(inout) :: Kpsi(:,:)
    complex(8) :: tmp
#endif
    integer,intent(in) :: k
    integer :: m, n, nb, ib, i, j, ierr
    real(8) :: ttmp(2)
    include 'mpif.h'

    if ( .not.init_flag ) call stop_program("stop@op_kinetic_allgatherv(1)")

    m =size(list_vec,1)
    n =size(psi,1)
    nb=size(psi,2)

    do ib=1,nb

       call watchb( ttmp, barrier='on' )

#ifdef _DRSDFT_
       call MPI_Allgatherv( psi(:,ib), n, MPI_REAL8, tpsi, &
                            ir_grid, id_grid, MPI_REAL8, comm, ierr )
#else
       call MPI_Allgatherv( psi(:,ib), n, MPI_COMPLEX16, tpsi, &
                            ir_grid, id_grid, MPI_COMPLEX16,comm,ierr )
#endif

       call watchb( ttmp, time_kine(:,1), barrier='on' )

       do i=1,n

          tmp=zero
          do j=1,m
             tmp = tmp + kin_coef(j,i,k)*tpsi( list_vec(j,i) )
          end do ! j

          Kpsi(i,ib) = Kpsi(i,ib) + tmp

       end do ! i

       call watchb( ttmp, time_kine(:,2), barrier='on' )

    end do ! ib

  end subroutine op_kinetic_allgatherv


end module kinetic_allgatherv_module
