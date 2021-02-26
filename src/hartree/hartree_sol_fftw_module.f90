module hartree_sol_fftw_module

  use hartree_variables, only: E_hartree, Vh
  use bb_module, only: bb
  use rgrid_module, only: Ngrid, Igrid, dV
  use ggrid_module, only: NGgrid, LLG, construct_ggrid, destruct_ggrid
  use grid_module, only: inner_product_grid
  use fftw_module, only: ML1_c, ML2_c, ML3_c0, N_ML3_c &
                       , zwork3_ptr0, zwork3_ptr1 &
                       , d1_to_z3_fftw, z3_to_d1_fftw &
                       , plan_forward, plan_backward
  use,intrinsic :: iso_c_binding

  implicit none

  private
  public :: init_hartree_sol_fftw
  public :: calc_hartree_sol_fftw

  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  integer,allocatable :: IGHT(:,:)
  real(8),allocatable :: GGHT(:)

  logical :: flag_init_done = .false.

contains

  subroutine init_hartree_sol_fftw( Ngrid )
    implicit none
    integer,intent(in) :: Ngrid(3)
    integer :: n,i,i1,i2,i3,NLG,ML
    real(8) :: const, g2

    if ( flag_init_done ) return

    call write_border( 0, ' init_hartree_fftw(start)' )

    call construct_Ggrid(0)

    ML = product( Ngrid(1:3) )
    NLG = size( LLG, 2 )

    n=0
    do i=1,NLG
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( all(LLG(1:3,i)==0) ) cycle
      if ( ML3_c0 <= i3 .and. i3 <= ML3_c0+N_ML3_c-1 ) then
        n=n+1
      end if
    end do
    NGHT=n

    allocate( LGHT(3,NGHT) ); LGHT=0
    allocate( IGHT(3,NGHT) ); IGHT=0
    allocate( GGHT(NGHT) ); GGHT=0.0d0

    const = 4.0d0*acos(-1.0d0)/ML

    n=0
    do i=1,NLG
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( all(LLG(1:3,i)==0) ) cycle
      if ( ML3_c0 <= i3 .and. i3 <= ML3_c0+N_ML3_c-1 ) then
        n=n+1
        LGHT(1,n)=i1+1
        LGHT(2,n)=i2+1
        LGHT(3,n)=i3+1-ML3_c0
        g2=( bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i) )**2 &
          +( bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i) )**2 &
          +( bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i) )**2
        GGHT(n)=const/g2
      end if
    end do

    call destruct_Ggrid

    flag_init_done = .true.

    call write_border( 0, ' init_hartree_fftw(end)' )

  end subroutine init_hartree_sol_fftw

  subroutine calc_hartree_sol_fftw( rho )
    implicit none
    real(8),intent(in) :: rho(:,:)
#ifdef _FFTW_
    integer :: i,i1,i2,i3,j1,j2,j3,n,ng,ns
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8),allocatable :: work(:)
    include 'fftw3-mpi.f03'

    call write_border( 1, " calc_hartree_sol_fftw(start)" )
    if ( .not.flag_init_done ) call stop_program('Call "init_hartree_fftw" first.')

    ng = size( rho, 1 )
    ns = size( rho, 2 )

    allocate( work(ng) )
    work(:)=rho(:,1)
    do n=2,ns
      work(:)=work(:)+rho(:,n)
    end do
    call d1_to_z3_fftw( work, zwork3_ptr0 )

    call fftw_mpi_execute_dft( plan_forward, zwork3_ptr0, zwork3_ptr1 )

    zwork3_ptr0(:,:,:)=z0
    do i=1,NGHT
      i1=LGHT(1,i)
      i2=LGHT(2,i)
      i3=LGHT(3,i)
      zwork3_ptr0(i1,i2,i3) = zwork3_ptr1(i1,i2,i3)*GGHT(i)
    end do

    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )

    call z3_to_d1_fftw( zwork3_ptr1, Vh )

    call inner_product_grid( work, Vh, 0.5d0*dV, E_hartree ) 

    deallocate( work )

    call write_border( 1, " calc_hartree_sol_fftw(end)" )
#endif
  end subroutine calc_hartree_sol_fftw


end module hartree_sol_fftw_module
