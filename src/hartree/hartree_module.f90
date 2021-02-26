module hartree_module

  use hartree_variables
  use hartree_sol_module
  use hartree_sol_1dffte_module
  use hartree_sol_ffte_module
  use hartree_mol_module
  use hartree_ene_module

  implicit none

  private
  public :: init_hartree
  public :: calc_hartree

  integer :: SYStype=0
  integer :: Md
  integer :: n1, n2

  character(5),parameter :: iswitch_fft = 'FFTE1'

contains

  subroutine init_hartree( Igrid, Ngrid, Md_in, SYStype_in )
    implicit none
    integer,intent(in) :: Igrid(2,0:3), Ngrid(0:3), Md_in, SYStype_in

    call write_border( 0, " init_hartree(start)" )

    n1      = Igrid(1,0)
    n2      = Igrid(2,0)
    Md      = Md_in
    SYStype = SYStype_in

    allocate( Vh(n1:n2) ); Vh=0.0d0

    select case( SYStype )
    case default

      select case( iswitch_fft )
      case( 'FFTE1' )
        call init_hartree_sol_ffte( Ngrid(1:3), Igrid(:,1:3) )
      case( 'FFTE2' )
        call init_hartree_sol_1dffte( Ngrid(1:3), Igrid(:,1:3) )
      end select

    case(1)

      call init_hartree_mol(Md)

    end select

    call write_border( 0, " init_hartree(end)" )

  end subroutine init_hartree


  subroutine calc_hartree( rho, vout )
    implicit none
    real(8),intent(in) :: rho(:,:)
    real(8),optional,intent(inout) :: vout(:)
    real(8),allocatable :: trho(:),tvh(:)
    integer :: i,s,ng,ns

    call write_border( 1, " calc_hartree(start)" )

    ng = size( rho, 1 )
    ns = size( rho, 2 )

    if ( present(vout) ) then
      allocate( tvh(ng) )
      tvh=Vh
    end if

    select case(SYStype)
    case default

      call calc_hartree_sol( rho, iswitch_fft )

    case(1)

      allocate( trho(ng) )

!$OMP parallel do
      do i=1,ng
        trho(i) = rho(i,1)
      end do
!$OMP end parallel do
      do s=2,ns
!$OMP parallel do
        do i=1,ng
          trho(i) = trho(i) + rho(i,s)
        end do
!$OMP end parallel do
      end do

      call calc_hartree_mol( n1, n2, 1, trho, Vh, E_hartree )
      call calc_hartree_ene( trho, Vh, E_hartree )

      deallocate( trho )

    end select

    if ( present(vout) ) then
      vout=Vh
      Vh=tvh
      deallocate( tvh )
    end if

    call write_border( 1, " calc_hartree(end)" )

  end subroutine calc_hartree

end module hartree_module
