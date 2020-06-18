MODULE force_module

  use parallel_module, only: disp_switch_parallel
  use atom_module, only: aa_atom, ki_atom, md_atom, Natom
  use symmetry_module, only: sym_force
  use force_sol_module
  use force_mol_module
  use io_tools_module, only: IOTools_findKeyword, IOTools_readIntegerKeyword
  use hpc_module
!  use force_correction_module

  implicit none

  PRIVATE
  PUBLIC :: calc_force, get_fmax_force
  PUBLIC :: init_force
  PUBLIC :: tim

  integer :: Ntim
  real(8),allocatable :: tim(:,:,:)

  integer :: SYStype   = 0
  logical :: flag_init = .true.

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_force( SYStype_in )
    implicit none
    integer,intent(IN) :: SYStype_in
    if ( .not.flag_init ) return
    flag_init = .false.
    call write_border( 0, " init_force(start)" )
    Ntim=maxval( md_atom )
    if ( Ntim <= 0 ) then
       Ntim=1
       md_atom(:)=1
    end if
    if ( disp_switch_parallel ) write(*,*) "Ntim=",Ntim
    allocate( tim(3,3,Ntim) ) ; tim=0.0d0
    tim(1,1,1:Ntim) = 1.0d0
    tim(2,2,1:Ntim) = 1.0d0
    tim(3,3,1:Ntim) = 1.0d0
    call read_force
    SYStype = SYStype_in
    call write_border( 0, " init_force(end)" )
  END SUBROUTINE init_force


  SUBROUTINE read_force
    implicit none
    integer :: i,j,k,mchk,unit
    logical :: exist_keyword
    real(8) :: dummy(2)
    character(8) :: cbuf
    call write_border( 0, " read_force(start)" )
!    call IOTools_readIntegerKeyword( "SYSTYPE", SYStype )
    call IOTools_findKeyword( "TIM", exist_keyword, unit_out=unit )
    if ( exist_keyword ) then
       mchk=0
       do k=1,Ntim
          do i=1,3
             read(unit,*) ( tim(i,j,k), j=1,3 )
             mchk=mchk+1
          end do ! j
       end do ! k
       if ( mchk /= 3*Ntim ) then
          write(*,*) "WARNING: Insufficient data for tim-matrix !!!"
          write(*,*) "mchk, 3*Ntim =",mchk,3*Ntim
       end if
       do k=1,Ntim
          write(*,'(1x,"tim(",i1,")",2x,3f10.5)') k,(tim(1,j,k),j=1,3)
          do i=2,3
             write(*,'(1x,6x,2x,3f10.5)') (tim(i,j,k),j=1,3)
          end do
       end do
    end if
    call mpi_bcast( tim, 9*Ntim, MPI_REAL8, 0, MPI_COMM_WORLD, i )
    call write_border( 0, " read_force(end)" )
  END SUBROUTINE read_force


  SUBROUTINE calc_force( MI, force, fmax, disp, unit )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    real(8),optional,intent(OUT) :: fmax
    logical,optional,intent(IN) :: disp
    integer,optional,intent(IN) :: unit

!    if ( flag_init ) call init_force

    select case( SYStype )
    case( 0 )
       call calc_force_sol( MI, force )
    case( 1 )
       call calc_force_mol( MI, force )
    end select

! --- correction term ---

    !call force_correction( force )

! --- constraint & symmetry ---

    call sym_force ( MI, ki_atom, aa_atom, force )

    if ( present(disp) ) call write_info_force( force, ki_atom, disp )
    if ( present(unit) ) then
       open( unit, file='forcelog_bare' )
       write(unit,'("# XYZ-components of atomic force (unit is in Hartree/bohr)")')
       write(unit,'("# Bare force")')
       call write_force_to_file( unit, force, ki_atom )
       close( unit )
    end if

    if ( Ntim > 0 ) then
       call constraint_force( MI, force )
       if ( present(disp) ) then
          if ( disp ) write(*,*) "(with constraint)"
          call write_info_force( force, ki_atom, disp )
       end if
       if ( present(unit) ) then
          open( unit, file='forcelog_constraint' )
          write(unit,'("# XYZ-components of atomic force (unit is in Hartree/bohr)")')
          write(unit,'("# Constraint force")')
          call write_force_to_file( unit, force, ki_atom )
          close( unit )
       end if
    end if

    call force_projection_hpc( MI, force )

    if ( present(fmax) ) call get_fmax_force( fmax, force )

  END SUBROUTINE calc_force


  SUBROUTINE constraint_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(INOUT) :: force(3,MI)
    integer :: i,j
    do i=1,MI
       j=md_atom(i)
       if ( j > 0 ) then
          force(:,i) = matmul( tim(:,:,j), force(:,i) )
       end if
    end do
  END SUBROUTINE constraint_force


  SUBROUTINE get_fmax_force( fmax, force_in )
    implicit none
    real(8),intent(OUT) :: fmax
    real(8),optional,intent(IN) :: force_in(:,:)
    real(8),allocatable :: force(:,:)
    if ( present(force_in) ) then
       call get_fmax_force_sub( fmax, force_in )
    else
       allocate( force(3,Natom) ) ; force=0.0d0
       call calc_force( Natom, force )
       call get_fmax_force_sub( fmax, force )
       deallocate( force )
    end if
  END SUBROUTINE get_fmax_force

  SUBROUTINE get_fmax_force_sub( fmax, force )
    implicit none
    real(8),intent(OUT) :: fmax
    real(8),intent(IN) :: force(:,:)
    real(8) :: ff
    integer :: a
    fmax=-1.0d10
    do a=1,Natom
       ff=force(1,a)**2+force(2,a)**2+force(3,a)**2
       fmax=max(fmax,ff)
    end do
    fmax=sqrt(fmax)
  END SUBROUTINE get_fmax_force_sub


  SUBROUTINE write_info_force( f, k, disp )
    implicit none
    real(8),intent(IN) :: f(:,:)
    integer,intent(IN) :: k(:)
    logical,intent(IN) :: disp
    integer :: N,a
    N=size( f, 2 )
    if ( disp ) then
       if (  N <= 11 ) then
          do a=1,N
             write(*,'(1x,i4,i3,3g21.12)') a,k(a),f(1:3,a)
          end do
       else
          do a=1,min(5,N)
             write(*,'(1x,i4,i3,3g21.12)') a,k(a),f(1:3,a)
          end do
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          do a=N-5,N
             write(*,'(1x,i4,i3,3g21.12)') a,k(a),f(1:3,a)
          end do
       end if
    end if
  END SUBROUTINE write_info_force


  SUBROUTINE write_force_to_file( unit, f, k )
    implicit none
    integer,intent(IN) :: unit
    real(8),intent(IN) :: f(:,:)
    integer,intent(IN) :: k(:)
    integer :: a
    do a=1,size( f, 2 )
       write(unit,'(1x,i3,3g21.12)') k(a),f(1:3,a)
    end do
  END SUBROUTINE write_force_to_file


END MODULE force_module
