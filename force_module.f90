MODULE force_module

  use parallel_module, only: disp_switch_parallel

!--kuchida_2015_0106
!  use atom_module, only: md_atom
  use atom_module, only: aa_atom, ki_atom, md_atom
  use symmetry_module  

  use force_sol_module
  use force_mol_module

  implicit none

  PRIVATE
  PUBLIC :: init_force, calc_force

  integer :: Ntim
  real(8),allocatable :: tim(:,:,:)

  integer :: SYStype
  logical :: flag_init = .true.

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_force( rank, unit, SYStype_in )
    implicit none
    integer,intent(IN) :: rank, unit, SYStype_in
    if ( disp_switch_parallel ) &
         write(*,'(1x,a60," init_force")') repeat("-",60)
    if ( flag_init ) flag_init=.false.
    SYStype = SYStype_in
    Ntim=0
    Ntim=maxval( md_atom )
    if ( rank == 0 ) then
       write(*,*) "SYStype=",SYStype
       write(*,*) "Ntim=",Ntim
    end if
    if ( Ntim <= 0 ) return
    allocate( tim(3,3,Ntim) ) ; tim=0.0d0
    tim(1,1,1) = 1.0d0
    tim(2,2,1) = 1.0d0
    tim(3,3,1) = 1.0d0
    call read_force( rank, unit )
  END SUBROUTINE init_force


  SUBROUTINE read_force( rank, unit )
    implicit none
    integer,intent(IN) :: rank, unit
    integer :: i,j,k,kr,ierr,mchk
    character(3) :: cbuf,ckey
    mchk=0
    if ( rank == 0 ) then
       write(*,'(1x,a40," read_force")') repeat("-",40)
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital( cbuf, ckey )
          if ( ckey == "TIM" ) then
             do j=1,Ntim
                do k=1,3
                   read(unit,*) ( tim(kr,k,j), kr=1,3 )
                   mchk=mchk+1
                end do ! k
             end do ! j
             exit
          end if
       end do ! i
999    continue
       if ( mchk /= 3*Ntim ) then
          write(*,*) "WARNING: Insufficient data for tim-matrix !!!"
          write(*,*) "mchk, 3*Ntim =",mchk,3*Ntim
       end if
       do j=1,Ntim
          write(*,'(1x,"tim(",i1,")",2x,3f10.5)') j,(tim(kr,1,j),kr=1,3)
          do k=2,3
             write(*,'(1x,6x,2x,3f10.5)') (tim(kr,k,j),kr=1,3)
          end do
       end do
    end if
    call mpi_bcast( tim, 9*Ntim, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
  END SUBROUTINE read_force


  SUBROUTINE calc_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)

    if ( flag_init ) then
       write(*,*) "init_force must be called first"
       stop "stop@calc_force"
    end if

    select case( SYStype )
    case( 0 )
       call calc_force_sol( MI, force )
    case( 1 )
       call calc_force_mol( MI, force )
    end select

! --- constraint & symmetry ---

!    call symforce
!--kuchida_2015_0106
    call sym_force ( MI, ki_atom, aa_atom, force )

    if ( Ntim > 0 ) call constraint_force( MI, force )

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


END MODULE force_module
