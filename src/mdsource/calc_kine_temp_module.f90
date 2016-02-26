MODULE calc_kine_temp_module

  use atom_module, only: Natom, ki_atom, zn_atom, md_atom
  use cpmd_variables, only: pmass, AMU, Ndof
  use force_module, only: tim

  implicit none

  PRIVATE
  PUBLIC :: calc_Ndof, get_Ndof
  PUBLIC :: calc_kine
  PUBLIC :: calc_temp

CONTAINS


  SUBROUTINE calc_Ndof   ! Degree of freedom
    implicit none
    real(8) :: v(3),u(3)
    integer :: i,m
    call random_number( v )
    Ndof=0
    do i=1,Natom
       m=md_atom(i)
       u(:) = matmul( tim(:,:,m), v(:) )
       Ndof = Ndof + count( u /= 0.0d0 )
    end do
    Ndof = Ndof - 3  ! Subtract center-of-mass dof
  END SUBROUTINE calc_Ndof


  SUBROUTINE get_Ndof( Ndof_out )
    implicit none
    integer,intent(OUT) :: Ndof_out
    if ( Ndof == 0 ) call calc_Ndof
    Ndof_out = Ndof
  END SUBROUTINE get_Ndof


  SUBROUTINE calc_kine( Velocity, kine )
    implicit none
    real(8),intent(IN)  :: Velocity(:,:)
    real(8),intent(OUT) :: kine
    integer :: i
    real(8) :: pm
    kine=0.0d0
    do i=1,Natom
       pm = pmass( zn_atom(ki_atom(i)) )
       kine = kine + sum( Velocity(:,i)**2 )*pm
    end do
    kine=kine*0.5d0*AMU
  END SUBROUTINE calc_kine

  SUBROUTINE calc_temp( Velocity, temp )
    implicit none
    real(8),intent(IN)  :: Velocity(:,:)
    real(8),intent(OUT) :: temp
    real(8),parameter :: TOJOUL=4.35975d-18 ! (J/hartree)
    real(8),parameter :: kB_J=1.380658d-23  ! (J/K)
    real(8),parameter :: kB=kB_J/TOJOUL     ! (hartree/K)
    real(8) :: kine
    call calc_kine( Velocity, kine )
    temp = kine/( 0.5d0*Ndof*kB )
  END SUBROUTINE calc_temp

END MODULE calc_kine_temp_module
