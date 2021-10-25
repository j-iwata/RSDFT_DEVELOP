module calc_kine_temp_module

  implicit none

  private
  public :: calc_Ndof, get_Ndof
  public :: calc_kine
  public :: calc_temp

  real(8),parameter :: TOJOUL=4.35975d-18 ! (J/hartree)
  real(8),parameter :: kB_J=1.380658d-23  ! (J/K)
  real(8),parameter :: kB=kB_J/TOJOUL     ! (hartree/K)

  integer :: N_degree_of_freedom=0

contains


  subroutine calc_Ndof ! Degree of freedom
    use atom_module, only: md_atom
    use force_module, only: tim, init_force
    implicit none
    real(8) :: v(3),u(3)
    integer :: i,m,natom,Ndof
    call init_force
    call random_number( v )
    natom=size(md_atom)
    Ndof=0
    do i=1,natom
      m=md_atom(i)
      u(:) = matmul( tim(:,:,m), v(:) )
      Ndof = Ndof + count( u /= 0.0d0 )
    end do
    Ndof = Ndof - 3  ! Subtract center-of-mass dof
    if ( Ndof <= 0 ) Ndof=Ndof+3
    N_degree_of_freedom=Ndof
  end subroutine calc_Ndof


  subroutine get_Ndof( Ndof_out )
    implicit none
    integer,intent(out) :: Ndof_out
    if ( N_degree_of_freedom == 0 ) call calc_Ndof
    Ndof_out = N_degree_of_freedom
  end subroutine get_Ndof


  subroutine calc_kine( Velocity, kine, temp )
    use atom_module, only: zn_atom, ki_atom
    use cpmd_variables, only: pmass, AMU
    implicit none
    real(8),intent(in)  :: Velocity(:,:)
    real(8),intent(out) :: kine
    real(8),optional,intent(out) :: temp
    integer :: i,natom,Ndof
    real(8) :: pm
    natom=size(Velocity,2)
    kine=0.0d0
    do i=1,natom
      pm = pmass( zn_atom(ki_atom(i)) )
      kine = kine + sum( Velocity(:,i)**2 )*pm
    end do
    kine=kine*0.5d0*AMU
    if ( present(temp) ) then
      call get_Ndof( Ndof )
      temp=kine/( 0.5d0*Ndof*kB )
    end if
  end subroutine calc_kine


  subroutine calc_temp( Velocity, temp )
    implicit none
    real(8),intent(in)  :: Velocity(:,:)
    real(8),intent(out) :: temp
    real(8) :: kine
    integer :: Ndof
    call calc_kine( Velocity, kine )
    call get_Ndof( Ndof )
    temp = kine/( 0.5d0*Ndof*kB )
  end subroutine calc_temp

end module calc_kine_temp_module
