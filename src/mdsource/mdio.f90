!-----------------------------------------------------------------------
!     IO routine for molecular dynamics
!-----------------------------------------------------------------------
subroutine mdio( io, tote0 )
  use cpmd_variables, only: Rion, Velocity, Force
  implicit none
  integer,intent(in) :: io
  real(8),intent(inout) :: tote0
  integer :: i,k,natom
  real(8) :: dum(3)
  character(50) :: aline
  call write_border( 0, 'mdio(start)' )
  natom=size(Rion,2)
  select case( io )
  case( 0 )
    open( 2, file='initial.dat', status='unknown' )
    read(2,'(a)') aline
    do i=1,natom
      read(2,'(3f24.16)') (Rion(k,i),k=1,3)
    end do
    read(2,'(a)') aline
    do i=1,natom
      read(2,'(3f24.16)') (Velocity(k,i),k=1,3)
    end do
    read(2,'(a)') aline
    do i=1,natom
      read(2,'(3f24.16)') (dum(k),k=1,3)
    end do
    read(2,*) tote0
    close( 2 )
  case( 1 )
    open( 2, file='final.dat', status='replace' )
    write(2,'(a)') "final coordinate"
    do i=1,natom
      write(2,'(3g24.16)') (Rion(k,i),k=1,3)
    end do
    write(2,'(a)') "final velocity"
    do i=1,natom
      write(2,'(3g24.16)') (Velocity(k,i),k=1,3)
    end do
    write(2,'(a)') "final force"
    do i=1,natom
      write(2,'(3g24.16)') (Force(k,i),k=1,3)
    end do
    write(2,*) tote0
    close( 2 )
  case( 2 )
    write(*,'(a)') "initial coordinate"
    do i=1,natom
      write(*,'(3f15.8)') (Rion(k,i),k=1,3)
    end do
    write(*,'(a)') "initial velocity"
    do i=1,natom
      write(*,'(3f15.8)') (Velocity(k,i),k=1,3)
    end do
    write(*,'(a)') "initial force"
    do i=1,natom
      write(*,'(3f15.8)') (Force(k,i),k=1,3)
    end do
    write(*,*) tote0
   end select
   call write_border( 0, 'mdio(end)' )
   return
end subroutine mdio
