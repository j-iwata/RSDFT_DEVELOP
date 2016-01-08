!-----------------------------------------------------------------------
!     Set initial velocity
!-----------------------------------------------------------------------
subroutine setv
  use atom_module, only: Natom
  use cpmd_variables, only: Velocity,temp
  implicit none
  integer i
  real(8) kine,rescale,pi2
  real(8),allocatable :: r(:)
  pi2=2.d0*acos(-1.d0)
  allocate( r(6*Natom) )
  call rnum(6*Natom,r)
  do i=1,Natom
     Velocity(1,i)=r(i*6-5)*sin(pi2*r(6*i-4))
     Velocity(2,i)=r(i*6-3)*sin(pi2*r(6*i-2))
     Velocity(3,i)=r(i*6-1)*sin(pi2*r(6*i  ))
  enddo
  deallocate(r)
  call calkin(kine)
  rescale=sqrt(1.5d0*Natom*temp/kine)
  Velocity=Velocity*rescale
  return
contains 
  subroutine rnum(N,randv)
    implicit none
    integer,intent(in)    :: N
    real(8),intent(inout) :: randv(N)
    real(8),parameter :: a=1953125.0d0
    real(8),parameter :: b=33554432.0d0
    real(8),parameter :: r0=177147.0d0
    real(8) :: r
    integer :: i
    r=r0
    do i=1,N
       r=mod(a*r,b)
       randv(i)=r/b
    enddo
    return
  end subroutine rnum
end subroutine setv
