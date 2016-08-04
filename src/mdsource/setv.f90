!-----------------------------------------------------------------------
!     Set initial velocity
!-----------------------------------------------------------------------
SUBROUTINE setv( temp, Velocity )

  use atom_module, only: Natom, md_atom
  use force_module, only: tim
  use calc_kine_temp_module, only: calc_temp

  implicit none
  real(8),intent(IN)  :: temp
  real(8),intent(OUT) :: Velocity(3,Natom)
  integer :: i,m
  real(8) :: temp_now,rescale,pi2,vtmp(3)
  real(8),allocatable :: r(:)

  pi2=2.d0*acos(-1.d0)
  allocate( r(6*Natom) )
  call rnum(6*Natom,r)
  do i=1,Natom
     vtmp(1)=r(i*6-5)*sin(pi2*r(6*i-4))
     vtmp(2)=r(i*6-3)*sin(pi2*r(6*i-2))
     vtmp(3)=r(i*6-1)*sin(pi2*r(6*i  ))
     m=md_atom(i)
     Velocity(:,i) = matmul( tim(:,:,m), vtmp(:) )
  end do
  deallocate(r)
  call calc_temp( Velocity, temp_now )
  rescale=sqrt(temp/temp_now)
  Velocity(:,:)=Velocity(:,:)*rescale
  call vcom( Velocity ) ! Center-of-mass velocity
  return

CONTAINS

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

END SUBROUTINE setv
