MODULE gradient_module

  use density_module, only: density
  use grid_module, only: grid
  use bc_module, only: www, bcset
  use fd_module, only: fd,construct_nabla_fd,destruct_nabla_fd
  use lattice_module, only: lattice,construct_aa_lattice,get_reciprocal_lattice

  implicit none

  PRIVATE
  PUBLIC :: gradient, construct_gradient, destruct_gradient

  type gradient
     real(8),allocatable :: gx(:),gy(:),gz(:)
     real(8),allocatable :: gg(:)
  end type gradient

CONTAINS

  SUBROUTINE construct_gradient( rgrid, rho, grad )
    implicit none
    type(grid),intent(IN)      :: rgrid
    type(density),intent(IN)   :: rho
    type(gradient),intent(OUT) :: grad
    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,s,m,Md,m0,m1
    real(8) :: g1,g2,g3,pi2,b(3,3)

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    call construct_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0d0*acos(-1.0d0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%SizeGrid(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%SizeGrid(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%SizeGrid(3) )

! ---

    m0 = rgrid%SubGrid(1,0)
    m1 = rgrid%SubGrid(2,0)

    allocate( grad%gx(m0:m1) ) ; grad%gx=0.0d0
    allocate( grad%gy(m0:m1) ) ; grad%gy=0.0d0
    allocate( grad%gz(m0:m1) ) ; grad%gz=0.0d0
    allocate( grad%gg(m0:m1) ) ; grad%gg=0.0d0

    www(:,:,:,:)=0.0d0
    do s=1,rho%nn
       i=m0-1
       do i3=rgrid%SubGrid(1,3),rgrid%SubGrid(2,3)
       do i2=rgrid%SubGrid(1,2),rgrid%SubGrid(2,2)
       do i1=rgrid%SubGrid(1,1),rgrid%SubGrid(2,1)
          i=i+1
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho%rho(i,s)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    i=m0-1
    do i3=rgrid%SubGrid(1,3),rgrid%SubGrid(2,3)
    do i2=rgrid%SubGrid(1,2),rgrid%SubGrid(2,2)
    do i1=rgrid%SubGrid(1,1),rgrid%SubGrid(2,1)
       g1=0.0d0
       g2=0.0d0
       g3=0.0d0
       do m=1,Md
          g1 = g1 - nabla%coef(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          g2 = g2 - nabla%coef(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          g3 = g3 - nabla%coef(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
       end do
       i=i+1
       grad%gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
       grad%gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
       grad%gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
    end do
    end do
    end do

    do i=m0,m1
       grad%gg(i) = grad%gx(i)**2 + grad%gy(i)**2 + grad%gz(i)**2
    end do

    call destruct_nabla_fd( nabla )

  END SUBROUTINE construct_gradient


  SUBROUTINE destruct_gradient( grad )
    implicit none
    type(gradient) :: grad
    deallocate( grad%gg )
    deallocate( grad%gz )
    deallocate( grad%gy )
    deallocate( grad%gx )
  END SUBROUTINE destruct_gradient


END MODULE gradient_module
