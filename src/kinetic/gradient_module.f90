module gradient_module

  use grid_module, only: grid
  use bc_module, only: www, bcset
  use fd_module, only: fd,construct_nabla_fd,destruct_nabla_fd
  use lattice_module, only: lattice,get_aa_lattice,get_reciprocal_lattice
  use basic_type_factory, only: GSArray
  use nabla_sym_module, only: isymmetry_nabla, init_nabla_sym, op_nabla_sym

  implicit none

  private
  public :: gradient, construct_gradient, destruct_gradient
  public :: gradient16, construct_gradient16, destruct_gradient16
  public :: calc_abc_gradient
  public :: calc_xyz_gradient

  integer,parameter :: DP=kind(0.0d0)
!#ifdef _NO_QPRECISION_
  integer,parameter :: QP=kind(0.0d0)
!#else
!  integer,parameter :: QP=kind(0.0q0)
!#endif

  type gradient
    real(DP),allocatable :: gx(:),gy(:),gz(:)
    real(DP),allocatable :: gg(:)
  end type gradient

  type gradient16
    real(QP),allocatable :: gx(:),gy(:),gz(:)
    real(QP),allocatable :: gg(:)
  end type gradient16

contains


  subroutine construct_gradient( rgrid, rho, grad, s_in )

    implicit none
    type( grid ),intent(in) :: rgrid
    type( GSArray ),intent(in) :: rho
    type( gradient ),intent(out) :: grad
    integer,optional,intent(in) :: s_in
    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,s,m,Md,m0,m1
    real(DP) :: g1,g2,g3,pi2,b(3,3)

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0d0*acos(-1.0d0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

! ---

    m0 = rgrid%g1%head
    m1 = rgrid%g1%tail

    allocate( grad%gx(m0:m1) ) ; grad%gx=0.0d0
    allocate( grad%gy(m0:m1) ) ; grad%gy=0.0d0
    allocate( grad%gz(m0:m1) ) ; grad%gz=0.0d0
    allocate( grad%gg(m0:m1) ) ; grad%gg=0.0d0

    www(:,:,:,:)=0.0d0
    do s=rho%s_range%head_global,rho%s_range%tail_global
      if ( present(s_in) .and. s/=s_in ) cycle
      i=m0-1
      do i3=rgrid%g3%z%head,rgrid%g3%z%tail
      do i2=rgrid%g3%y%head,rgrid%g3%y%tail
      do i1=rgrid%g3%x%head,rgrid%g3%x%tail
        i=i+1
        www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho%val(i,s)
      end do
      end do
      end do
    end do

    call bcset(1,1,Md,0)

    i=m0-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
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

  end subroutine construct_gradient


  subroutine destruct_gradient( grad )
    implicit none
    type(gradient) :: grad
    deallocate( grad%gg )
    deallocate( grad%gz )
    deallocate( grad%gy )
    deallocate( grad%gx )
  end subroutine destruct_gradient


  SUBROUTINE construct_gradient16( rgrid, rho, grad )

    implicit none
    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN)   :: rho
    type( gradient16 ),intent(OUT) :: grad
    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,s,m,Md,m0,m1
    integer :: a1,a2,a3,b1,b2,b3
    real(QP) :: g1,g2,g3,pi2,b(3,3)
    real(QP),allocatable :: w(:,:,:)
    real(8),allocatable :: work(:),gtmp(:,:)

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0_QP*acos(-1.0_QP)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

! ---

    m0 = rgrid%g1%head
    m1 = rgrid%g1%tail
    allocate( grad%gx(m0:m1) ) ; grad%gx=0.0_QP
    allocate( grad%gy(m0:m1) ) ; grad%gy=0.0_QP
    allocate( grad%gz(m0:m1) ) ; grad%gz=0.0_QP
    allocate( grad%gg(m0:m1) ) ; grad%gg=0.0_QP

    call init_nabla_sym

    if ( isymmetry_nabla /= 0 ) then

       allocate( gtmp(m0:m1,3) ) ; gtmp=0.0d0
       allocate( work(m0:m1)   ) ; work=0.0d0

       do s=rho%s_range%head_global,rho%s_range%tail_global
          work(:) = work(:) + rho%val(:,s)
       end do

       call op_nabla_sym( 1, work, gtmp(:,1) )
       call op_nabla_sym( 2, work, gtmp(:,2) )
       call op_nabla_sym( 3, work, gtmp(:,3) )

       grad%gx(:) = b(1,1)*gtmp(:,1) + b(1,2)*gtmp(:,2) + b(1,3)*gtmp(:,3)
       grad%gy(:) = b(2,1)*gtmp(:,1) + b(2,2)*gtmp(:,2) + b(2,3)*gtmp(:,3)
       grad%gz(:) = b(3,1)*gtmp(:,1) + b(3,2)*gtmp(:,2) + b(3,3)*gtmp(:,3)

       deallocate( work )
       deallocate( gtmp )

    else

       www(:,:,:,:)=0.0d0
       do s=rho%s_range%head_global,rho%s_range%tail_global
          i=m0-1
          do i3=rgrid%g3%z%head,rgrid%g3%z%tail
          do i2=rgrid%g3%y%head,rgrid%g3%y%tail
          do i1=rgrid%g3%x%head,rgrid%g3%x%tail
             i=i+1
             www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho%val(i,s)
          end do
          end do
          end do
       end do

       call construct_nabla_fd( nabla )
       Md = nabla%md

       call bcset(1,1,Md,0)

       a1 = rgrid%g3%x%head - Md
       b1 = rgrid%g3%x%tail + Md
       a2 = rgrid%g3%y%head - Md
       b2 = rgrid%g3%y%tail + Md
       a3 = rgrid%g3%z%head - Md
       b3 = rgrid%g3%z%tail + Md
       allocate( w(a1:b1,a2:b2,a3:b3) ) ; w=0.0_QP

       w=www(:,:,:,1)

       i=m0-1
       do i3=rgrid%g3%z%head,rgrid%g3%z%tail
       do i2=rgrid%g3%y%head,rgrid%g3%y%tail
       do i1=rgrid%g3%x%head,rgrid%g3%x%tail
          g1=0.0_QP
          g2=0.0_QP
          g3=0.0_QP
          do m=1,Md
             g1 = g1 - nabla%coef(m)*( w(i1-m,i2,i3) - w(i1+m,i2,i3) )
             g2 = g2 - nabla%coef(m)*( w(i1,i2-m,i3) - w(i1,i2+m,i3) )
             g3 = g3 - nabla%coef(m)*( w(i1,i2,i3-m) - w(i1,i2,i3+m) )
          end do
          i=i+1
          grad%gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
          grad%gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
          grad%gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
       end do
       end do
       end do

       call destruct_nabla_fd( nabla )

    end if

    do i=m0,m1
       grad%gg(i) = grad%gx(i)**2 + grad%gy(i)**2 + grad%gz(i)**2
    end do

  END SUBROUTINE construct_gradient16


  SUBROUTINE destruct_gradient16( grad )
    implicit none
    type(gradient16) :: grad
    deallocate( grad%gg )
    deallocate( grad%gz )
    deallocate( grad%gy )
    deallocate( grad%gx )
  END SUBROUTINE destruct_gradient16


  subroutine calc_xyz_gradient( ixyz, rgrid, func, grad )

    implicit none
    integer,intent(in) :: ixyz
    type( grid ),intent(in) :: rgrid
    real(DP),intent(in)  :: func(:)
    real(DP),intent(out) :: grad(:)
    type(fd) :: nabla
    type(lattice) :: aa,bb
    integer :: i,i1,i2,i3,m,Md
    real(DP) :: g1,g2,g3,pi2,b(3,3),c(3)

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0d0*acos(-1.0d0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

    select case( ixyz )
    case( 1 ); c(:) = b(1,:)
    case( 2 ); c(:) = b(2,:)
    case( 3 ); c(:) = b(3,:)
    end select

! ---

    www(:,:,:,:)=0.0d0

    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
      i=i+1
      www(i1,i2,i3,1) = www(i1,i2,i3,1) + func(i)
    end do
    end do
    end do

    call bcset(1,1,Md,0)

    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
      g1=0.0d0
      g2=0.0d0
      g3=0.0d0
      do m=1,Md
        g1 = g1 - nabla%coef(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
        g2 = g2 - nabla%coef(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
        g3 = g3 - nabla%coef(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
      end do
      i=i+1
      grad(i) = c(1)*g1 + c(2)*g2 + c(3)*g3
    end do
    end do
    end do

    call destruct_nabla_fd( nabla )

  end subroutine calc_xyz_gradient


  SUBROUTINE calc_abc_gradient( iabc, rgrid, func, grad )
    implicit none
    integer,intent(IN) :: iabc
    type( grid ),intent(IN) :: rgrid
    real(DP),intent(IN)  :: func(:)
    real(DP),intent(OUT) :: grad(:)
    call init_nabla_sym
    if ( isymmetry_nabla == 1 ) then
       call calc_abc_gradient_sym( iabc, func, grad )
    else
       call calc_abc_gradient_org( iabc, rgrid, func, grad )
    end if
  END SUBROUTINE calc_abc_gradient


  SUBROUTINE calc_abc_gradient_sym( iabc, func, grad )
    implicit none
    integer,intent(IN) :: iabc
    real(DP),intent(IN)  :: func(:)
    real(DP),intent(OUT) :: grad(:)
    call op_nabla_sym( iabc, func, grad )
  END SUBROUTINE calc_abc_gradient_sym


  SUBROUTINE calc_abc_gradient_org( iabc, rgrid, func, grad )

    implicit none
    integer,intent(IN) :: iabc
    type( grid ),intent(IN) :: rgrid
    real(DP),intent(IN)  :: func(:)
    real(DP),intent(OUT) :: grad(:)
    type(fd) :: nabla
    integer :: i,i1,i2,i3,m,Md
    real(8) :: c

! ---

    call construct_nabla_fd( nabla )

    Md = nabla%md

! ---

    www(:,:,:,:)=0.0d0

    i=0
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       www(i1,i2,i3,1) = func(i)
    end do
    end do
    end do

    grad(:)=0.0d0

    select case( iabc )
    case( 1 )

       call bcset(1,1,Md,1)

       do m=1,Md
          c=nabla%coef(m)
          i=0
          do i3=rgrid%g3%z%head,rgrid%g3%z%tail
          do i2=rgrid%g3%y%head,rgrid%g3%y%tail
          do i1=rgrid%g3%x%head,rgrid%g3%x%tail
             i=i+1
             grad(i) = grad(i) - c*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          end do
          end do
          end do
       end do

    case( 2 )

       call bcset(1,1,Md,3)

       do m=1,Md
          c=nabla%coef(m)
          i=0
          do i3=rgrid%g3%z%head,rgrid%g3%z%tail
          do i2=rgrid%g3%y%head,rgrid%g3%y%tail
          do i1=rgrid%g3%x%head,rgrid%g3%x%tail
             i=i+1
             grad(i) = grad(i) - c*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          end do
          end do
          end do
       end do

    case( 3 )

       call bcset(1,1,Md,5)

       do m=1,Md
          c=nabla%coef(m)
          i=0
          do i3=rgrid%g3%z%head,rgrid%g3%z%tail
          do i2=rgrid%g3%y%head,rgrid%g3%y%tail
          do i1=rgrid%g3%x%head,rgrid%g3%x%tail
             i=i+1
             grad(i) = grad(i) - c*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
          end do
          end do
          end do
       end do

    end select

    call destruct_nabla_fd( nabla )

  END SUBROUTINE calc_abc_gradient_org


end module gradient_module
