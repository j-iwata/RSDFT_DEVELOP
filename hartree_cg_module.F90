MODULE hartree_module

  use hartree_mol_module

  implicit none

  PRIVATE
  PUBLIC :: E_hartree,Vh,calc_hartree

  real(8) :: E_hartree
  real(8),allocatable :: Vh(:)
  integer :: SYStype=0

CONTAINS


  SUBROUTINE calc_hartree(n1,n2,n3,rho,SYStype_in)
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer,optional,intent(IN) :: SYStype_in

    if ( .not.allocated(Vh) ) then
       allocate( Vh(n1:n2) )
       Vh=0.d0
    end if

    if ( present(SYStype_in) ) then
       SYStype=SYStype_in
    end if

    select case(SYStype)
    case default
       call calc_hartree_sol_cg(n1,n2,n3,rho)
    case(1)
       call calc_hartree_mol(n1,n2,n3,rho,Vh,E_hartree)
    end select

  END SUBROUTINE calc_hartree


  SUBROUTINE calc_hartree_sol_cg(n1,n2,n3,rho)
    use rgrid_module, only: Ngrid,Igrid,dV
    use ggrid_module
    use parallel_module
    use electron_module
    use kinetic_module
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
    integer,parameter :: max_iter=2000
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,ispin
    integer :: a1,a2,a3,b1,b2,b3
    real(8),parameter :: tol=1.d-24
    real(8) :: Eh0,pi4,gg_0,gg_1,sum0,bet,alp
    real(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork(:),zwork1(:)
    real(8),allocatable :: w(:,:,:)
    integer :: ML1,ML2,ML3,ML,iter
    real(8),allocatable :: rhot(:),q(:),b(:),p(:),g(:),x(:)
    real(8) :: c,Vcell,pi2,sum2,sum1,sb(2),rb(2)

    pi2 = 2.d0*acos(-1.d0)
    pi4 = 4.d0*acos(-1.d0)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    a1=Igrid(1,1)
    b1=Igrid(2,1)
    a2=Igrid(1,2)
    b2=Igrid(2,2)
    a3=Igrid(1,3)
    b3=Igrid(2,3)

    Vcell = ML*dV

    allocate( rhot(n1:n2) )
    rhot(n1:n2) = rho(n1:n2,1)
    do ispin=2,n3
       rhot(n1:n2) = work(n1:n2) + rho(n1:n2,ispin)
    end do

    allocate( x(n1:n2) )
    allocate( b(n1:n2) )
    allocate( g(n1:n2) )
    allocate( p(n1:n2) )
    allocate( q(n1:n2) )

    c=pi4*Nelectron/Vcell
    b(:) = -pi4*rhot(:)+c

    x(:) = Vh(:)

    allocate( zwork(n1:n2),zwork1(n1:n2) )
    zwork(:)=x(:)
    zwork1(:)=0.d0

    call op_kinetic(0,zwork,zwork1,n1,n2,1,1)

    g(n1:n2) = b(n1:n2) - (-2.d0)*zwork1(n1:n2)

    p(n1:n2) = g(n1:n2)

!    sum0=sum( g(n1:n2)*g(n1:n2) )*dV
!    call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)

    gg_1 = 1.d10
    sum1 = 1.d10
    sum2 = 1.d10

    do iter=1,max_iter

!       zwork(:)=p(:)
!       zwork1(:)=0.d0
!       call op_kinetic(0,zwork,zwork1,n1,n2,1,1)
!       q(:) = (-2.d0)*zwork1(:)

!       sum0 = sum( g(n1:n2)*q(n1:n2) )*dV
!       call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)

!       alp = sum1/sum2

       gg_0 = gg_1

       sum0 = sum( g(n1:n2)*g(n1:n2) )
       call mpi_allreduce(sum0,gg_1,1,mpi_real8,mpi_sum,comm_grid,ierr)

!       if ( myrank==0 ) write(*,'(1x,i6,3g15.5)') iter,gg_1,sum1,sum2

!       if ( gg_1 < tol ) exit
       if ( sum1 < tol ) exit

       if ( iter > 1 ) then
          bet=gg_1/gg_0
          p(n1:n2)=g(n1:n2)+bet*p(n1:n2)
       end if

       zwork(:)=p(:)
       zwork1(:)=0.d0
       call op_kinetic(0,zwork,zwork1,n1,n2,1,1)
       q(:) = (-2.d0)*zwork1(:)

       sum0 = sum( p(n1:n2)*q(n1:n2) )
       call mpi_allreduce(sum0,alp,1,mpi_real8,mpi_sum,comm_grid,ierr)
       alp=gg_1/alp

       Vh(:)=x(:)

       x(n1:n2) = x(n1:n2) + alp*p(n1:n2)

       g(n1:n2) = g(n1:n2) - alp*q(n1:n2)

       sb(1)=sum( (Vh-x)**2 )/ML
       sb(2)=0.5d0*sum(x*rhot)*dV
       call mpi_allreduce(sb,rb,2,mpi_real8,mpi_sum,comm_grid,ierr)
       sum1=rb(1)
       sum2=rb(2)

!       sum0 = sum( g(n1:n2)*g(n1:n2) )*dV
!       call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)
!       if ( myrank==0 ) write(*,*) iter,sum2,sum(x*rhot)*dV

!       if ( sum2 < tol ) exit

!       bet=sum2/sum1
!       sum1=sum2

!       p(n1:n2) = bet*p(n1:n2) + g(n1:n2)

    end do

    Vh(n1:n2) = x(n1:n2)

    Eh0=0.5d0*sum( Vh(:)*rhot(:) )*dV
    call mpi_allreduce(Eh0,E_hartree,1,mpi_real8,mpi_sum,comm_grid,ierr)

    deallocate( zwork1 )
    deallocate( zwork  )
    deallocate( q,p,g,b,x )
    deallocate( rhot )

  END SUBROUTINE calc_hartree_sol_cg


END MODULE hartree_module
