MODULE hartree_mol_module

  use rgrid_mol_module, only: LL,KK
  use rgrid_module, only: dV,Hgrid,Ngrid
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_mol

  logical :: flag_prep = .true.
  integer :: MEO=2, Lmax_ME=4
  integer :: lmmax_ME
  real(8),allocatable :: shf1(:,:),shf2(:,:)
  real(8),allocatable :: lap(:)
  integer :: NMadv
  integer,allocatable :: adv(:),Mdv(:),Ixyz(:,:)

CONTAINS

!--------1---------2---------3---------4---------5---------6---------7--
!
! Hartree potential ( with real(8) density )
! ( The reference of this C.-G. Poisson solver is
!                  Phys. Rev. C17, 1682 (1978). )
!
  SUBROUTINE calc_hartree_mol(n1,n2,n3,rho,Vh,Eh)
    use bc_module
    use fd_module
    use kinetic_module, only: Md
    use atom_module
    implicit none
    integer,intent(IN)    :: n1,n2,n3
    real(8),intent(IN)    :: rho(n1:n2,n3)
    real(8),intent(INOUT) :: Vh(n1:n2)
    real(8),intent(OUT)   :: Eh
    integer,parameter :: maxiter = 2000
    integer :: i,j,k,k1,k2,ix,iy,iz,m,n,lm,L,a,iter,m1,m2,ierr
    integer :: ML0,M_max
    real(8) :: sum0,sum1,sum2,ak,ck,b0,const,pi,pi4,x,y,z,r,c
    real(8) :: s(4),t(4),E_0,plm(3),H
    real(8),parameter :: ep=1.d-16
    real(8),parameter :: zero=0.d0
    real(8),allocatable :: rholm(:,:),rholm_0(:,:)
    real(8),allocatable :: Vh0(:),zk0(:)
    real(8),allocatable :: tk(:),zk(:),qk(:),sk(:)
    real(8),allocatable :: tn(:)
    real(8),allocatable :: rholm0(:),rholm1(:)
    real(8),allocatable :: coef(:),clm(:,:)
    logical :: flag_alloc(2)

    if ( flag_prep ) then

       allocate( lap(-Md:Md) )
       call get_coef_laplacian_fd(Md,lap)

       lmmax_ME = (Lmax_ME+1)**2

       select case(MEO)
       case(1)
          call prep1_hartree_mol
       case(2)
          call prep2_hartree_mol
       end select

       flag_prep = .false.

    end if

    Eh  = 0.d0

    ML0 = n2-n1+1
    m1  = 1
    m2  = size(KK,2)
    pi  = acos(-1.d0)
    pi4 = 4.d0*pi
    E_0 = 0.d0
    H   = Hgrid(1)

    www(:,:,:,:) = zero


    allocate( tn(n1:n2) )
    allocate( tk(n1:n2) )

    tk(n1:n2) = rho(n1:n2,1)
    do j=2,n3
       tk(n1:n2) = tk(n1:n2) + rho(n1:n2,j)
    end do

    call mpi_allreduce(tk,tn,ML0,MPI_REAL8,MPI_SUM,comm_spin,ierr)

    sum0 = sum( tn(n1:n2)*tn(n1:n2) )*dV

    call mpi_allreduce(sum0,b0,1,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( b0 == 0.d0 ) then
       Vh(n1:n2) = 0.d0
       deallocate( tk,tn )
       return
    end if

    allocate( zk(n1:n2)  )
    allocate( qk(n1:n2)  )
    allocate( Vh0(n1:n2) )

    Vh0(n1:n2) = Vh(n1:n2)

!
! - Boundary condition 1
!
    do i=n1,n2
       www( LL(1,i),LL(2,i),LL(3,i),1 ) = Vh(i)
    end do
    call bcset(1,1,Md,0)
!
! --- Boundary condition 2 (multipole expansion) ---
!

    select case( MEO )
!
! Single center expansion
!
    case(1)

       allocate( rholm0(lmmax_ME),rholm1(lmmax_ME) )

       do lm=1,lmmax_ME
          rholm0(lm)=sum( tn(n1:n2)*shf1(n1:n2,lm) )*dV
       end do

       call mpi_allreduce &
            (rholm0,rholm1,lmmax_ME,mpi_real8,mpi_sum,comm_grid,ierr)

       do i=m1,m2
          ix=KK(1,i)
          iy=KK(2,i)
          iz=KK(3,i)
          www(ix,iy,iz,1) = sum( shf2(1:lmmax_ME,i)*rholm1(1:lmmax_ME) )
       end do

       deallocate( rholm1,rholm0 )
!
! Multi center expansion
!
    case(2)

       allocate( rholm_0(lmmax_ME,Natom),rholm(lmmax_ME,Natom) )

       rholm_0(:,:) = 0.d0
       rholm(:,:)   = 0.d0
       M_max        = Lmax_ME

       allocate( sk(0:Lmax_ME) )
       allocate( coef(0:M_max) )
       allocate( clm(0:Lmax_ME,0:M_max) )

! Pmm=(-1)^m*(2m-1)!!*(1-x^2)^(m/2)

       coef(:)=1.d0
       do m=1,M_max
          do i=1,2*m-1,2
             coef(m)=coef(m)*i
          end do
          coef(m)=coef(m)*(-1)**m
       end do

       Clm(:,:)=1.d0
       do L=0,Lmax_ME
       do m=0,L
          k1=l-m
          k2=l+m
          if ( k1 < k2 ) then
             do k=k2,k1+1,-1
                Clm(L,m)=Clm(L,m)*k
             end do
             Clm(L,m)=1.d0/Clm(L,m)
          else if ( k1 > k2 ) then
             do k=k1,k2+1,-1
                Clm(L,m)=Clm(L,m)*k
             end do
          end if
          Clm(L,m)=sqrt(Clm(L,m)*(2*l+1)/pi)*0.5d0
          if ( m/=0 ) clm(l,m)=clm(l,m)*sqrt(2.d0)
       end do
       end do

       do n=1,NMadv
          a=adv(n)
          do j=1,Mdv(n)
             i=Ixyz(j,n)
             x=LL(1,i)*H-aa_atom(1,a)
             y=LL(2,i)*H-aa_atom(2,a)
             z=LL(3,i)*H-aa_atom(3,a)
             r=sqrt(x*x+y*y+z*z)
             if ( r == 0.d0 ) then
                rholm_0(1,a)=rholm_0(1,a)+tn(i)*Clm(0,0)
                cycle
             end if
             sk(0)=tn(i)
             do L=1,Lmax_ME
                sk(L)=tn(i)*r**L
             end do
             ck=z/r ; if ( abs(ck) > 1.d0 ) ck=sign(1.d0,ck)
             if ( abs(x) < 1.d-10 ) then
                ak=0.5d0*pi
                if ( y < 0.d0 ) ak=ak+pi
             else
                ak=atan(y/x)
                if ( x < 0.d0 ) ak=ak+pi
             end if
             lm=0
             do m=0,0
                plm(1)=coef(m) !*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(0,0)*sk(0)
                plm(2)=ck*(2*m+1)*plm(1)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)+plm(2)*Clm(1,0)*sk(1)
                do L=m+2,Lmax_ME
                   plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                   lm=lm+1
                   rholm_0(lm,a)=rholm_0(lm,a)+plm(3)*Clm(L,0)*sk(L)
                   plm(1)=plm(2)
                   plm(2)=plm(3)
                end do
             end do
             do m=1,M_max-1
                plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m
                plm(2)=ck*(2*m+1)*plm(1)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)+plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)-plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m
                do L=m+2,Lmax_ME
                   plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                   lm=lm+1
                   rholm_0(lm,a)=rholm_0(lm,a)+plm(3)*Clm(L,m)*sk(L)*cos(ak*m)
                   lm=lm+1
                   rholm_0(lm,a)=rholm_0(lm,a)-plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m
                   plm(1)=plm(2)
                   plm(2)=plm(3)
                end do
             end do
             do m=M_max,M_max
                plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)
                lm=lm+1
                rholm_0(lm,a)=rholm_0(lm,a)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m
             end do
          end do ! j
       end do ! n
       rholm_0(:,:)=rholm_0(:,:)*dV
       call mpi_allreduce(rholm_0,rholm,lmmax_ME*Natom,mpi_real8,mpi_sum,comm_grid,ierr)

       do a=1,Natom
          do i=m1,m2
             ix=KK(1,i) ; iy=KK(2,i) ; iz=KK(3,i)
             x=ix*H-aa_atom(1,a)
             y=iy*H-aa_atom(2,a)
             z=iz*H-aa_atom(3,a)
             r=sqrt(x*x+y*y+z*z)
             do L=0,Lmax_ME
                sk(L)=pi4/((2.d0*L+1.d0)*r**(L+1))
             end do
             ck=z/r
             if( abs(x) < 1.d-10 )then
                ak=0.5d0*pi
                if ( y < 0.d0 ) ak=ak+pi
             else
                ak=atan(y/x)
                if ( x < 0.d0 ) ak=ak+pi
             end if
             lm=0
             do m=0,0
                plm(1)=coef(m) !*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(0,0)*sk(0)*rholm(lm,a)
                plm(2)=ck*(2*m+1)*plm(1)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(2)*Clm(1,0)*sk(1)*rholm(lm,a)
                do L=m+2,Lmax_ME
                   plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                   lm=lm+1
                   www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(3)*Clm(L,0)*sk(L)*rholm(lm,a)
                   plm(1)=plm(2)
                   plm(2)=plm(3)
                end do
             end do
             do m=1,M_max-1
                plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a)

                plm(2)=ck*(2*m+1)*plm(1)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)*rholm(lm,a)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m*rholm(lm,a)

                do L=m+2,Lmax_ME
                   plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                   lm=lm+1
                   www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(3)*Clm(L,m)*sk(L)*cos(ak*m)*rholm(lm,a)
                   lm=lm+1
                   www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m*rholm(lm,a)
                   plm(1)=plm(2)
                   plm(2)=plm(3)
                end do
             end do
             do m=M_max,M_max
                plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a)
                lm=lm+1
                www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a)
             end do
          end do ! i
       end do ! a

       deallocate( clm,coef,sk )
       deallocate( rholm,rholm_0 )

    end select ! MEO

!
! --- C.-G. minimization ---
!

    const = 3.d0*lap(0)/H**2

    zk(n1:n2) = -const*Vh(n1:n2) - pi4*tn(n1:n2)

    do j=1,Md
       c=lap(j)/H**2
       do i=n1,n2
          ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
          zk(i)=zk(i)-c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
                         +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
                         +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
       end do
    end do

    www(:,:,:,1)=0.d0

    qk(:)=zk(:)
    sum0=sum(zk(n1:n2)*zk(n1:n2))*dV

    call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)

!
! --- Iteration ---
!
    Iteration : do iter=1,maxiter

       do i=n1,n2
          ix=LL(1,i)
          iy=LL(2,i)
          iz=LL(3,i)
          www(ix,iy,iz,1)=qk(i)
       end do

       call bcset(1,1,Md,0)

       tk(n1:n2) = const*qk(n1:n2)

       do j=1,Md
          c=lap(j)/H**2
          do i=n1,n2
             ix=LL(1,i)
             iy=LL(2,i)
             iz=LL(3,i)
             tk(i)=tk(i)+c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
     &                      +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
     &                      +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
          end do
       end do

       sum0 = sum( zk(n1:n2)*tk(n1:n2) )*dV

       call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)

       Vh0(n1:n2) = Vh(n1:n2)

       ak=sum1/sum2
       Vh(n1:n2) = Vh(n1:n2) + ak*qk(n1:n2)
       zk(n1:n2) = zk(n1:n2) - ak*tk(n1:n2)

!
! Conv. Check
!
       E_0  = Eh

       s(1) = sum( zk(n1:n2)*zk(n1:n2) )*dV
       s(2) = 0.5d0*sum( Vh(n1:n2)*tn(n1:n2) )*dV
       s(3) = sum( (Vh(n1:n2)-Vh0(n1:n2))**2 )/Ngrid(0)
       s(4) = maxval( abs(Vh(n1:n2)-Vh0(n1:n2)) )

       call mpi_allreduce(s,t,4,mpi_real8,mpi_sum,comm_grid,ierr)

       sum2 = t(1)
       Eh   = t(2)
       if ( t(4) <= ep ) exit

       ck=sum2/sum1
       sum1=sum2

       qk(n1:n2) = ck*qk(n1:n2) + zk(n1:n2)

    end do Iteration

    deallocate( Vh0 )
    deallocate( qk )
    deallocate( zk )
    deallocate( tk )
    deallocate( tn )

    www(:,:,:,:) = zero

    return

  END SUBROUTINE calc_hartree_mol

!--------1---------2---------3---------4---------5---------6---------7--
!
! Spherical Harmonic Funciton
!
  SUBROUTINE prep1_hartree_mol
    integer :: i,j,k,lm,L,a,m,n,n1,n2,ML0
    integer :: m1,m2
    logical :: flag_alloc(2)
    real(8),parameter :: eps=1.d-20
    real(8) :: r,x,y,z,const,H,pi4

    INTERFACE
       FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
       END FUNCTION Ylm
    END INTERFACE

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)
    m1  = 1
    m2  = size(KK,2)
    H   = Hgrid(1)
    pi4 = 4.d0*acos(-1.d0)

    allocate( shf1(n1:n2,lmmax_ME) ) ; shf1=0.d0
    allocate( shf2(lmmax_ME,m1:m2) ) ; shf2=0.d0

    do i=n1,n2
       x=LL(1,i)*H ; y=LL(2,i)*H ; z=LL(3,i)*H
       r=sqrt(x*x+y*y+z*z)
       lm=0
       do L=0,Lmax_ME
          do m=-L,L
             lm=lm+1
             shf1(i,lm)=Ylm(x,y,z,L,m)*r**L
          end do
       end do
    end do

    do i=m1,m2
       x=KK(1,i)*H ; y=KK(2,i)*H ; z=KK(3,i)*H
       r=sqrt(x*x+y*y+z*z)
       lm=0
       do L=0,Lmax_ME
          const=pi4/(2.d0*L+1.d0)
          do m=-L,L
             lm=lm+1
             shf2(lm,i)=Ylm(x,y,z,L,m)/r**(L+1)*const
          end do
       end do
    end do

    return
  END SUBROUTINE prep1_hartree_mol

!--------1---------2---------3---------4---------5---------6---------7--
!
! Division of Area for Multipole Expansion (MEO=2)
!
  SUBROUTINE prep2_hartree_mol
    use atom_module
    real(8),allocatable :: ra(:)
    real(8) :: x,y,z,r2,H
    integer :: i,a,m,n,i1,i2,i3,n1,n2,ierr,maxMdv
    integer,allocatable :: itmp(:),jtmp(:)

    n1 = idisp(myrank)+1
    n2 = idisp(myrank)+ircnt(myrank)
    H  = Hgrid(1)

    allocate( itmp(n1:n2) ) ; itmp=0
    allocate( jtmp(Natom) ) ; jtmp=0

    allocate( ra(n1:n2) )
    ra=1.d10
    do a=1,Natom
       do i=n1,n2
          x=LL(1,i)*H-aa_atom(1,a)
          y=LL(2,i)*H-aa_atom(2,a)
          z=LL(3,i)*H-aa_atom(3,a)
          r2=x*x+y*y+z*z
          if ( r2 < ra(i) ) then
             ra(i)=r2
             itmp(i)=a
          end if
       end do
    end do
    deallocate( ra )

    do a=1,Natom
       jtmp(a)=count(itmp==a) ! # of grid points near the atom a
    end do

    NMadv=count( jtmp>0 )     ! # of atoms (or # of regions) in my rank
    maxMdv=maxval(jtmp)       ! max # of grid points around each atom

!    if ( DISP_SWITCH ) then
       write(*,*) "NMadv,maxMdv,ML,sum(jtmp)=",NMadv,maxMdv,Ngrid(0),sum(jtmp)
!    end if

    allocate( Ixyz(maxMdv,NMadv) )
    allocate( Mdv(NMadv) )
    allocate( adv(NMadv) )

    n=0
    do a=1,Natom
       if ( jtmp(a) > 0 ) then
          n=n+1
          Mdv(n)=jtmp(a)
          adv(n)=a
       end if
    end do

    if ( n /= NMadv ) stop "n/=NMadv!!!"

    do n=1,NMadv
       a=adv(n)
       m=0
       do i=n1,n2
          if ( itmp(i) == a ) then
             m=m+1
             Ixyz(m,n)=i
          end if
       end do
    end do

    deallocate( jtmp,itmp )

    return
  END SUBROUTINE prep2_hartree_mol


END MODULE hartree_mol_module
