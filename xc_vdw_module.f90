MODULE xc_vdw_module

  use grid_module, only: grid
  use xc_variables, only: xc
  use density_module, only: density
  use ggrid_module, only: NMGL, MGL, GG, NGgrid, LLG, MG_0, MG_1 &
                         ,construct_ggrid, destruct_ggrid
  use gradient_module, only: gradient, construct_gradient, destruct_gradient
  use parallel_module, only: comm_grid

  implicit none

  PRIVATE
  PUBLIC :: init_xc_vdw, calc_xc_vdw

  real(8) :: Qmax=5.0d0
  real(8) :: Qmin=0.0d0
  integer :: NumQgrid=20
  real(8) :: SizeQgrid
  real(8),allocatable :: Qgrid(:)
  real(8) :: Dmax=48.6d0

  real(8),allocatable :: phiG(:,:,:)
  real(8),allocatable :: bmat(:,:),cmat(:,:),dmat(:,:)

  real(8),parameter :: zero_density = 1.d-10
  logical :: flag_init=.false.

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_xc_vdw
    implicit none
    integer :: i,j

    if ( flag_init ) return

    call write_border(30," init_xc_vdw(start)")

    SizeQgrid = Qmax / NumQgrid
    allocate( Qgrid(0:NumQgrid) ) ; Qgrid=0.0d0
    do i=1,NumQgrid
       Qgrid(i) = Qmin + i*SizeQgrid
    end do

    allocate( phiG(NMGL,0:NumQgrid,0:NumQgrid) )
    phiG=0.0d0

    do j=1,NumQgrid
    do i=j,NumQgrid

       call calc_phiG( Qgrid(i), Qgrid(j), phiG(1,i,j) )

       if ( i /= j ) phiG(:,j,i) = phiG(:,i,j)

    end do ! i
    end do ! j

    call init1_xc_vdw ! coefficients of 3rd spline

    flag_init = .true.

    call write_border(30," init_xc_vdw(end)")

  END SUBROUTINE init_xc_vdw


  SUBROUTINE calc_phiG( qi, qj, phiG )
    implicit none
    real(8),intent(IN)  :: qi,qj
    real(8),intent(OUT) :: phiG(NMGL)
    integer :: nr,ig,i
    real(8) :: rmax,dr,r,pi4,G,sum0
    real(8),allocatable :: phi(:)

    rmax = Dmax/max( qi, qj )
    nr   = 400
    dr   = rmax/nr

    pi4  = 4.0d0*acos(-1.0d0)

    allocate( phi(0:nr) ) ; phi=0.0d0

    do i=1,nr
       r = i*dr
       call calc_phi( qi*r, qj*r, phi(i) )
    end do ! ir

    phiG=0.0d0

    do ig=1,NMGL

       G=sqrt( GG(ig) )

       if ( G <= 1.d-9 ) then

          sum0=0.0d0
          do i=1,nr
             r = i*dr
             sum0 = sum0 + r*r*phi(i)
          end do
          phiG(ig)=pi4*sum0*dr

       else

          sum0=0.0d0
          do i=1,nr
             r = i*dr
             sum0 = sum0 + r*phi(i)*sin(G*r)
          end do
          phiG(ig)=pi4/G*sum0*dr

       end if

    end do ! ig

    deallocate( phi )

  END SUBROUTINE calc_phiG


  SUBROUTINE calc_phi( di, dj, phi )
    implicit none
    real(8),intent(IN)  :: di,dj
    real(8),intent(OUT) :: phi
    real(8) :: rmax,dr,dt,r,t,a,b,gamma,pi,WWW,TTT
    real(8) :: ct,st,sa,sb,ca,cb,hai,haj,hbi,hbj,vai,vaj,vbi,vbj
    integer :: nr, nt, ir, it

    pi = acos(-1.0d0)
    gamma = 4.0d0*pi/9.0d0

    nr   = 100
    rmax = 20.0d0
    dr   = rmax/nr

    nt = 9 ! to use a symmetry(a <-> b), odd number should be given
    dt = 0.5d0*pi / nt

    phi = 0.0d0

    do it=1,(nt-1)/2
       t  = it*dt
       ct = cos(t)
       st = sin(t)
    do ir=0,nr
       r = ir*dr
       a = r*ct
       b = r*st

       ca = cos(a)
       cb = cos(b)
       sa = 1.0d0 ; if ( a /= 0.0d0 ) sa=sin(a)/a
       sb = 1.0d0 ; if ( b /= 0.0d0 ) sb=sin(b)/b

       WWW = 2.0d0*( (3.0d0-a*a)*cb*sa + (3.0d0-b*b)*ca*sb &
                   + (r*r-3.0d0)*sa*sb - 3.0d0*ca*cb )

       if ( di == 0.0d0 ) then
          hai = 0.0d0
          hbi = 0.0d0
       else
          hai = 1.0d0 - exp( -gamma*(a/di)**2 )
          hbi = 1.0d0 - exp( -gamma*(b/di)**2 )
       end if
       if ( dj == 0.0d0 ) then
          haj = 0.0d0
          hbj = 0.0d0
       else
          haj = 1.0d0 - exp( -gamma*(a/dj)**2 )
          hbj = 1.0d0 - exp( -gamma*(b/dj)**2 )
       end if

       if ( a == 0.0d0 ) then
          vai = di*di/(2.0d0*gamma)
          vaj = dj*dj/(2.0d0*gamma)
       else
          vai = a*a/(2.0d0*hai)
          vaj = a*a/(2.0d0*haj)
       end if
       if ( b == 0.0d0 ) then
          vbi = di*di/(2.0d0*gamma)
          vbj = dj*dj/(2.0d0*gamma)
       else
          vbi = b*b/(2.0d0*hbi)
          vbj = b*b/(2.0d0*hbj)
       end if

       TTT = 0.5d0*( 1.0d0/(vai+vbi) + 1.0d0/(vaj+vbj) ) &
                  *( 1.0d0/( (vai+vaj)*(vbi+vbj) ) &
                    +1.0d0/( (vai+vbj)*(vaj+vbi) ) )

       phi = phi + r*WWW*TTT

    end do ! ir
    end do ! it

    phi = 2.0d0*(2.0d0/pi**2)*phi*dt*dr

  END SUBROUTINE calc_phi


  SUBROUTINE init1_xc_vdw
    implicit none
    real(8),allocatable :: mat1(:,:),mat2(:,:),work(:)
    integer,allocatable :: ipiv(:)
    real(8) :: h
    integer :: m,n,i,info

    n = NumQgrid
    h = SizeQgrid
    m = n+1

    allocate( bmat(0:n,0:n) ) ; bmat=0.0d0
    allocate( cmat(0:n,0:n) ) ; cmat=0.0d0
    allocate( dmat(0:n,0:n) ) ; dmat=0.0d0

    allocate( mat1(0:n,0:n) ) ; mat1=0.0d0
    allocate( mat2(0:n,0:n) ) ; mat2=0.0d0
    allocate( ipiv(0:n)     ) ; ipiv=0
    allocate( work(64*m)    ) ; work=0.0d0

    mat1(0,0)=1.0d0
    do i=1,n-1
       mat1(i,i-1)=h
       mat1(i,i  )=4.0d0*h
       mat1(i,i+1)=h
    end do
    mat1(n,n)=1.0d0

    do i=1,n-1
       mat2(i,i-1)= 1.0d0
       mat2(i,i  )=-2.0d0
       mat2(i,i+1)= 1.0d0
    end do
    mat2(0:n,0:n)=mat2(0:n,0:n)*3.0d0/h

    call dgetrf( m,m,mat1,m,ipiv,info)
    call dgetri( m,mat1,m,ipiv,work,size(work),info)

    cmat(0:n,0:n) = matmul( mat1(0:n,0:n), mat2(0:n,0:n) )

    mat1=0.0d0
    do i=0,n-1
       mat1(i,i  )=-1.0d0
       mat1(i,i+1)= 1.0d0
    end do

    dmat(0:n,0:n) = matmul( mat1(0:n,0:n), cmat(0:n,0:n) )/(3.0d0*h)

    bmat(0:n,0:n) = mat1(0:n,0:n)/h - cmat(0:n,0:n)*h - dmat(0:n,0:n)*h**2

    deallocate( work )
    deallocate( ipiv )
    deallocate( mat2 )
    deallocate( mat1 )

  END SUBROUTINE init1_xc_vdw


  SUBROUTINE calc_xc_vdw( rgrid, rho, vdw )
    implicit none
    type(grid) :: rgrid
    type(density) :: rho
    type(xc) :: vdw
    type(gradient) :: grad
    logical :: disp_sw
    integer :: i,j,m0,m1,mm
    real(8) :: rho_tmp(2), trho, edx_lda(2), edc_lda
    real(8) :: kf,pi,onethr,fouthr,ss
    real(8),allocatable :: q0(:), pol_spline(:,:)
    real(8),parameter :: Zab=-0.8491d0
    complex(8),allocatable :: theta(:,:)

    call write_border(30," calc_xc_vdw(start)")

    pi     = acos(-1.0d0)
    onethr = 1.0d0/3.0d0
    fouthr = 4.0d0/3.0d0
    m0     = rgrid%SubGrid(1,0)
    m1     = rgrid%SubGrid(2,0)

    call construct_gradient( rgrid, rho, grad )

    allocate( q0(m0:m1) ) ; q0=0.0d0

    do i=m0,m1

       rho_tmp(1) = rho%rho(i,1)
       rho_tmp(2) = rho%rho(i,rho%nn)
       trho       = sum(rho_tmp)

       call calc_edx_lda( rho%nn, rho_tmp, edx_lda )
       call calc_edc_lda( rho%nn, rho_tmp, edc_lda )

       kf = ( 3.0d0*pi*pi*trho )**onethr

       ss = grad%gg(i)/(2.0d0*kf*trho)**2

       q0(i) = -fouthr*pi*(edx_lda(1)+edc_lda) - Zab/9.0d0*ss*kf

    end do ! i

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) write(*,*) "q0(min,max)=",minval(q0),maxval(q0)

! ---

    allocate( pol_spline(m0:m1,0:NumQgrid) ) ; pol_spline=0.0d0

    call calc_pol_spline( m0, m1, q0, pol_spline )

! ---

    mm = max( m1-m0+1, MG_1-MG_0+1 )

    allocate( theta(mm,0:NumQgrid) ) ; theta=(0.0d0,0.0d0)

    do j=0,NumQgrid
       do i=m0,m1
          theta(i-m0+1,j) = sum( rho%rho(i,1:rho%nn) )*pol_spline(i,j)
       end do
    end do

    call fft_theta( rgrid, theta )

! ---

    call calc_vdw_energy( rgrid, theta, vdw%Ec )

! ---

    deallocate( theta )
    deallocate( pol_spline )
    deallocate( q0 )

    call destruct_gradient( grad )

    call write_border(30," calc_xc_vdw(end)")

  END SUBROUTINE calc_xc_vdw


  SUBROUTINE calc_edx_lda( nn, rho, edx_lda )
    implicit none
    integer,intent(IN)  :: nn
    real(8),intent(IN)  :: rho(nn)
    real(8),intent(OUT) :: edx_lda(nn)
    real(8) :: onethr, thrfou, thrpi, factor
    integer :: s

    onethr = 1.0d0/3.0d0
    thrfou = 3.0d0/4.0d0
    thrpi  = 3.0d0/acos(-1.0d0)

    factor = 1.0d0
    if ( nn == 2 ) factor = 2.0d0**onethr

    do s=1,nn
       edx_lda(s) = -factor*thrfou*( thrPi*rho(s) )**onethr
    end do

  END SUBROUTINE calc_edx_lda


  SUBROUTINE calc_edc_lda( nn, rho, edc_lda )
    implicit none
    integer,intent(IN)  :: nn
    real(8),intent(IN)  :: rho(nn)
    real(8),intent(OUT) :: edc_lda
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
    real(8) :: one,two,fouthr,pi,const1,const2,onethr,ThrFouPi
    real(8) :: trho,rhoa,rhob,zeta,fz,kf,rs,rssq,ec_U,ec_P,alpc

    edc_lda = 0.0d0

    trho = sum( rho(1:nn) )
    rhoa = rho(1)
    rhob = rho(nn)

    if ( trho <= zero_density ) return

    one      = 1.0d0
    two      = 2.0d0
    fouthr   = 4.0d0/3.0d0
    onethr   = 1.0d0/3.0d0
    pi       = acos(-1.0d0)
    ThrFouPi = 3.0d0/(4.0d0*pi)
    const1   = two**fouthr-two
    const2   = 9.0d0/4.0d0*(two**onethr-one)

    zeta = ( rhoa - rhob )/trho
    fz   = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1
    kf   = ( 3.0d0*Pi*Pi*trho )**onethr
    rs   = ( ThrFouPi/trho )**onethr
    rssq = sqrt(rs)
    ec_U = -two*A00*( one + alp10*rs ) &
         *log( one + one/( two*A00 &
         *(bt10*rssq + bt20*rs + bt30*rs*rssq + bt40*rs*rs) ) )
    ec_P = -two*A01*( one + alp11*rs ) &
         *log( one + one/( two*A01 &
         *(bt11*rssq + bt21*rs + bt31*rs*rssq + bt41*rs*rs) ) )
    alpc = -two*A02*( one + alp12*rs ) &
         *log( one + one/( two*A02 &
         *(bt12*rssq + bt22*rs + bt32*rs*rssq + bt42*rs*rs) ) )

    edc_lda = ec_U - alpc*fz*const2*(one-zeta**4) &
         + ( ec_P - ec_U )*fz*zeta**4

  END SUBROUTINE calc_edc_lda


  SUBROUTINE calc_pol_spline( m0, m1, q0, pol )
    implicit none
    integer,intent(IN)  :: m0,m1
    real(8),intent(IN)  :: q0(m0:m1)
    real(8),intent(OUT) :: pol(m0:m1,0:NumQgrid)
    integer :: i,j,k
    real(8) :: q

    pol(:,:)=0.0d0

    do k=m0,m1

       do i=0,NumQgrid
          if ( q0(k) < Qgrid(i) ) exit
       end do
       i=i-1

       q=q0(k)-Qgrid(i)
       do j=0,NumQgrid
          pol(k,j) = bmat(i,j)*q + cmat(i,j)*q**2 + dmat(i,j)*q**3
       end do
       pol(k,i) = pol(k,i) + 1.0d0

! ---------------------------- Wu & Gygi ---
       do j=0,NumQgrid
          pol(k,j) = Qgrid(j)*pol(k,j)/q0(k)
       end do
! ------------------------------------------

    end do ! k

  END SUBROUTINE calc_pol_spline


  SUBROUTINE fft_theta( rgrid, theta )
    implicit none
    type(grid),intent(IN) :: rgrid
    complex(8),intent(INOUT) :: theta(:,0:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    integer :: ML,ML1,ML2,ML3,mm,i1,i2,i3,i,j,info
    complex(8),allocatable :: work0(:,:,:), work1(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    ML  = rgrid%NumGrid(0)
    ML1 = rgrid%NumGrid(1)
    ML2 = rgrid%NumGrid(2)
    ML3 = rgrid%NumGrid(3)

    allocate( work0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work0=zero
    allocate( work1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; work1=zero

    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call construct_ggrid(1)

    do j=0,NumQgrid

       work1(:,:,:)=zero
       i=0
       do i3=rgrid%SubGrid(1,3),rgrid%SubGrid(2,3)
       do i2=rgrid%SubGrid(1,2),rgrid%SubGrid(2,2)
       do i1=rgrid%SubGrid(1,1),rgrid%SubGrid(2,1)
          i=i+1
          work1(i1,i2,i3) = theta(i,j)
       end do ! i1
       end do ! i2
       end do ! i3

       call MPI_ALLREDUCE( work1, work0, size(work0) &
            , MPI_COMPLEX16, MPI_SUM, comm_grid, info )

       call fft3fx(ML1,ML2,ML3,ML,work0,work1,wsavex,wsavey,wsavez &
            ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

       theta(:,j) = zero
       do i=MG_0,MG_1
          i1=mod( LLG(1,i)+ML1, ML1 )
          i2=mod( LLG(2,i)+ML2, ML2 )
          i3=mod( LLG(3,i)+ML3, ML3 )
          theta(i-MG_0+1,j) = work0(i1,i2,i3)
       end do

    end do ! j

    call destruct_ggrid

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( work1 )
    deallocate( work0 )

  END SUBROUTINE fft_theta


  SUBROUTINE calc_vdw_energy( rgrid, theta, Ec )
    implicit none
    complex(8),intent(IN) :: theta(:,0:)
    type(grid),intent(IN) :: rgrid
    real(8),intent(OUT) :: Ec
    integer :: a,b,i,j,info
    complex(8) :: z
    real(8) :: c,E0

    z=(0.0d0,0.0d0)

    do b=0,NumQgrid
    do a=b,NumQgrid

       c=1.0d0 ; if ( a /= b ) c=2.0d0

       do i=MG_0,MG_1
          j=i-MG_0+1
          z = z + c*theta(j,a)*conjg(theta(j,b))*phiG(MGL(i),a,b)
       end do ! i

    end do ! a
    end do ! b

    E0 = 0.5d0*z * rgrid%NumGrid(0)*rgrid%dV
    call MPI_ALLREDUCE(E0,Ec,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,info)

    write(*,*) "Ec_vdw=",Ec

  END SUBROUTINE calc_vdw_energy


END MODULE xc_vdw_module
