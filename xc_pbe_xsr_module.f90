MODULE xc_pbe_xsr_module

  use gradient_module
  use grid_module, only: grid, get_map_3d_to_1d
  use xc_variables, only: xcene, xcpot
  use fd_module, only: fd, construct_nabla_fd
  use lattice_module
  use parallel_module
  use expint_module
  use gradient_module
  use BasicTypeFactory

  implicit none

  PRIVATE
  PUBLIC :: calc_pbe_xsr

  integer,parameter :: DP=kind(0.0d0)
  integer,parameter :: QP=kind(0.0q0)

  real(DP),parameter :: zero_density = 1.d-10
  real(DP),allocatable :: nab(:)
  real(DP),allocatable :: vx(:,:)
  real(DP) :: Ex
  real(DP) :: q(3,3)
  integer,allocatable :: LLL(:,:,:)
  integer :: Md, ML1,ML2,ML3
  integer :: SYStype=0

  real(DP) :: mu=0.21951d0
  real(DP) :: Kp=0.804d0
  real(DP) :: beta=0.066725d0
  real(DP) :: omega=0.11d0    !(HSE06)

CONTAINS


  SUBROUTINE calc_pbe_xsr( rgrid, rho, ene, pot )

    implicit none

    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN) :: rho
    type( xcene ) :: ene
    type( xcpot ) :: pot

    type(gradient) :: grad
    type(fd) :: nabla
    type(lattice) :: aa, bb
    integer :: mm,m1,m2,n1,n2,i,i1,i2,i3,s,j,ierr
    integer :: m,k1,k2,k3,j1,j2,j3
    real(DP) :: Pi, cm
    real(DP) :: sb(1),rb(1)

    real(DP),parameter :: A = 1.0161144d0
    real(DP),parameter :: B =-0.37170836d0
    real(DP),parameter :: C =-0.077215461d0
    real(DP),parameter :: D = 0.57786348d0
    real(DP),parameter :: E =-0.051955731d0
    real(DP),parameter :: ah1=0.00979681d0, ah2=0.0410834d0
    real(DP),parameter :: ah3=0.187440d0, ah4=0.00120824d0, ah5=0.0347188d0
    real(DP),parameter :: EGc1=-0.02628417880d0,EGc2=-0.07117647788d0
    real(DP),parameter :: EGc3=0.08534541323d0
    real(DP),parameter :: EGscut=0.08d0

    real(DP) :: ss,s1,kf,vx_lda,ex_lda,factor,trho
    real(DP) :: Fs, Gs, Hs, aG, bG, EGs
    real(DP) :: dHds, dFds, dEGds, daGds, dbGds
    real(DP) :: cc(0:9),dd(0:11),Iy(32),Jy(0:9),Igauss(0:15)
    real(DP) :: A1,H1,H2,H3,erfc_ha,FF,GG,Fx,dFxds,dFxdn,x,expint_ha
    real(DP) :: DHs,DHs2,DHs3,dsdDHs,srpi,DHs92,A12
    real(DP),allocatable :: rtmp(:),rrrr(:,:)

    INTERFACE
       FUNCTION bberf(x)
         real(8) :: bberf,x
       END FUNCTION bberf
    END INTERFACE

! ---

    beta = mu*3.0d0/acos(-1.0d0)**2

    call construct_gradient( rgrid, rho, grad )

    call construct_nabla_fd( nabla )

    Md = nabla%md

    if ( .not.allocated(nab) ) allocate( nab(-Md:Md) )
    nab(:) = nabla%coef(:)

    call construct_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )

    Pi=acos(-1.0d0)

    q(1:3,1)=aa%Length(1)*bb%LatticeVector(1:3,1)/( 2*Pi*rgrid%spacing(1) )
    q(1:3,2)=aa%Length(2)*bb%LatticeVector(1:3,2)/( 2*Pi*rgrid%spacing(2) )
    q(1:3,3)=aa%Length(3)*bb%LatticeVector(1:3,3)/( 2*Pi*rgrid%spacing(3) )

    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    call get_map_3d_to_1d( LLL )

! ---

    cc(0) = 1.0d0
    cc(1) =-1.128223947d0
    cc(2) = 1.452736266d0
    cc(3) =-1.243162299d0
    cc(4) = 0.971824836d0
    cc(5) =-0.568861080d0
    cc(6) = 0.246880515d0
    cc(7) =-0.065032364d0
    cc(8) = 0.008401793d0
    cc(9) = 1.455915450d0

    cc(1) =-1.128223946706117d0
    cc(2) = 1.452736265762971d0
    cc(3) =-1.243162299390327d0
    cc(4) = 0.971824836115601d0
    cc(5) =-0.568861079687373d0
    cc(6) = 0.246880514820192d0
    cc(7) =-0.065032363850763d0
    cc(8) = 0.008401793031216d0
    cc(9) = 1.455915450052607d0

! ---

    m1 = pot%xc%g_range%head
    m2 = pot%xc%g_range%tail
    n1 = pot%xc%s_range%head
    n2 = pot%xc%s_range%tail
    mm = pot%xc%g_range%size_global

    allocate( vx(m1:m2,n1:n2) ) ; vx=0.0d0
    allocate( rtmp(m1:m2) ) ; rtmp=0.0d0
    allocate( rrrr(mm,3) ) ; rrrr=0.0d0

    Ex = 0.0d0

    factor = 1.0d0
    if ( rho%s_range%size_global == 2 ) factor = 2.0d0

    do s=n1,n2

       do i=m1,m2

          trho = factor*rho%val(i,s)
          if ( trho <= zero_density ) cycle

          kf = (3.0d0*Pi*Pi*trho)**(1.0d0/3.0d0)

          ex_lda = -3.0d0/(4.0d0*Pi)*kf
          vx_lda = -1.0d0/Pi*kf

          ss = grad%gg(i)/( 2.0d0*trho*kf )**2
          s1 = sqrt(ss)
          if ( s1 > 8.3d0 ) then
             s1 = 8.572844d0 - 18.796223d0/ss
             ss = s1*s1
          end if

          Hs = (ah1*ss+ah2*ss**2)/(1.0d0+ah3*ss**2+ah4*s1*ss**2+ah5*ss**3)

!         Fs = ( Hs*(16.d0*A*A+36.d0*(B-A*D))+9.d0*(-4.d0/27.d0) )/(36.d0*C)
          Fs = 6.4753871d0*Hs + 0.47965830d0

          aG = sqrt(Pi)/(16.d0*(D+Hs*ss)**3.5d0)*( 15.d0*E &
               +6.d0*C*(1.d0+Fs*ss)*(D+Hs*ss)+4.d0*B*(D+Hs*ss)**2 &
               +8.d0*A*(D+Hs*ss)**3 ) &
               -3.d0*Pi*sqrt(A)/4.d0*exp(9.d0*Hs*ss/(4.d0*A)) &
               *(1.0d0-bberf(1.5d0*s1*sqrt(Hs/A)))
          bG = 15.d0*sqrt(Pi)*ss/(16.d0*(D+Hs*ss)**3.5d0)

          Gs  = -(0.75d0*Pi+aG)/(bG*E)
          EGs = -(0.75d0*Pi+aG)/bG
          if ( s1 <= EGscut ) EGs = EGc1 + EGc2*ss + EGc3*ss*ss

!          dHds = ( (2*ah1*s1+4*ah2*ss*s1)*(1+ah3*ss*ss+ah4*ss*ss*s1+ah5*ss**3) &
!                  -(ah1*ss+ah2*ss*ss)*(4*ah3*ss*s1+5*ah4*ss*ss+6*ah5*ss*ss*s1) ) &
!                /(1+ah3*ss*ss+ah4*ss*ss*s1+ah5*ss**3)**2

          dHds=((2.d0*ah1*s1+4.d0*ah2*ss*s1)*(1.d0+ah3*ss*ss+ah4*ss*ss*s1+ah5*ss**3) &
               -(ah1*ss+ah2*ss*ss)*(4.d0*ah3*ss*s1+5.d0*ah4*ss*ss+6.d0*ah5*ss*ss*s1)) &
               /(1.d0+ah3*ss*ss+ah4*ss*ss*s1+ah5*ss**3)**2.d0

          dFds = 6.4753871d0*dHds

!          daGds = sqrt(Pi)*7.0d0/16.0d0*Hs*s1/(D+ss*Hs)**4.5d0 &
!                *(15*E+6*C*(1+ss*Fs)*(D+ss*Hs)+4*B*(D+ss*Hs)**2+8*A*(D+ss*Hs)**3) &
!                +0.25d0*s1*sqrt(Pi)/(D+ss*Hs)**3.5d0 &
!                *( 3*C*Fs*(D+ss*Hs)+3*C*Hs*(1+ss*Fs)+4*B*Hs*(D+ss*Hs) &
!                  +12*A*Hs*(D+ss*Hs)**2 ) &
!                -27.0d0/8.0d0*Pi*Hs*s1/sqrt(A)*exp(2.25d0*Hs*ss/A) &
!                *(1.0d0-bberf(1.5d0*s1*sqrt(Hs/A))) &
!                +2.25d0*sqrt(Pi*Hs)

!          dbGds = 15.0d0/8.0d0*sqrt(Pi)*s1/(D+ss*Hs)**3.5d0 &
!                - 105.0d0/16.0d0*sqrt(Pi)*ss*s1*Hs/(D+ss*Hs)**4.5d0

!          dEGds = -daGds/bG + (4.0d0*Pi/3.0d0+aG)/(bG*bG)*dbGds

          srpi=sqrt(Pi)
          DHs=D+Hs*ss
          DHs2=DHs*DHs
          DHs3=DHs2*DHs
          DHs92=DHs2*DHs2*sqrt(DHs)
          A12=sqrt(A)
          dsdDHs=2.0d0*s1*Hs+dHds*ss
          daGds=1.d0/32.d0*srpi*((36.d0*(2.d0*Hs+dHds*s1)/sqrt(Hs)) &
               +(-8.d0*A*dsdDHs*DHs3-105.d0*dsdDHs*E &
               -30.d0*C*dsdDHs*DHs*(1.d0+ss*Fs) &
               +12.d0*DHs2*(-B*dsdDHs+C*s1*(dFds*s1+2.d0*Fs)))/DHs92 &
               -(54.d0*exp(2.25d0*Hs*ss/A)*srpi*s1*(2.d0*Hs+dHds*s1) &
               *(1.d0-bberf(1.5d0*sqrt(Hs/A)*s1)))/A12)
          dbGds=15.d0*srpi*s1*(4.d0*DHs-7.d0*dsdDHs*s1)/(32.d0*DHs92)
          dEGds=(-4.d0*daGds*bG+dbGds*(4.d0*aG+3.d0*Pi))/(4.d0*bG*bG)

          if ( s1 <= EGscut ) dEGds = 2.0d0*EGc2*s1 + 4.0d0*EGc3*ss*s1

! ---

          x = omega/kf
          do j=0,8
             dd(j) = cc(j)*x**j
          end do
          dd(9) = cc(9)*x**2

          A1 = 4.0d0/9.0d0*A
          H1 = ss*Hs
          H2 = ss*Hs + dd(9)

          expint_ha = -expint(1,H2/A1)
          erfc_ha = 1.0d0 - bberf(sqrt(H2/A1))
          Jy(0) = 0.5d0*Pi/sqrt(A1)*exp(H2/A1)*erfc_ha
          Jy(1) =-0.5d0/A1*exp(H2/A1)*expint_ha
          Jy(2) =-0.5d0*Pi/sqrt(A1)/A1 &
               *( erfc_ha*exp(H2/A1) - sqrt(A1/(H2*Pi)) )
          Jy(3) = 0.5d0/A1*( 1.0d0/H2 + exp(H2/A1)*expint_ha/A1 )
          Jy(4) = 0.5d0*Pi/sqrt(A1)/A1**2*( exp(H2/A1)*erfc_ha &
               +sqrt(A1/(Pi*H2))*( 0.5d0*A1/H2 - 1.0d0 ) )
          Jy(5) =-0.5d0/A1**2*( expint_ha*exp(H2/A1)/A1 + (H2-A1)/H2**2 )
          Jy(6) =-0.5d0*Pi/sqrt(A1)/A1**2*( exp(H2/A1)/A1*erfc_ha &
               -1.0d0/(A1*sqrt(H2/A1)*sqrt(Pi)) &
               +0.5d0/H2**1.5d0*sqrt(A1/Pi) &
               -3.0d0/4.0d0*A1/H2**2.5d0*sqrt(A1/Pi) )
          Jy(7) = 0.5d0/A1**3*( exp(H2/A1)*expint_ha/A1 &
               +1.0d0/H2 + (2.0d0*A1*A1-A1*H2)/H2**3 )
          Jy(8) =-0.5d0*Pi/A1**3.5d0*( exp(H2/A1)*erfc_ha/A1 &
                                     - 1.0d0/sqrt(Pi*A1*H2) ) &
                 -0.5d0*sqrt(Pi)/A1**3.5d0*( 0.5d0/H2**1.5d0 &
                 - 0.75d0*A1/H2**2.5d0 + 1.875d0*A1**2/H2**3.5d0 )
          Jy(9) =-0.5d0/A1**3*( exp(H2/A1)*expint_ha/A1**2 &
                               + 1.0d0/(A1*H2) - 1.0d0/H2**2 ) &
                 +1.0d0/A1**2*( 3.0d0*A1/H2**4 - 1.0d0/H2**3 )

          Iy( 1) =-0.5d0*A*log(1.0d0+D/H2)
          Iy( 2) = 0.0d0
          Iy( 3) = A*A1*Jy(1)
          Iy( 4) =-A*dd(1)*Jy(0)
          Iy( 5) =-A*dd(2)*Jy(1)
          Iy( 6) =-A*dd(3)*Jy(2)
          Iy( 7) =-A*dd(4)*Jy(3)
          Iy( 8) =-A*dd(5)*Jy(4)
          Iy( 9) =-A*dd(6)*Jy(5)
          Iy(10) =-A*dd(7)*Jy(6)
          Iy(11) =-A*dd(8)*Jy(7)

          H3 = D + H2

          Igauss(0) = 0.5d0*sqrt(Pi/H3)
          Igauss(1) = 0.5d0/H3
          do j=2,15
             Igauss(j) = (j-1)/(2.0d0*H3)*Igauss(j-2)
          end do

          do j=1,8
             Iy(11+j) = A*dd(j)*Igauss(j-1)
          end do

          FF=1.0d0+ss*Fs
          GG=E+ss*EGs
          Iy(20) = ( B                               )*Igauss( 1)
          Iy(21) = ( B*dd(1)                         )*Igauss( 2)
          Iy(22) = ( B*dd(2) + C*FF                  )*Igauss( 3)
          Iy(23) = ( B*dd(3) + C*FF*dd(1)            )*Igauss( 4)
          Iy(24) = ( B*dd(4) + C*FF*dd(2) + GG       )*Igauss( 5)
          Iy(25) = ( B*dd(5) + C*FF*dd(3) + GG*dd(1) )*Igauss( 6)
          Iy(26) = ( B*dd(6) + C*FF*dd(4) + GG*dd(2) )*Igauss( 7)
          Iy(27) = ( B*dd(7) + C*FF*dd(5) + GG*dd(3) )*Igauss( 8)
          Iy(28) = ( B*dd(8) + C*FF*dd(6) + GG*dd(4) )*Igauss( 9)
          Iy(29) = (           C*FF*dd(7) + GG*dd(5) )*Igauss(10)
          Iy(30) = (           C*FF*dd(8) + GG*dd(6) )*Igauss(11)
          Iy(31) = (                        GG*dd(7) )*Igauss(12)
          Iy(32) = (                        GG*dd(8) )*Igauss(13)

          Fx = -8.0d0/9.0d0*sum( Iy(1:32) )

          Ex = Ex + trho*ex_lda*Fx

! --- potential

          Iy(:) = 0.0d0

! --- dFx/ds

          do j=1,9
             Iy(j) = -A*dd(j-1)*Jy(j)
          end do
          do j=1,9
             Iy(10) = Iy(10) + A*dd(j-1)*Igauss(j)
          end do
          do j=1,9
             Iy(11) = Iy(11) + B*dd(j-1)*Igauss(j+2)
          end do
          do j=1,9
             Iy(12) = Iy(12) + C*(1.0d0+ss*Fs)*dd(j-1)*Igauss(j+4)
          end do
          do j=1,9
             Iy(13) = Iy(13) + (E+ss*EGs)*dd(j-1)*Igauss(j+6)
          end do
          Iy(1:13) = Iy(1:13)*( -2.0d0*s1*Hs - ss*dHds )

          do j=1,9
             Iy(14) = Iy(14) + C*( 2.0d0*s1*Fs + ss*dFds )*dd(j-1)*Igauss(j+2)
          end do
          do j=1,9
             Iy(15) = Iy(15) + (2.0d0*s1*EGs+ss*dEGds)*dd(j-1)*Igauss(j+4)
          end do

          dFxds = -8.0d0/9.0d0*sum(Iy(1:15))

! --- dFx/dn

          dd(1) = -cc(1)*x
          do j=2,8
             dd(j) = ( 2.0d0*cc(9)*cc(j-2) - j*cc(j) )*x**j
          end do
          dd( 9) = 2.0d0*cc(9)*cc(7)*x**9
          dd(10) = 2.0d0*cc(9)*cc(8)*x**10 
          dd(11) = cc(9)*x**2

          do j=1,10
             Iy(j+15) = -A*dd(j)*Jy(j-1)
          end do

          do j=1,10
             Iy(26) = Iy(26) + A*dd(j)*Igauss(j-1)
          end do
          do j=1,10
             Iy(27) = Iy(27) + B*dd(j)*Igauss(j+1)
          end do
          do j=1,10
             Iy(28) = Iy(28) + C*(1.0d0+ss*Fs)*dd(j)*Igauss(j+3)
          end do
          do j=1,10
             Iy(29) = Iy(29) + (E+ss*EGs)*dd(j)*Igauss(j+5)
          end do

          dFxdn = -8.0d0/9.0d0*Pi*Pi/kf**3 * sum(Iy(16:29))

! ---

          vx(i,s) = vx_lda*Fx + trho*ex_lda*dFxdn - 4.0d0/3.0d0*s1*ex_lda*dFxds
          rtmp(i) = -3.0d0/( sqrt(grad%gg(i))*8.0d0*Pi )*dFxds

       end do ! i

       rrrr(m1:m2,1) = rtmp(m1:m2)*grad%gx(m1:m2)
       rrrr(m1:m2,2) = rtmp(m1:m2)*grad%gy(m1:m2)
       rrrr(m1:m2,3) = rtmp(m1:m2)*grad%gz(m1:m2)
       do i=1,3
          call mpi_allgatherv(rrrr(m1,i),ir_grid(myrank_g),mpi_real8 &
               ,rrrr(1,i),ir_grid,id_grid,mpi_real8,comm_grid,ierr)
       end do

       select case( SYStype )
       case default

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1
             i=LLL(i1,i2,i3)
             do m=-Md,Md
                cm=nab(m)*sign(1,m)
                j1=i1+m
                k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
                j1=j1-k1*ML1
                j =LLL(j1,i2,i3)
! The potential vx is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
                if ( m1 <= j .and. j <= m2 ) then
                   vx(j,s) = vx(j,s) + cm*( rrrr(i,1)*q(1,1) &
                                           +rrrr(i,2)*q(2,1) &
                                           +rrrr(i,3)*q(3,1) )
                end if
                j2=i2+m
                k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
                j2=j2-k2*ML2
                j =LLL(i1,j2,i3)
                if ( m1 <= j .and. j <= m2 ) then
                   vx(j,s) = vx(j,s) + cm*( rrrr(i,1)*q(1,2) &
                                           +rrrr(i,2)*q(2,2) &
                                           +rrrr(i,3)*q(3,2) )
                end if
                j3=i3+m
                k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
                j3=j3-k3*ML3
                j =LLL(i1,i2,j3)
                if ( m1 <= j .and. j <= m2 ) then
                   vx(j,s) = vx(j,s) + cm*( rrrr(i,1)*q(1,3) &
                                           +rrrr(i,2)*q(2,3) &
                                           +rrrr(i,3)*q(3,3) )
                end if
             end do ! m
          end do ! i1
          end do ! i2
          end do ! i3

       end select

    end do ! s

! ---

    sb(1)=Ex*rgrid%VolumeElement
    call MPI_ALLREDUCE( sb, rb, 1, MPI_REAL8, MPI_SUM, comm_grid, i )

    Ex = rb(1)/factor

    ene%Ex = Ex
    pot%x%val(:,:) = vx(:,:)

! ---

    deallocate( rrrr )
    deallocate( rtmp )
    deallocate( vx   )
    deallocate( LLL  )
    call destruct_gradient( grad )

  END SUBROUTINE calc_pbe_xsr


END MODULE xc_pbe_xsr_module
