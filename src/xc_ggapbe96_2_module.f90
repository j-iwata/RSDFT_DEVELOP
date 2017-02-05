MODULE xc_ggapbe96_2_module

  use gradient_module
  use grid_module, only: grid, get_map_3d_to_1d_grid
  use xc_variables, only: xcpot, xcene
  use fd_module, only: fd, construct_nabla_fd
  use lattice_module, only: lattice, get_aa_lattice, get_reciprocal_lattice
  use parallel_module
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: calc_GGAPBE96_2

  integer,parameter :: DP=kind(0.0d0)
!#ifdef _NO_QPRECISION_
  integer,parameter :: QP=kind(0.0d0)
!#else
!  integer,parameter :: QP=kind(0.0q0)
!#endif

  real(QP),parameter :: zero_density = 1.e-10_QP
  real(QP),allocatable :: nab(:)
  real(QP),allocatable :: vx(:,:),vc(:,:)
  real(QP) :: Ex,Ec
  real(QP) :: b(3,3)
  integer,allocatable :: LLL(:,:,:)
  integer :: Md, ML1,ML2,ML3
  integer :: SYStype=0

  real(8) :: mu=0.21951d0
  real(8) :: Kp=0.804d0
  real(8) :: beta=0.066725d0

CONTAINS


  SUBROUTINE calc_GGAPBE96_2( rgrid, rho, ene, pot, mu_in, Kp_in )

    implicit none
    type( grid ),intent(IN) :: rgrid
    type( GSArray ),intent(IN) :: rho
    type( xcene ) :: ene
    type( xcpot ) :: pot
    real(8),optional,intent(IN) :: mu_in, Kp_in
    type(gradient16) :: grad16
    type(fd) :: nabla
    type(lattice) :: aa, bb
    integer :: m1,m2,n1,n2,i,i1,i2,i3
    real(QP) :: Pi
    real(DP) :: sb(2),rb(2)

! ---

    if ( present(mu_in) ) mu=mu_in
    if ( present(Kp_in) ) Kp=Kp_in
    beta = mu*3.0d0/acos(-1.0d0)**2

    call construct_gradient16( rgrid, rho, grad16 )

    call construct_nabla_fd( nabla )

    Md = nabla%md

    if ( .not.allocated(nab) ) allocate( nab(-Md:Md) )
    nab(:) = nabla%coef(:)

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )

    Pi=acos(-1.0_QP)
    b(1:3,1)=aa%Length(1)*bb%LatticeVector(1:3,1)/( 2*Pi*rgrid%spacing(1) )
    b(1:3,2)=aa%Length(2)*bb%LatticeVector(1:3,2)/( 2*Pi*rgrid%spacing(2) )
    b(1:3,3)=aa%Length(3)*bb%LatticeVector(1:3,3)/( 2*Pi*rgrid%spacing(3) )

    ML1 = rgrid%g3%x%size_global
    ML2 = rgrid%g3%y%size_global
    ML3 = rgrid%g3%z%size_global

    call get_map_3d_to_1d_grid( rgrid, LLL )

! ---

    m1 = pot%xc%g_range%head
    m2 = pot%xc%g_range%tail
    n1 = rho%s_range%head
    n2 = rho%s_range%tail

    allocate( vx(m1:m2,n1:n2) ) ; vx=0.0_QP

    call calc_PBE_x( rgrid, rho, grad16 )

    allocate( vc(m1:m2,n1:n2) ) ; vc=0.0_QP

    call calc_PBE_c( rgrid, rho, grad16 )

! ---

    sb(1)=Ex*rgrid%VolumeElement
    sb(2)=Ec*rgrid%VolumeElement
    call MPI_ALLREDUCE( sb, rb, 2, MPI_REAL8, MPI_SUM, comm_grid, i )

    ene%Ex  = rb(1)
    ene%Ec  = rb(2)
    ene%Exc = rb(1)+rb(2)

    n1 = pot%xc%s_range%head
    n2 = pot%xc%s_range%tail

    pot%xc%val(:,:)  = vx(:,n1:n2) + vc(:,n1:n2)
    if ( allocated(pot%x%val) ) pot%x%val(:,:) = vx(:,n1:n2)
    if ( allocated(pot%c%val) ) pot%c%val(:,:) = vc(:,n1:n2)

! ---

    deallocate( vc )
    deallocate( vx )
    deallocate( LLL )
    call destruct_gradient16( grad16 )

  END SUBROUTINE calc_GGAPBE96_2


  SUBROUTINE calc_PBE_x( rgrid, rho, grad16 )

    implicit none
    type( grid ) :: rgrid
    type( GSArray ) :: rho
    type( gradient16 ) :: grad16
!   real(8),parameter :: mu=0.21951d0, Kp=0.804d0
    integer :: i,ispin,m,ierr
    real(QP) :: trho,cm
    real(QP),allocatable :: rrrr(:,:),rtmp(:)
    real(QP) :: kf, vx_lda, ex_lda, Fx, Pi, g2, factor
    real(QP) :: onethr,const1,const2
    integer :: i1,i2,i3,j,j1,j2,j3,k1,k2,k3
    integer :: mm,m1,m2
    real(8),allocatable :: f(:), Gf(:)

    Pi = acos(-1.0_QP)

    factor = 1.0_QP
    if ( rho%s_range%size_global == 2 ) factor = 2.0_QP

    onethr = 1.0_QP/3.0_QP
    const1 = 3.0_QP*Pi*Pi
    const2 = 3.0_QP/(4.0_QP*Pi)

    m1 = rho%g_range%head
    m2 = rho%g_range%tail
    allocate( rtmp(m1:m2) ) ; rtmp=0.0_QP

    Ex = 0.0_QP

    do ispin=rho%s_range%head,rho%s_range%tail

       do i=rho%g_range%head,rho%g_range%tail

          trho = factor*rho%val(i,ispin)

          if ( trho <= zero_density ) cycle

          kf = (const1*trho)**onethr

          ex_lda = -const2*kf
          vx_lda = -kf/Pi

          g2 = grad16%gg(i)

          Fx = 1.0_QP + Kp - 4.0_QP*Kp*Kp*(trho*kf)**2 &
                            /( 4.0_QP*Kp*(trho*kf)**2 + mu*g2 )

          Ex = Ex + trho*ex_lda*Fx

          vx(i,ispin) = vx(i,ispin) &
               + Fx*vx_lda + ( 24.0_QP*Pi*Kp*Kp*mu*trho**3*g2 ) &
                            /( 4.0_QP*Kp*(trho*kf)**2 + mu*g2 )**2

          rtmp(i) = -18.0_QP*Pi*Kp*Kp*mu*trho**4 &
                    /( 4.0_QP*Kp*(trho*kf)**2 + mu*g2 )**2

       end do ! i

       allocate( rrrr(m1:m2,3)  ) ; rrrr=0.0_QP
       rrrr(m1:m2,1) = rtmp(m1:m2)*grad16%gx(m1:m2)
       rrrr(m1:m2,2) = rtmp(m1:m2)*grad16%gy(m1:m2)
       rrrr(m1:m2,3) = rtmp(m1:m2)*grad16%gz(m1:m2)

       allocate( f(m1:m2)  ) ; f=0.0d0
       allocate( Gf(m1:m2) ) ; Gf=0.0d0

       do j=1,3
          f(:) = rrrr(m1:m2,1)*b(1,j) + rrrr(m1:m2,2)*b(2,j) + rrrr(m1:m2,3)*b(3,j)
          call calc_abc_gradient( j, rgrid, f, Gf )
          vx(:,ispin) = vx(:,ispin) - Gf(:)
       end do

       deallocate( Gf )
       deallocate( f  )
       deallocate( rrrr )

    end do ! ispin

    Ex = Ex/factor

    deallocate( rtmp )

  END SUBROUTINE calc_PBE_x


  SUBROUTINE calc_PBE_c( rgrid, rho, grad16 )

    implicit none
    type( grid ) :: rgrid
    type( GSArray ) :: rho
    type( gradient16 ) :: grad16
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
!   real(8),parameter :: C1=2.14611945579107d0,C2=0.031091d0
    real(QP) :: C1,C2
    real(QP) :: const1, const2, factor
    integer :: i,j,ispin,m,i1,i2,i3,j1,j2,j3,k1,k2,k3,ierr
    integer :: m1,m2
    real(QP) :: trho, rhoa, rhob
    real(QP) :: kf, rs, ec_U, ec_P, ec_lda, phi, g2
    real(QP) :: dac_dn, dfz_dz, deU_dn, deP_dn, H, A, T, alpc, fz
    real(QP) :: dec_dz, dphi_dz, dH_dA, dH_dT, tmp, dH_dphi, Ai
    real(QP) :: dA_dn, dec_dn, Hs, zeta, dT_dphi
    real(QP) :: dz_dn(2), cm, rssq
    real(QP) :: dac_drs, deP_drs, deU_drs, drs_dn
    real(QP),allocatable :: rrrr(:,:), rtmp(:)
    real(QP) :: Pi, one, two, fouthr, onethr, zero
    real(QP) :: sevthr, twothr, thrtwo, ThrFouPi
    real(8),allocatable :: f(:),Gf(:)

    const1 = 2.0_QP**(4.0_QP/3.0_QP)-2.0_QP
    const2 = 9.0_QP*(2.0_QP**(1.0_QP/3.0_QP)-1.0_QP)/4.0_QP

    factor = 1.0_QP
    if ( rho%s_range%size_global == 1 ) factor = 0.5_QP
 
    Pi       = acos(-1.0_QP)
    zero     = 0.0_QP
    one      = 1.0_QP
    two      = 2.0_QP
    fouthr   = 4.0_QP/3.0_QP
    onethr   = 1.0_QP/3.0_QP
    ThrFouPi = 3.0_QP/(4.0_QP*Pi)
    thrtwo   = 3.0_QP/2.0_QP
    twothr   = 2.0_QP/3.0_QP
    sevthr   = 7.0_QP/3.0_QP

    C2 = ( 1.0_QP-log(2.0_QP) )/Pi**2 ! "gamma" in PBE paper
    C1 = beta/C2

    Ec = 0.0_QP

    m1 = rho%g_range%head
    m2 = rho%g_range%tail
    allocate( rtmp(m1:m2) ) ; rtmp=0.0_QP

    do i=rho%g_range%head,rho%g_range%tail

       rhoa = rho%val(i,1)*factor
       rhob = rho%val(i,rho%s_range%size_global)*factor
       trho = rhoa + rhob

       if ( trho <= zero_density ) cycle

       zeta = ( rhoa - rhob )/trho
       if ( abs(zeta) > one ) zeta=sign(one,zeta)

       fz = ( (one+zeta)**fouthr + (one-zeta)**fouthr - two )*const1

       kf = ( 3.0_QP*Pi*Pi*trho )**onethr

       rs = ( ThrFouPi/trho )**onethr

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

       ec_lda = ec_U - alpc*fz*const2*(one-zeta**4) &
            + ( ec_P - ec_U )*fz*zeta**4

       phi = 0.5_QP*( (one+zeta)**twothr + (one-zeta)**twothr )

       g2 = grad16%gg(i)

       T = Pi*g2/( 16.0_QP*phi**2*kf*trho**2 )

       A = C1/( exp( -ec_lda/(C2*phi**3) ) - one )

       Hs = C2 * phi**3 * log( one + C1*T*(one+A*T)/(one+A*T+(A*T)**2) )

       Ec = Ec + trho*( ec_lda + Hs )

       drs_dn = -4.0_QP*Pi*rs**4/9.0_QP

       tmp = bt10*rssq + bt20*rs + bt30*rs*rssq + bt40*rs*rs
       deU_drs = alp10*ec_U/( one + alp10*rs ) &
            +A00*( one + alp10*rs )/rssq &
            *( bt10 + two*bt20*rssq + 3.0_QP*bt30*rs + 4.0_QP*bt40*rs*rssq ) &
            /( two*A00*tmp*tmp+tmp )

       tmp = bt11*rssq + bt21*rs + bt31*rs*rssq + bt41*rs*rs
       deP_drs = alp11*ec_P/( one + alp11*rs ) &
            +A01*( one + alp11*rs )/rssq &
            *( bt11 + two*bt21*rssq + 3.0_QP*bt31*rs + 4.0_QP*bt41*rs*rssq ) &
            /( two*A01*tmp*tmp+tmp )

       tmp = bt12*rssq + bt22*rs + bt32*rs*rssq + bt42*rs*rs
       dac_drs = alp12*alpc/( one + alp12*rs ) &
            +A02*( one + alp12*rs )/rssq &
            *( bt12 + two*bt22*rssq + 3.0_QP*bt32*rs + 4.0_QP*bt42*rs*rssq ) &
            /( two*A02*tmp*tmp+tmp )

       deU_dn = deU_drs * drs_dn
       deP_dn = deP_drs * drs_dn
       dac_dn = dac_drs * drs_dn


       dfz_dz = fouthr*( (one+zeta)**onethr - (one-zeta)**onethr )*const1

       dec_dz = -alpc*dfz_dz*const2*(one-zeta**4) &
               + 4.0_QP*alpc*fz*const2*zeta**3 &
               +(ec_P-ec_U)*( dfz_dz*zeta**4 + fz*4.0_QP*zeta**3 )

       if ( abs(zeta) == one ) then
          dphi_dz = zero
       else
          dphi_dz = ( (one+zeta)**(-onethr)-(one-zeta)**(-onethr) )/3.0_QP
       end if

       tmp = one + A*T + (A*T)**2

       dH_dA = -phi**3*C1*C2*A*T**3*(two+A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dT =  phi**3*C1*C2*(one+two*A*T)/(tmp**2+C1*T*(one+A*T)*tmp)

       dH_dphi = 3.0_QP*Hs/phi

       dz_dn(1)   = 2.0_QP*rhob/trho**2
       dz_dn(rho%s_range%size_global) =-2.0_QP*rhoa/trho**2

       do ispin=rho%s_range%head,rho%s_range%tail

          dec_dn = deU_dn - dac_dn*fz*const2*(one-zeta**4) &
               +(deP_dn-deU_dn)*fz*zeta**4 + dec_dz*dz_dn(ispin)

          dA_dn = A*(C1+A)/(C1*C2*phi**3) &
               *( dec_dn - 3.0_QP*ec_lda/phi*dphi_dz*dz_dn(ispin) )

          vc(i,ispin) = vc(i,ispin) + ec_lda + Hs + trho*dec_dn &
               + trho*dH_dA*dA_dn &
               + dH_dT*(-sevthr*T - trho*two*T/phi*dphi_dz*dz_dn(ispin) ) &
               + trho*dH_dphi*dphi_dz*dz_dn(ispin)

       end do ! ispin

       rtmp(i) = dH_dT*Pi/(8.0_QP*kf*trho*phi**2)

    end do ! i

    allocate( rrrr(m1:m2,3)  ) ; rrrr=0.0_QP
    rrrr(m1:m2,1) = rtmp(m1:m2)*grad16%gx(m1:m2)
    rrrr(m1:m2,2) = rtmp(m1:m2)*grad16%gy(m1:m2)
    rrrr(m1:m2,3) = rtmp(m1:m2)*grad16%gz(m1:m2)

    allocate( f(m1:m2)  ) ; f=0.0d0
    allocate( Gf(m1:m2) ) ; Gf=0.0d0

    do j=1,3
       f(:) = rrrr(m1:m2,1)*b(1,j) + rrrr(m1:m2,2)*b(2,j) + rrrr(m1:m2,3)*b(3,j)
       call calc_abc_gradient( j, rgrid, f, Gf )
       do ispin=rho%s_range%head,rho%s_range%tail
          vc(:,ispin) = vc(:,ispin) - Gf(:)
       end do
    end do

    deallocate( Gf )
    deallocate( f  )
    deallocate( rrrr )
    deallocate( rtmp )

    return
  END SUBROUTINE calc_PBE_c


END MODULE xc_ggapbe96_2_module
