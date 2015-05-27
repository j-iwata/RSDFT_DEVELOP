MODULE xc_ggaak13_module

  use aa_module
  use bb_module
  use bc_module
  use parallel_module
  use fd_module
  use xc_ggapbe96_mol_module


  implicit none

  PRIVATE
  PUBLIC :: init_GGAAK13, calc_GGAAK13

  integer :: Igrid(2,0:3)
  integer :: ML_0, ML_1, MSP_0, MSP_1, MSP, comm
  real(8) :: dV
  logical :: flag_init = .true.

  real(8),allocatable :: nab(:)

  real(8),allocatable :: gx(:),gy(:),gz(:)
  real(8),parameter :: zero_density = 1.d-10

  real(8) :: b(3,3)

  integer,allocatable :: LLL2(:,:,:)
  integer :: Md,ML,ML1,ML2,ML3
  real(8) :: Hgrid(3)
  integer :: SYStype

CONTAINS


  SUBROUTINE init_GGAAK13 &
       ( Igrid_in,MSP_0_in,MSP_1_in,MSP_in,comm_in,dV_in &
        ,Md_in,Hgrid_in,Ngrid_in,SYStype_in )
    implicit none
    integer,intent(IN) :: Igrid_in(2,0:3),MSP_0_in,MSP_1_in,MSP_in,comm_in
    real(8),intent(IN) :: dV_in, Hgrid_in(3)
    integer,intent(IN) :: Md_in, Ngrid_in(0:3), SYStype_in
    real(8) :: aaL(3), Pi

    if ( .not.flag_init ) return

    Igrid(:,:) = Igrid_in(:,:)
    ML_0       = Igrid(1,0)
    ML_1       = Igrid(2,0)
    MSP_0      = MSP_0_in
    MSP_1      = MSP_1_in
    MSP        = MSP_in
    comm       = comm_in
    dV         = dV_in
    flag_init  = .false.
    Md         = Md_in
    Hgrid(:)   = Hgrid_in(:)
    ML  = Ngrid_in(0)
    ML1 = Ngrid_in(1)
    ML2 = Ngrid_in(2)
    ML3 = Ngrid_in(3)

    SYStype = SYStype_in

    aaL(1) = sqrt( sum(aa(1:3,1)**2) )
    aaL(2) = sqrt( sum(aa(1:3,2)**2) )
    aaL(3) = sqrt( sum(aa(1:3,3)**2) )
    Pi     = acos(-1.0d0)
    b(:,:) = 0.0d0
    b(1:3,1)=aaL(1)*bb(1:3,1)/(2.0d0*Pi)/Hgrid(1)
    b(1:3,2)=aaL(2)*bb(1:3,2)/(2.0d0*Pi)/Hgrid(2)
    b(1:3,3)=aaL(3)*bb(1:3,3)/(2.0d0*Pi)/Hgrid(3)

    allocate( nab(-Md:Md) ) ; nab=0.d0
    call get_coef_nabla_fd( Md, nab )

    flag_init = .false.

  END SUBROUTINE init_GGAAK13


  SUBROUTINE calc_GGAAK13 &
       ( rho, Exc_out, Vxc_out, Ex_out, Ec_out, Vx_out, Vc_out )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: Exc_out
    real(8),optional,intent(OUT) :: Vxc_out(ML_0:,MSP_0:)
    real(8),optional,intent(OUT) :: Ex_out, Ec_out
    real(8),optional,intent(OUT) :: Vx_out(ML_0:,MSP_0:),Vc_out(ML_0:,MSP_0:)
    integer :: i1,i2,i3,i,irank,j1,j2,j3,l1,l2,l3,m1,m2,m3,ierr
    real(8) :: s0(2),s1(2), Ex_part, Ec_part
    real(8),allocatable :: vxc_tmp(:,:)

    if ( flag_init ) then
       write(*,*) "Call INIT_GGAAK13 first"
       stop "stop@calc_GGAAK13(in xc_ggaak13_module)"
    end if

    Ex_part = 0.0d0
    Ec_part = 0.0d0

    Exc_out = 0.0d0
    if ( present(Ex_out)  ) Ex_out =0.0d0
    if ( present(Ec_out)  ) Ec_out =0.0d0
    if ( present(Vxc_out) ) Vxc_out=0.0d0
    if ( present(Vx_out)  ) Vx_out =0.0d0
    if ( present(Vc_out)  ) Vc_out =0.0d0

    select case( SYStype )
    case default

    allocate( LLL2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL2=0
    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       l1=pinfo_grid(1,irank) ; m1=pinfo_grid(2,irank)+l1-1
       l2=pinfo_grid(3,irank) ; m2=pinfo_grid(4,irank)+l2-1
       l3=pinfo_grid(5,irank) ; m3=pinfo_grid(6,irank)+l3-1
       do j3=l3,m3
       do j2=l2,m2
       do j1=l1,m1
          i=i+1
          LLL2(j1,j2,j3)=i
       end do
       end do
       end do
    end do
    end do
    end do

    call construct_gradient( rho )

    case( 1,2 )

       m1 = ( ML1-1 )/2 + Md
       m2 = ( ML2-1 )/2 + Md
       m3 = ( ML3-1 )/2 + Md

       allocate( LLL2(-m1:m1,-m2:m2,-m3:m3) ) ; LLL2=0

       call get_LLL_mol( m1,m2,m3,LLL2 )

       allocate( gx(ML_0:ML_1) ) ; gx=0.0d0
       allocate( gy(ML_0:ML_1) ) ; gy=0.0d0
       allocate( gz(ML_0:ML_1) ) ; gz=0.0d0

       call construct_gradient_mol( Md,nab,rho,gx,gy,gz )

    end select

! --

    allocate( vxc_tmp(ML_0:ML_1,1:MSP) ) ; vxc_tmp=0.0d0

! -- Exchange --

    call calc_GGAAK13_x( rho, vxc_tmp, Ex_part )

    if ( present(Vx_out) ) then
       Vx_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

    if ( present(Vxc_out) ) then
       Vxc_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

! -- Correlation --

! call calc_GGAPBE96_c2( rho, vxc_tmp, Ec_part )
    call calc_GGAAK13_c92( rho, vxc_tmp, Ec_part )

    if ( present(Vc_out) ) then
       Vc_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

    if ( present(Vxc_out) ) then
       Vxc_out(ML_0:ML_1,MSP_0:MSP_1) = Vxc_out(ML_0:ML_1,MSP_0:MSP_1) &
                                      + vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

! --

    s0(1) = Ex_part*dV
    s0(2) = Ec_part*dV
    call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,comm,ierr)

    s0(1) = s1(1)
    call mpi_allreduce(s0,s1,1,mpi_real8,mpi_sum,comm_spin,ierr)

    Exc_out = s1(1) + s1(2)

    if ( present(Ex_out) ) Ex_out=s1(1)
    if ( present(Ec_out) ) Ec_out=s1(2)

! --

    deallocate( vxc_tmp )

    deallocate( gz, gy, gx )

    deallocate( LLL2 )

  END SUBROUTINE calc_GGAAK13


  SUBROUTINE construct_gradient(rho)
    implicit none
    real(8),intent(IN) :: rho(ML_0:,:)
    integer :: i,i1,i2,i3,s,m
    real(8) :: g1,g2,g3
    allocate( gx(ML_0:ML_1) ) ; gx=0.0d0
    allocate( gy(ML_0:ML_1) ) ; gy=0.0d0
    allocate( gz(ML_0:ML_1) ) ; gz=0.0d0
    www(:,:,:,:)=0.0d0
    do s=1,MSP
       i=ML_0-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho(i,s)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       g1=0.0d0
       g2=0.0d0
       g3=0.0d0
       do m=1,Md
          g1 = g1 - nab(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          g2 = g2 - nab(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          g3 = g3 - nab(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
       end do
       i=i+1
       gx(i) = b(1,1)*g1 + b(1,2)*g2 + b(1,3)*g3
       gy(i) = b(2,1)*g1 + b(2,2)*g2 + b(2,3)*g3
       gz(i) = b(3,1)*g1 + b(3,2)*g2 + b(3,3)*g3
    end do
    end do
    end do
  END SUBROUTINE construct_gradient


  SUBROUTINE calc_GGAAK13_x( rho, vex, Ex )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vex(ML_0:,:)
    real(8),intent(OUT) :: Ex
    !real(8),parameter :: mu=0.21951d0, Kp=0.804d0
    integer :: i,ispin,ierr
    real(8) :: rhoa,rhob,trho,cm
    real(8),allocatable :: rrrr(:,:),rtmp(:)
    real(8) :: kf, vx_lda, ex_lda, Fx, Pi, g2, factor, s, muge, B1, dFx_ds ! s, muge, B1 and dFx_ds are added
    integer :: i1,i2,i3,j,j1,j2,j3,k1,k2,k3,m,m1,m2,m3

    Pi = acos(-1.0d0)
    muge = 10.d0/81.d0
    B1 = (3.d0/5.d0)*muge + 8.d0*Pi/15.d0

    factor = 1.0d0
    if ( MSP == 2 ) factor = 2.0d0

    vex = 0.0d0
    Ex  = 0.0d0

    allocate( rtmp(ML_0:ML_1) ) ; rtmp=0.0d0
    allocate( rrrr(ML,3)      ) ; rrrr=0.0d0

    do ispin=MSP_0,MSP_1

       do i=ML_0,ML_1

          trho = factor*rho(i,ispin)

          if ( trho <= zero_density ) cycle

          kf = (3.0d0*Pi*Pi*trho)**(1.0d0/3.0d0)

          ex_lda = -3.0d0/(4.0d0*Pi)*kf
          vx_lda = -1.0d0/Pi*kf

          g2 = gx(i)*gx(i) + gy(i)*gy(i) + gz(i)*gz(i)
          ! calc s
          s = dsqrt(g2)/(2.d0*kf*trho)
          !Fx = 1.0d0 + Kp - 4.0d0*Kp*Kp*(trho*kf)**2/( 4.0d0*Kp*(trho*kf)**2 + mu*g2 )
          Fx = 1.d0 + B1*s*dlog(1.d0+s) + (muge-B1)*s*dlog(1.d0+dlog(1.d0+s))
          Ex = Ex + trho*ex_lda*Fx
          dFx_ds =  B1*( dlog(1.d0+s)+s/(1.d0+s) ) + (muge-B1)*( dlog(1.d0+dlog(1.d0+s))+s/(1.d0+s)/(1.d0+dlog(1.d0+s)) )
          vex(i,ispin) = vex(i,ispin) + Fx*vx_lda + 3.d0*s*kf/(4.0d0*Pi)*dFx_ds
          !rtmp(i) = -18.0d0*Pi*Kp*Kp*mu*trho**4/( 4.0d0*Kp*(trho*kf)**2 + mu*g2 )**2
          rtmp(i) = -3.d0/(8.d0*Pi*dsqrt(g2))*dFx_ds
       end do ! i

       rrrr(ML_0:ML_1,1) = rtmp(ML_0:ML_1)*gx(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,1),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       rrrr(ML_0:ML_1,2) = rtmp(ML_0:ML_1)*gy(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,2),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       rrrr(ML_0:ML_1,3)=rtmp(ML_0:ML_1)*gz(ML_0:ML_1)
       call mpi_allgatherv(rrrr(ML_0,3),ir_grid(myrank_g),mpi_real8 &
            ,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       select case( SYStype )
       case default

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1
          i=LLL2(i1,i2,i3)
          do m=-Md,Md
             cm=nab(m)*sign(1,m)
             j1=i1+m
             k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
             j1=j1-k1*ML1
             j =LLL2(j1,i2,i3)
! The potential vex is calculated at j-th grid point rather than i-th.
! This is because non-transposed nabla matrix Dij is used (See XC.doc).
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,1) &
                                                  +rrrr(i,2)*b(2,1) &
                                                  +rrrr(i,3)*b(3,1) )
             end if
             j2=i2+m
             k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
             j2=j2-k2*ML2
             j =LLL2(i1,j2,i3)
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,2) &
                                                  +rrrr(i,2)*b(2,2) &
                                                  +rrrr(i,3)*b(3,2) )
             end if
             j3=i3+m
             k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
             j3=j3-k3*ML3
             j =LLL2(i1,i2,j3)
             if ( ML_0 <= j .and. j <= ML_1 ) then
                vex(j,ispin) = vex(j,ispin) + cm*( rrrr(i,1)*b(1,3) &
                                                  +rrrr(i,2)*b(2,3) &
                                                  +rrrr(i,3)*b(3,3) )
             end if
          end do ! m
       end do ! i1
       end do ! i2
       end do ! i3

       case( 1,2 )

          m1 = (ML1-1)/2 + Md
          m2 = (ML2-1)/2 + Md
          m3 = (ML3-1)/2 + Md

          call calc_ve_mol &
               ( ML_0, ML_1, m1,m2,m3, Md, nab, vex(:,ispin), LLL2, rrrr )

       end select

    end do ! ispin

    Ex = Ex/factor

    deallocate( rrrr )
    deallocate( rtmp )

  END SUBROUTINE calc_GGAAK13_x
  

  SUBROUTINE calc_GGAAK13_c92( rho, vco, Ec )
  !LDAPW92 correlation
    implicit none

    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vco(ML_0:,:),Ec

    real(8) :: trho,zeta,d0,d1,d2,d3,d4,d5,d6,d7,d8,fx,fpx,rs,epsilon_c
    real(8) :: de_dr,de_dz,dec_dr(2),dac_dr,Q0(3),Q1(3),Qp1(3)
    real(8),parameter :: A(1:3)=(/0.031091d0,0.015545d0,0.016887d0/)
    real(8),parameter :: alpha1(1:3)=(/0.21370d0,0.20548d0,0.11125d0/)
    real(8),parameter :: beta1(1:3)=(/7.5957d0,14.1189d0,10.357d0/)
    real(8),parameter :: beta2(1:3)=(/3.5876d0,6.1977d0,3.6231d0/)
    real(8),parameter :: beta3(1:3)=(/1.6382d0,3.3662d0,0.88026d0/)
    real(8),parameter :: beta4(1:3)=(/0.49294d0,0.62517d0,0.49671d0/)
    real(8),parameter :: fpp0=1.709921d0
    real(8) :: sgn(2),ecd(2),ac
    integer :: i
    
    d0 = 0.5d0*MSP
    d1 = 4.d0/3.d0
    d2 = 1.d0/( 2.d0*(2.d0**(1.d0/3.d0)-1.d0) )
    d3 = 3.d0/(4.d0*acos(-1.d0))
    d4 = 1.d0/3.d0
    d5 = 1.d0/(4.d0*acos(-1.d0))
    d6 = 1.d0/2.d0
    d7 = 3.d0/2.d0
    d8 = 9.d0*acos(-1.d0)/4.d0

    sgn(1) = 1.d0
    sgn(2) = -1.d0

    Ec = 0.0d0

    do i=ML_0,ML_1
       trho = d0*( rho(i,1)+rho(i,MSP) )

       if ( trho <= 0.0d0 ) cycle

       zeta = ( rho(i,1)-rho(i,MSP) )/trho

       fx = ( abs(1.d0+zeta)**d1 + abs(1.d0-zeta)**d1 - 2.d0 )*d2
       fpx = d1*( abs(1.d0+zeta)**d4 - abs(1.d0-zeta)**d4 )*d2
       rs = (d3/trho)**d4

       ecd(1) = -2.d0*A(1)*(1+alpha1(1)*rs)*dlog(1.d0+1.d0/(2.d0*A(1) &
        *(beta1(1)*rs**d6+beta2(1)*rs+beta3(1)*rs**d7+beta4(1)*rs**2.d0)))
       ecd(2) = -2.d0*A(2)*(1+alpha1(2)*rs)*dlog(1.d0+1.d0/(2.d0*A(2) &
        *(beta1(2)*rs**d6+beta2(2)*rs+beta3(2)*rs**d7+beta4(2)*rs**2.d0)))
       
       ac = 2.d0*A(3)*(1+alpha1(3)*rs)*dlog(1.d0+1.d0/(2.d0*A(3) &
        *(beta1(3)*rs**d6+beta2(3)*rs+beta3(3)*rs**d7+beta4(3)*rs**2.d0)))
       
       epsilon_c = ecd(1)+ac*(fx/fpp0)*(1.d0-zeta**4)+(ecd(2)-ecd(1))*fx*zeta**4
       Ec = Ec + trho*epsilon_c

       Q0 = -2.d0*A*(1.d0+alpha1*rs)
       Q1 = 2.d0*A*(beta1*rs**d6+beta2*rs+beta3*rs**d7+beta4*rs**2.d0)
       Qp1 = A*(beta1*rs**(-d6)+2.d0*beta2+3.d0*beta3*rs**d6+4.d0*beta4*rs)
       dec_dr(1) = -2.d0*A(1)*alpha1(1)*dlog(1.d0+1.d0/Q1(1)) &
                   -Q0(1)*Qp1(1)/(Q1(1)**2+Q1(1))
       dec_dr(2) = -2.d0*A(2)*alpha1(2)*dlog(1.d0+1.d0/Q1(2)) &
                   -Q0(2)*Qp1(2)/(Q1(2)**2+Q1(2))
       dac_dr = 2.d0*A(3)*alpha1(3)*dlog(1.d0+1.d0/Q1(3)) &
                   -Q0(3)*Qp1(3)/(Q1(3)**2+Q1(3))
       de_dr = dec_dr(1)*(1.d0-fx*zeta**4)+dec_dr(2)*fx*zeta**4 &
                +dac_dr*(fx/fpp0)*(1.d0-zeta**4.d0)
       de_dz = 4.d0*zeta**3*fx*(ecd(2)-ecd(1)-ac/fpp0) &
                +fpx*( ecd(2)*zeta**4-ecd(1)*zeta**4+(1-zeta**4)*ac/fpp0 )

       vco(i,MSP_0) = epsilon_c - d4*rs*de_dr - (zeta-sgn(MSP_0))*de_dz

       vco(i,MSP_1) = epsilon_c - d4*rs*de_dr - (zeta-sgn(MSP_1))*de_dz

    end do ! i

    return
  END SUBROUTINE calc_GGAAK13_c92

END MODULE xc_ggaak13_module
