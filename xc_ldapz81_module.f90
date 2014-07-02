MODULE xc_ldapz81_module

  implicit none

  PRIVATE
  PUBLIC :: init_LDAPZ81, calc_LDAPZ81

  integer :: ML_0, ML_1, MSP_0, MSP_1, MSP, comm
  real(8) :: dV
  logical :: flag_init = .true.

CONTAINS


  SUBROUTINE init_LDAPZ81 &
       (ML_0_in,ML_1_in,MSP_0_in,MSP_1_in,MSP_in,comm_in,dV_in)
    implicit none
    integer,intent(IN) :: ML_0_in,ML_1_in,MSP_0_in,MSP_1_in,MSP_in,comm_in
    real(8),intent(IN) :: dV_in
    ML_0      = ML_0_in
    ML_1      = ML_1_in
    MSP_0     = MSP_0_in
    MSP_1     = MSP_1_in
    MSP       = MSP_in
    comm      = comm_in
    dV        = dV_in
    flag_init = .false.
  END SUBROUTINE init_LDAPZ81


  SUBROUTINE calc_LDAPZ81 &
       ( rho, Exc_out, Vxc_out, Ex_out, Ec_out, Vx_out, Vc_out )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: Exc_out
    real(8),optional,intent(OUT) :: Vxc_out(ML_0:,MSP_0:)
    real(8),optional,intent(OUT) :: Ex_out, Ec_out
    real(8),optional,intent(OUT) :: Vx_out(ML_0:,MSP_0:),Vc_out(ML_0:,MSP_0:)
    real(8),allocatable :: vxc_tmp(:,:)
    real(8) :: Ex_part, Ec_part, s0(2),s1(2)
    integer :: ierr
    include 'mpif.h'

    if ( flag_init ) then
       write(*,*) "Call INIT_LDAPZ81 first"
       stop "stop@calc_LDAPZ81(in xc_ldapz81_module)"
    end if

    allocate( vxc_tmp(ML_0:ML_1,1:MSP) ) ; vxc_tmp=0.0d0

! -- Exchange --

    call calc_LDAPZ81_x( rho, vxc_tmp, Ex_part )

    if ( present(Vx_out) ) then
       Vx_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

    if ( present(Vxc_out) ) then
       Vxc_out(ML_0:ML_1,MSP_0:MSP_1) = vxc_tmp(ML_0:ML_1,MSP_0:MSP_1)
    end if

! -- Correlation --

    call calc_LDAPZ81_c( rho, vxc_tmp, Ec_part )

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

    Exc_out = s1(1) + s1(2)

    if ( present(Ex_out) ) Ex_out=s1(1)
    if ( present(Ec_out) ) Ec_out=s1(2)

! --

    deallocate( vxc_tmp )

  END SUBROUTINE calc_LDAPZ81


  SUBROUTINE calc_ldapz81_x( rho, vex, Ex )
    implicit none
    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vex(ML_0:,:), Ex
    real(8) :: onetwo,onethr,thrfou,fouthr,thrPi
    real(8) :: cnst,exd
    integer :: i,s

    onethr = 1.0d0/3.0d0
    thrfou = 3.0d0/4.0d0
    fouthr = 4.0d0/3.0d0
    thrPi  = 3.0d0/acos(-1.0d0)

    vex = 0.0d0
    Ex  = 0.0d0

    cnst = 1.0d0
    if ( MSP == 2 ) cnst = (2.0d0)**onethr

    do s=1,MSP

       do i=ML_0,ML_1

          exd = -cnst*thrfou*( thrPi*rho(i,s) )**onethr
          Ex = Ex + rho(i,s)*exd
          vex(i,s) = fouthr*exd

       end do

    end do

  END SUBROUTINE calc_ldapz81_x


  SUBROUTINE calc_ldapz81_c( rho, vco, Ec )
    implicit none

    real(8),intent(IN)  :: rho(ML_0:,:)
    real(8),intent(OUT) :: vco(ML_0:,:),Ec

    real(8),parameter :: gam(1:2)=(/-0.1423d0,-0.0843d0/)
    real(8),parameter :: bet1(1:2)=(/1.0529d0,1.3981d0/)
    real(8),parameter :: bet2(1:2)=(/0.3334d0,0.2611d0/)
    real(8),parameter :: A(1:2)=(/0.0311d0,0.01555d0/)
    real(8),parameter :: B(1:2)=(/-0.048d0,-0.0269d0/)
    real(8),parameter :: C(1:2)=(/0.002d0,0.0007d0/)
    real(8),parameter :: D(1:2)=(/-0.0116d0,-0.0048d0/)
    real(8) :: onethr,fouthr,onesix,ThrFouPi
    real(8) :: c0,factor
    real(8) :: s0(2),s1(2)
    real(8) :: trho,rs,rssq,rsln,rhoa,rhob
    real(8) :: f,dfdrhoa,dfdrhob
    real(8) :: ecd(2),ecdz,mu(2)
    integer :: i

    ThrFouPi = 3.0d0/( 4.0d0*acos(-1.0d0) )
    onethr   = 1.0d0/3.0d0
    fouthr   = 4.0d0/3.0d0
    onesix   = 1.0d0/6.0d0
    c0       = 1.d0/( 2.0d0**fouthr-2.0d0 )

    factor = 1.0d0
    if ( MSP == 1 ) factor=0.5d0

    Ec = 0.0d0

    do i=ML_0,ML_1

       rhoa = rho(i,1  )*factor
       rhob = rho(i,MSP)*factor
       trho = rhoa + rhob

       if ( trho <= 0.0d0 ) cycle

       rs = (ThrFouPi/trho)**onethr

       if ( rs >= 1.0d0 ) then

          rssq = sqrt(rs)

          ecd(1) = gam(1)/( 1.0d0 + bet1(1)*rssq + bet2(1)*rs )
          ecd(2) = gam(2)/( 1.0d0 + bet1(2)*rssq + bet2(2)*rs )

          f=c0*((2.0d0*rhoa)**fouthr+(2.0d0*rhob)**fouthr-2.0d0*trho**fouthr)
          Ec = Ec + trho*ecd(1) &
               + gam(2)/(  trho**onethr &
                         + bet1(2)*(trho*ThrFouPi)**onesix &
                         + bet2(2)*ThrFouPi**onethr )*f &
               - gam(1)/(  trho**onethr &
                         + bet1(1)*(trho*ThrFouPi)**onesix &
                         + bet2(1)*ThrFouPi**onethr )*f

          f = f/trho**fouthr

          mu(1) = ecd(1) + ecd(1)*onethr*( 0.5d0*bet1(1)/rssq + bet2(1) ) &
                          /( 1.0d0/rs + bet1(1)/rssq + bet2(1) )
          mu(2) = ecd(2) + ecd(2)*onethr*( 0.5d0*bet1(2)/rssq + bet2(2) ) &
                          /( 1.0d0/rs + bet1(2)/rssq + bet2(2) )

       else if ( rs < 1.0d0 ) then

          if ( rs <= 0.0d0 ) stop "calc_ldapz81"

          rsln = log(rs)

          ecd(1) = A(1)*rsln + B(1) + C(1)*rs*rsln + D(1)*rs
          ecd(2) = A(2)*rsln + B(2) + C(2)*rs*rsln + D(2)*rs

          f=c0*((2.0d0*rhoa)**fouthr+(2.0d0*rhob)**fouthr-2.0d0*trho**fouthr) &
               /trho**fouthr

          ecdz = ecd(1) + ( ecd(2) - ecd(1) )*f

          Ec = Ec + trho*ecdz

          mu(1) = ecd(1) - onethr*( A(1) + C(1)*rs*(1.0d0+rsln) + rs*D(1) )
          mu(2) = ecd(2) - onethr*( A(2) + C(2)*rs*(1.0d0+rsln) + rs*D(2) )

       end if

       dfdrhoa = c0*fouthr*2.0d0*rhob &
            *( (2.0d0*rhoa)**onethr - (2.0d0*rhob)**onethr )/trho**fouthr
       dfdrhob =-c0*fouthr*2.0d0*rhoa &
            *( (2.0d0*rhoa)**onethr - (2.0d0*rhob)**onethr )/trho**fouthr

       vco(i,1)   = mu(1) + ( mu(2)-mu(1) )*f + ( ecd(2)-ecd(1) )*dfdrhoa
       vco(i,MSP) = mu(1) + ( mu(2)-mu(1) )*f + ( ecd(2)-ecd(1) )*dfdrhob

    end do ! i

    return
  END SUBROUTINE calc_ldapz81_c


END MODULE xc_ldapz81_module
