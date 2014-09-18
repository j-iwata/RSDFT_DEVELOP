MODULE localpot2_xc_module

  use parallel_module
  use rgrid_module
  use localpot2_variables, only: Igrid_dense, Ngrid_dense
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_xc

  integer :: MSP=1

CONTAINS


  SUBROUTINE localpot2_xc( n_in, vout, exc )
    implicit none
    real(8),intent(IN) :: n_in(:,:,:)
    real(8),intent(OUT) :: vout(:,:,:), exc
    real(8) :: d0
    integer :: ierr,m1,m2,m3

    exc  = 0.0d0
    vout = 0.0d0

    m1=Ngrid_dense(1)
    m2=Ngrid_dense(2)
    m3=Ngrid_dense(3)

    call calc_ldapz81_x( n_in, vout, exc )
    call calc_ldapz81_c( n_in, vout, exc )

    d0=Exc*dV*Ngrid(0)/(m1*m2*m3)
    call mpi_allreduce(d0,Exc,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

  END SUBROUTINE localpot2_xc


  SUBROUTINE calc_ldapz81_x( rho, vex, Ex )
    implicit none
    real(8),intent(IN)  :: rho(:,:,:)
    real(8),intent(INOUT) :: vex(:,:,:), Ex
    real(8) :: onetwo,onethr,thrfou,fouthr,thrPi
    real(8) :: cnst,exd,trho
    integer :: i1,i2,i3,s

    onethr = 1.0d0/3.0d0
    thrfou = 3.0d0/4.0d0
    fouthr = 4.0d0/3.0d0
    thrPi  = 3.0d0/acos(-1.0d0)

    cnst = 1.0d0
    if ( MSP == 2 ) cnst = (2.0d0)**onethr

    do s=1,MSP

       do i3=1,Igrid_dense(2,3)-Igrid_dense(1,3)+1
       do i2=1,Igrid_dense(2,2)-Igrid_dense(1,2)+1
       do i1=1,Igrid_dense(2,1)-Igrid_dense(1,1)+1

          trho = rho(i1,i2,i3)

          exd = -cnst*thrfou*( thrPi*trho )**onethr
          Ex = Ex + trho*exd
          vex(i1,i2,i3) = vex(i1,i2,i3) + fouthr*exd

       end do
       end do
       end do

    end do

  END SUBROUTINE calc_ldapz81_x


  SUBROUTINE calc_ldapz81_c( rho, vco, Ec )
    implicit none

    real(8),intent(IN)  :: rho(:,:,:)
    real(8),intent(INOUT) :: vco(:,:,:),Ec

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
    integer :: i1,i2,i3

    ThrFouPi = 3.0d0/( 4.0d0*acos(-1.0d0) )
    onethr   = 1.0d0/3.0d0
    fouthr   = 4.0d0/3.0d0
    onesix   = 1.0d0/6.0d0
    c0       = 1.d0/( 2.0d0**fouthr-2.0d0 )

    factor = 1.0d0
    if ( MSP == 1 ) factor=0.5d0

    do i3=1,Igrid_dense(2,3)-Igrid_dense(1,3)+1
    do i2=1,Igrid_dense(2,2)-Igrid_dense(1,2)+1
    do i1=1,Igrid_dense(2,1)-Igrid_dense(1,1)+1

       rhoa = rho(i1,i2,i3)*factor ; rhoa=abs(rhoa)
       rhob = rho(i1,i2,i3)*factor ; rhob=abs(rhob)
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

       vco(i1,i2,i3) = vco(i1,i2,i3) + &
            mu(1) + ( mu(2)-mu(1) )*f + ( ecd(2)-ecd(1) )*dfdrhoa
!       vco(i1,i2,i3) = mu(1) + ( mu(2)-mu(1) )*f + ( ecd(2)-ecd(1) )*dfdrhob

    end do ! i1
    end do ! i2
    end do ! i3

    return
  END SUBROUTINE calc_ldapz81_c


END MODULE localpot2_xc_module
