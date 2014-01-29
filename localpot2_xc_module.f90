MODULE localpot2_xc_module

  use parallel_module
  use rgrid_module
  use localpot2_variables
  use watch_module

  implicit none

CONTAINS

  SUBROUTINE localpot2_xc(m1_0,m2_0,m3_0,n_in,vout,exc)
    implicit none
    integer,intent(IN) :: m1_0,m2_0,m3_0
    real(8),intent(IN) :: n_in(m1_0,m2_0,m3_0)
    real(8),intent(OUT) :: vout(m1_0,m2_0,m3_0),exc

    real(8),parameter :: a0=0.4581652932831429d0 ,da0=0.119086804055547d0
    real(8),parameter :: a1=2.217058676663745d0  ,da1=0.6157402568883345d0
    real(8),parameter :: a2=0.7405551735357053d0 ,da2=0.1574201515892867d0
    real(8),parameter :: a3=0.01968227878617998d0,da3=0.003532336663397157d0
    real(8),parameter :: b1=1.0d0                ,db1=0.0d0
    real(8),parameter :: b2=4.504130959426697d0  ,db2=0.2673612973836267d0
    real(8),parameter :: b3=1.110667363742916d0  ,db3=0.2052004607777787d0
    real(8),parameter :: b4=0.02359291751427506d0,db4=0.004200005045691381d0

    real(8) :: trho,zeta,d0,d1,d2,d3,d4,d5,fx,rs,epsilon_xc
    real(8) :: a0z,a1z,a2z,a3z,b1z,b2z,b3z,b4z,factor(2)
    real(8) :: de_dr,dr_dn,de_df,df_dz,ct0,ct1,et0,et1
    integer :: i,ierr,i1,i2,i3,m1,m2,m3,i10,i20,i30

    call watch(ct0,et0)

    m1=Ngrid_dense(1)
    m2=Ngrid_dense(2)
    m3=Ngrid_dense(3)

    i10=Igrid_dense(1,1)-1
    i20=Igrid_dense(1,2)-1
    i30=Igrid_dense(1,3)-1

    Exc = 0.d0
    vout = 0.d0

    d0 = 1.0d0
    d1 = 4.d0/3.d0
    d2 = 1.d0/( 2.d0*(2.d0**(1.d0/3.d0)-1.d0) )
    d3 = 3.d0/(4.d0*acos(-1.d0))
    d4 = 1.d0/3.d0
    d5 = 1.d0/(4.d0*acos(-1.d0))

    factor(1) = 2.d0
    factor(2) =-2.d0

    do i3=Igrid_dense(1,3),Igrid_dense(2,3)
    do i2=Igrid_dense(1,2),Igrid_dense(2,2)
    do i1=Igrid_dense(1,1),Igrid_dense(2,1)

       trho = d0*( n_in(i1-i10,i2-i20,i3-i30) )

       if ( trho <= 0.d0 ) cycle

       zeta = 0.0d0

       fx = ( abs(1.d0+zeta)**d1 + abs(1.d0-zeta)**d1 - 2.d0 )*d2

       rs = (d3/trho)**d4

       a0z = a0 + da0*fx
       a1z = a1 + da1*fx
       a2z = a2 + da2*fx
       a3z = a3 + da3*fx

       b1z = b1 + db1*fx
       b2z = b2 + db2*fx
       b3z = b3 + db3*fx
       b4z = b4 + db4*fx

       epsilon_xc = -( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                    /( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )

       Exc = Exc + trho*epsilon_xc

       de_df = ( -( da0 + da1*rs + da2*rs**2 + da3*rs**3 ) &
                 *( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 ) &
                 +( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                 *( db1*rs + db2*rs**2 + db3*rs**3 + db4*rs**4 ) &
               )/( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )**2

       df_dz = d1*d2*( abs(1.d0+zeta)**d4 - abs(1.d0-zeta)**d4 )

       de_dr = ( -( a1z + 2*a2z*rs + 3*a3z*rs**2 ) &
                 *( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 ) &
                 +( a0z + a1z*rs + a2z*rs**2 + a3z*rs**3 ) &
                 *( b1z + 2*b2z*rs + 3*b3z*rs**2 + 4*b4z*rs**3 ) &
               )/( b1z*rs + b2z*rs**2 + b3z*rs**3 + b4z*rs**4 )**2

       dr_dn = -d5/(trho*rs)**2

       vout(i1-i10,i2-i20,i3-i30) = epsilon_xc &
            + trho*( de_dr*dr_dn &
                   + de_df*df_dz*factor(1)*n_in(i1-i10,i2-i20,i3-i30)/trho**2 )

    end do
    end do
    end do

    d0=Exc*dV*Ngrid(0)/(m1*m2*m3)
    call mpi_allreduce(d0,Exc,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    call watch(ct1,et1) ; write(*,*) "localpot2_xc",ct1-ct0,et1-et0

  END SUBROUTINE localpot2_xc

END MODULE localpot2_xc_module
