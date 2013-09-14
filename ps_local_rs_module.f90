MODULE ps_local_rs_module

  use rgrid_module
  use ggrid_module
  use atom_module
  use strfac_module
  use bb_module
  use density_module
  use electron_module
  use pseudopot_module
  use parallel_module
!
  use aa_module
  use esm_rgrid_module
  use esm_rshell_module

  use simc_module

  implicit none

  PRIVATE
  PUBLIC :: init_ps_local_rs,construct_ps_local_rs,construct_ps_density_longloc &
            ,rho_ps,construct_ps_initrho_rs

  real(8),allocatable :: Rcloc(:),vqlg(:,:),vqlgl(:,:),vqls(:,:)
  integer,allocatable :: NRcloc(:)
  real(8),allocatable :: Vion(:)
  real(8),allocatable :: rho_ps(:)

  INTERFACE
     FUNCTION bberf(x)
       real(8) :: bberf,x
     END FUNCTION bberf
  END INTERFACE

CONTAINS


  SUBROUTINE init_ps_local_rs
    implicit none
    integer :: i,ig,ik,iorb,MMr,NRc,MKI
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell
    real(8),allocatable :: vshort(:),tmp(:)

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.d0)
    const = 4.d0*Pi/Vcell

    MMr=maxval(Mr)
    allocate( vqls(MMr,MKI)   ) ; vqls=0.d0
    allocate( vqlg(NMGL,MKI)  ) ; vqlg=0.d0
    allocate( vqlgl(NMGL,MKI) ) ; vqlgl=0.d0
    allocate( NRcloc(MKI)     ) ; NRcloc=0
    allocate( Rcloc(MKI)      ) ; Rcloc=0.d0

    allocate( vshort(MMr) )

    do ik=1,MKI

       MMr = Mr(ik)

       Rc=0.d0
       NRc=0
       do iorb=1,norb(ik)
          Rc=max( Rc, Rps(iorb,ik) )
          NRc=max( NRc, NRps(iorb,ik) )
       end do

       if ( Rc<1.d-8 ) Rc=5.d0
       if ( NRc==0 ) then
          do i=1,MMr
             if ( rad(i,ik)>Rc ) then
                NRc=i
                exit
             end if
          end do
       end if

       call simc(rad(1,ik),vql(1,ik),Rc,Zps(ik),parloc(1,ik),MMr)

       p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
       p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))

       do i=1,MMr
          r=rad(i,ik)
          if ( r<1.d-9 ) then
             vlong=-2.d0*Zps(ik)/sqrt(Pi)*(p1*p2+p3*p4)
          else
             vlong=-Zps(ik)/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
          end if
          vshort(i)=vql(i,ik)-vlong
          vqls(i,ik)=vql(i,ik)-vlong
          if( NRcloc(ik)==0 )then
             if( abs(vshort(i))<1.d-8 ) then
                NRcloc(ik)=i
                Rcloc(ik)=r
             end if
          end if
       end do

       if ( NRcloc(ik)==0 ) then
          Rcloc(ik)=Rc
          NRcloc(ik)=NRc
       end if

       allocate( tmp(MMr) )

       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,MMr
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*rab(i,ik)
             end do
          else
             do i=1,MMr
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                       +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          vqlg(ig,ik)=sum0*const
       end do

       p1=-Zps(ik)*parloc(1,ik) ; p2=0.25d0/parloc(2,ik)
       p3=-Zps(ik)*parloc(3,ik) ; p4=0.25d0/parloc(4,ik)
       do ig=1,NMGL
          G2=GG(ig)
          if ( G2 == 0.d0 ) then
             vqlgl(ig,ik)=-(p1*p2+p3*p4)*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
          else
             vqlgl(ig,ik)=(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
          end if
!          vqlg(ig,ik)=Vloc(sqrt(G2),1)/Vcell
       end do
       deallocate( tmp )

    end do ! ik

    deallocate( vshort )

  END SUBROUTINE init_ps_local_rs

  SUBROUTINE simp(f,s,m)
    implicit none
    integer,intent(IN)  :: m
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,n,nn,nmax
    n=size(f) ; nmax=int(n/m)*m
    do i=0,m
       nmax=nmax+i ; if ( nmax>=n ) exit
    end do
    allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
    select case(m)
    case default
       s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
    case(2)
       s=0.d0
       do i=1,nmax-2,2
          s = s + g(i) + 4.d0*g(i+1) + g(i+2)
       end do
       s=s/3.d0
    case(4)
       s=0.d0
       do i=1,nmax-4,4
          s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
       end do
       s=s*2.d0/45.d0
    case(6)
       s=0.d0
       do i=1,nmax-6,6
          s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3) &
               +27*g(i+4)+216*g(i+5)+41*g(i+6)
       end do
       s=s/140.d0
    end select
    deallocate( g )
    return
  END SUBROUTINE simp


  SUBROUTINE construct_ps_local_rs(vout)
    use aa_module
    use esm_rgrid_module
    use esm_rshell_module
    implicit none
    real(8),intent(OUT) :: vout(*)
    integer :: i,a,i1,i2,i3,j1,j2,j3,j,ik,ir,mm,m1,m2,MMr,NML
    real(8) :: H1,H2,H3,x1,x2,x3,x,y,z,r,chk_sum,tmp_sum,tmp_sum0
    real(8) :: p1,p2,p3,p4,err,err0,v,v0,xion,yion,zion
    real(8) :: c1,c2,c3,d1,d2,d3
    integer,allocatable :: map_rgrd(:,:),map2_rgrd(:,:)
    integer :: imember,ishell,ierr

    NML = ML1_ESM - ML0_ESM + 1

    vout(1:NML) = 0.d0

    H1 = Hgrid(1)
    H2 = Hgrid(2)
    H3 = Hgrid(3)

    c1 = 1.d0/Ngrid(1)
    c2 = 1.d0/Ngrid(2)
    c3 = 1.d0/Ngrid(3)

    call construct_lattice_rshell

!    allocate( map_rgrd(3,n1:n2) ) ; map_rgrd=0
!    allocate( map2_rgrd(3,m1:m2) ) ; map2_rgrd=0
!    call construct_esm_rgrid(ML_ESM,map_rgrd,MK_ESM,map2_rgrd)

    do a=1,Natom

       ik=ki_atom(a)
       p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
       p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))
       MMr=Mr(ik)

       xion=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
       yion=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
       zion=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)

       chk_sum = 0.d0

       do ishell=1,Mshell
       do imember=1,num_member_shell(ishell)

          j1 = lattice_shell(1,imember,ishell)
          j2 = lattice_shell(2,imember,ishell)
          j3 = lattice_shell(3,imember,ishell)

          x1 = xion + j1*aa(1,1)+j2*aa(1,2)+j3*aa(1,3)
          x2 = yion + j1*aa(2,1)+j2*aa(2,2)+j3*aa(2,3)
          x3 = zion + j1*aa(3,1)+j2*aa(3,2)+j3*aa(3,3)

          do i=ML0_ESM,ML1_ESM

             d1 = c1*LL(1,i)
             d2 = c2*LL(2,i)
             d3 = c3*LL(3,i)
             x  = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3 - x1
             y  = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3 - x2
             z  = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3 - x3
             r  = sqrt(x*x+y*y+z*z)

             do ir=1,MMr
                if ( r < rad(ir,ik) ) exit
             end do
             if ( ir > MMr ) then
                v0=0.d0
                cycle
             end if
             err0=1.d10
             do mm=1,20
                m1=max(1,ir-mm)
                m2=min(ir+mm,MMr)
                call polint(rad(m1,ik),vqls(m1,ik),m2-m1+1,r,v,err)
                if ( abs(err) < err0 ) then
                   v0=v
                   err0=abs(err)
                   if ( err0 < 1.d-8 ) exit
                end if
             end do

             j = i - ML0_ESM + 1
             vout(j) = vout(j) + v0

          end do ! i

       end do ! imember

          tmp_sum0=sum( abs(vout(1:NML))**2 )
          call mpi_allreduce(tmp_sum0,tmp_sum,1,mpi_real8,mpi_sum,comm_grid,ierr)
          if ( myrank == 0 ) write(*,'(1x,"ishell",i8,2g20.10)') ishell,tmp_sum,chk_sum
          if ( abs(tmp_sum-chk_sum) < 1.d-10 ) then
             exit
          else
             chk_sum=tmp_sum
          end if

       end do ! ishell

    end do ! a

  END SUBROUTINE construct_ps_local_rs


  SUBROUTINE construct_ps_density_longloc
    implicit none
    integer :: ik,j,ierr,ML_0,ML_1
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3,pi,q1,q2,q3,q4,c4
    real(8) :: x,y,z,rr,rrmin,chk_conv,const,x2
    real(8),allocatable :: p1(:),p2(:),p3(:),p4(:)
    integer :: a,i,i1,i2,i3,j1,j2,j3,k1,k2,k3,ishell,imember
    integer,allocatable :: map_grid(:,:),map2_grid(:,:)

    pi = acos(-1.d0)

    allocate( rho_ps(ML0_ESM:ML1_ESM) ) ; rho_ps=0.d0

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)
    const = 1.d0/pi**1.5d0

    allocate( p1(Nelement),p2(Nelement) )
    allocate( p3(Nelement),p4(Nelement) )

    do ik=1,Nelement
       p2(ik)= parloc(2,ik)
       p1(ik)=-Zps(ik)*parloc(1,ik)*const*p2(ik)**1.5d0
       p4(ik)= parloc(4,ik)
       p3(ik)=-Zps(ik)*parloc(3,ik)*const*p4(ik)**1.5d0
    end do

!    allocate( map_grid(3,ML_ESM)  ) ; map_grid=0
!    allocate( map2_grid(3,MK_ESM) ) ; map2_grid=0
!    call construct_esm_rgrid(ML_ESM,map_grid,MK_ESM,map2_grid)

    rho_ps(:) = 0.d0
    chk_conv  = 0.d0

    do ishell=1,Mshell
    do imember=1,num_member_shell(ishell)

       j1 = lattice_shell(1,imember,ishell)
       j2 = lattice_shell(2,imember,ishell)
       j3 = lattice_shell(3,imember,ishell)

       do a=1,Natom

          ik=ki_atom(a)
          a1=aa_atom(1,a) + j1
          a2=aa_atom(2,a) + j2
          a3=aa_atom(3,a) + j3
          q1=p1(ik)
          q2=p2(ik)
          q3=p3(ik)
          q4=p4(ik)

          do i=ML0_ESM,ML1_ESM
             d1=LL(1,i)*c1-a1
             d2=LL(2,i)*c2-a2
             d3=LL(3,i)*c3-a3
             x =aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3
             y =aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3
             z =aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3
             rr=x*x+y*y+z*z
             rho_ps(i) = rho_ps(i) + q1*exp(-q2*rr)+q3*exp(-q4*rr)
          end do

       end do ! a

    end do ! imember

       x=sum(rho_ps)*dV
       call mpi_allreduce(x,y,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       if ( myrank == 0 ) then
          write(*,'(1x,i5,"sum(rho_ps)=",4g20.8)') ishell, &
               y,minval(rho_ps),maxval(rho_ps),y-chk_conv
       end if
       if ( abs(y-chk_conv) < 1.d-12 ) exit
       chk_conv=y

    end do ! ishell

    deallocate( p4,p3,p2,p1 )

  END SUBROUTINE construct_ps_density_longloc


  SUBROUTINE construct_ps_initrho_rs(n1,n2,nsp,initrho)
    implicit none
    integer,intent(IN) :: n1,n2,nsp
    real(8),intent(OUT) :: initrho(n1:n2,nsp)
    integer :: MMr,ik,ishell,imember,j1,j2,j3,a,ir
    integer :: i,mm,m1,m2,ierr
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3,x,y,z,r,v0,err0
    real(8) :: v,err,chk_conv,pi4
    real(8),allocatable :: cdt(:,:)

    initrho=0.d0

    pi4 = 4.d0*acos(-1.d0)

    if ( myrank == 0 ) then
       do ik=1,Nelement
          write(*,*) "sum(cdd*rab)",sum(cdd(:,ik)*rab(:,ik)),ik
       end do
    end if

    mm=size(cdd,1)
    allocate( cdt(mm,Nelement) ) ; cdt=0.d0

    do ik=1,Nelement
       do ir=2,Mr(ik)
          cdt(ir,ik) = cdd(ir,ik)/pi4/rad(ir,ik)**2
       end do
       cdt(1,ik) = (cdt(3,ik)-cdt(2,ik))/(rad(3,ik)-rad(2,ik)) &
                 * ( rad(1,ik) - rad(2,ik) ) + cdt(2,ik)
    end do
!    if ( myrank == 0 ) then
!       rewind 10
!       do ir=1,minval(Mr)
!          write(10,*) ir,(rad(ir,ik),cdt(ir,ik),ik=1,Nelement)
!       end do
!    end if

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    chk_conv=0.d0

    do ishell=1,Mshell
    do imember=1,num_member_shell(ishell)

       j1 = lattice_shell(1,imember,ishell)
       j2 = lattice_shell(2,imember,ishell)
       j3 = lattice_shell(3,imember,ishell)

       do a=1,Natom

          ik = ki_atom(a)
          MMr= Mr(ik)
          a1 = aa_atom(1,a) + j1
          a2 = aa_atom(2,a) + j2
          a3 = aa_atom(3,a) + j3

          do i=ML0_ESM,ML1_ESM
             d1=LL(1,i)*c1-a1
             d2=LL(2,i)*c2-a2
             d3=LL(3,i)*c3-a3
             x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3
             y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3
             z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3
             r = sqrt(x*x+y*y+z*z)
             do ir=1,MMr
                if ( r < rad(ir,ik) ) exit
             end do
             if ( ir > MMr ) then
                v0=0.d0
                cycle
             end if
             err0=1.d10
             do mm=1,20
                m1=max(1,ir-mm)
                m2=min(ir+mm,MMr)
                call polint(rad(m1,ik),cdt(m1,ik),m2-m1+1,r,v,err)
                if ( abs(err) < err0 ) then
                   v0=v
                   err0=abs(err)
                   if ( err0 < 1.d-8 ) exit
                end if
             end do
             initrho(i,1) = initrho(i,1) + v0
          end do

       end do ! a

    end do ! imember

       x=sum(initrho)*dV
       call mpi_allreduce(x,y,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       if ( myrank == 0 ) then
          write(*,'(1x,i5," sum(initrho)=",4g20.8)') ishell, &
               y,minval(initrho),maxval(initrho),y-chk_conv
       end if
       if ( abs(y-chk_conv) < 1.d-12 ) exit
       chk_conv=y

    end do ! ishell

    deallocate( cdt )

  END SUBROUTINE construct_ps_initrho_rs

END MODULE ps_local_rs_module
